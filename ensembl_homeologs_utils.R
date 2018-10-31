library("biomaRt")
library("tidyr")
library("dplyr")
library("tibble")


#Find the homeologs of a given gene.
#Returns a vector with the homeologs of the input gene (including itself excep if it's a singleton)
find_homeologs <- function(gene_id){
  print(gene_id)
  plants_ensembl <- useMart("plants_mart", dataset = "taestivum_eg_gene" ,host="http://plants.ensembl.org")
  homeologs_vector <- getBM(attributes=c("taestivum_eg_homoeolog_gene"), filters = "ensembl_gene_id", values = gene_id, plants_ensembl) %>%
    pull(taestivum_eg_homoeolog_gene)
  #Add the homeologs of the homeologs to the list
  homeologs_vector <- c(homeologs_vector,lapply(homeologs_vector,
                                                function(gene) getBM(attributes=c("taestivum_eg_homoeolog_gene"),
                                                                     filters = "ensembl_gene_id",
                                                                     values = gene,
                                                                     plants_ensembl) %>%
                                                  pull(taestivum_eg_homoeolog_gene)) %>% unlist())
  #Remove duplicates
  homeologs_vector <- homeologs_vector %>% unique()
  return(homeologs_vector)
}

#Check if ALL the homeologs of a given gene appear in a list of gene ids (e.g. DEGs).
#Returns TRUE or FALSE
all_homeologs_appear <- function(gene_id, gene_vector, consider_singletons=F){
  gene_homeologs <- find_homeologs(gene_id)
  if(consider_singletons){
    return(length(gene_homeologs) == lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length())
  }else{ #Don't consider genes which do not have any homeologs, aka singletons
    if(length(gene_homeologs)>0){
      return(all(gene_homeologs%in%gene_vector))
    }else{
      return(FALSE)
    }
  }
}

#Check if ANY of the homeologs of a given gene appear in a list of gene ids (e.g. DEGs).
#Returns TRUE or FALSE
any_homeologs_appear <- function(gene_id, gene_vector, consider_singletons=F){
  gene_homeologs <- find_homeologs(gene_id)
  if(consider_singletons){
    if(length(gene_homeologs)==0){
      return(TRUE)
    }else{
      return(lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length() > 1)
    }
  }else{ #Don't consider genes which do not have any homeologs, aka singletons
    return(lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length() > 1)
  }
}

#Check if a given gene has any homeologs or is a singleton
#Returns TRUE if it is a singleton and FALSE otherwise
is_singleton <- function(gene_id){
  if(find_homeologs(gene_id) %>% length() == 0) return(T) else return(F)
}

#Given a gene extract the baseMean, log2FoldChange and padj of itself and the homeologs
#Returns a data frame with those values
retrieve_gene_results <- function(gene_vector, DESeq_results, merge_homeologs=F){
  gene_in_results <- gene_vector %in% rownames(DESeq_results)
  
  #Create a data frame of NAs for those homeologs that do not appear in the results
  #If all the homeologs appear in the results this will be empty and it will not affect the results
  NA_results <- data.frame(matrix(nrow=length(gene_vector[!gene_in_results]), ncol = 3))
  rownames(NA_results) <- gene_vector[!gene_in_results]
  colnames(NA_results) <- c("baseMean", "log2FoldChange", "padj")  
  gene_results <- rbind(as.data.frame(DESeq_results)[gene_vector[gene_in_results], c("baseMean", "log2FoldChange", "padj")], NA_results)
  if(merge_homeologs){ #Combine all homeologs in one row
    homeologs_row <- gene_results %>%
      rownames_to_column() %>% 
      #arrange(padj) %>% #Order homeologs by p value
      arrange(rowname) %>% 
      format.data.frame() %>% #Make sure there won't be problems with decimals and scientific numbers
      mutate(group=1) %>% #Dummy group id to put the rows together
      mutate(merged_cols=paste(rowname, baseMean, log2FoldChange, padj, sep=";")) %>%
      group_by(group) %>%
      summarise(merged_cols_and_rows=paste0(merged_cols, collapse = ";")) %>%
      dplyr::select(merged_cols_and_rows) %>%
      as.character() 
    return(paste(length(gene_vector), homeologs_row, sep=";"))
  }else{
    return(gene_results %>% as.data.frame())
  }
}

homeolog_results_as_data_frame <- function(gene_vector, DESeq_results){
  homeolog_results_list <- lapply(gene_vector,
                                  function(gene) {
                                    gene_homeologs <- find_homeologs(gene)
                                    retrieve_gene_results(gene_homeologs, DESeq_results, T)}) %>%
                            unlist() %>%
                            strsplit(split=";")
  n.obs <- sapply(homeolog_results_list, length)
  seq.max <- seq_len(max(n.obs))
  homeolog_results_df <- t(sapply(homeolog_results_list, "[", i=seq.max)) %>% replace_na("") %>% as.data.frame()
  homeolog_results_df <- homeolog_results_df %>% distinct()
  colnames(homeolog_results_df) <- c("total_homeologs",rep(c("Gene id", "baseMean", "log2FoldChange", "padj"), max(n.obs)/4))
  #print(homeolog_results_df)
  return(homeolog_results_df)
}

genes_same_FC_sign <- function(gene_vector, DESeq_results){
  genes_FC <- retrieve_gene_results(gene_vector, DESeq_results, F) %>% dplyr::select(log2FoldChange)
  if(anyNA(genes_FC)){
    return(FALSE)
  }else{
    return(genes_FC %>% sign() %>% colMeans %>% abs() == 1)
  }
}