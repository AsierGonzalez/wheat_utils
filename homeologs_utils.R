library("tidyr")
library("dplyr")

#Parse the file with the IWGSC homeologs to create a vector of homeologs
parse_iwgsc_homeologs <- function(){
  
  iwgsc_homeologs <- read.table(file.path("//salt","wheat_rnaseq","NRgeneV1","Triticum_aestivum_V1_PGSB.homeologous_gene_groups.txt"), sep="\t", header=T, comment="", as.is=T, quote = "\"")
  homeologs <- iwgsc_homeologs %>%
                unite(Homeologs, A, B, D, sep = ",") %>% #Merge the columns A, B and D homeologs
                select(Homeologs) %>%
                distinct() %>% #Remove duplicates
                apply(1, function(x) gsub("[0-9A-Za-z]+LC", "",x)) %>% #Delete LC genes
                grep("TraesCS", ., value = T) %>% #Select rows that contain genes
                gsub(",+",",",.) %>% #Remove consecutive commas caused by LC genes            
                gsub(",$","",.) %>% #Remove trailing commas
                gsub("^,","",.) %>% #Remove leading commas
                gsub("01G", "02G", .) #Change gene ids so that they match those of the DEGs
    
  return(homeologs)
}

#Find the homeologs of a given gene
#Returns a vector with the homeologs of the input gene (including itself)
find_homeologs <- function(gene_id, homeologs){
  homeologs_vector <- grep(gene_id, homeologs, value = T) %>% strsplit(",") %>% unlist()
  if(is.null(homeologs_vector)){
    homeologs_vector <- gene_id
  }
  return(homeologs_vector)
}

#Check if ALL the homeologs of a given gene appear in a list of gene ids (e.g. DEGs)
#Returns TRUE or FALSE
all_homeologs_appear <- function(gene_id, homeologs, gene_vector, consider_singletons=F){
  gene_homeologs <- find_homeologs(gene_id, homeologs)
  if(consider_singletons){
    return(length(gene_homeologs) == lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length())
  }else{ #Don't consider genes which do not have any homeologs, aka singletons
    if(length(gene_homeologs)>1){
      return(all(gene_homeologs%in%gene_vector))
    }else{
      return(FALSE)
    }
  }
}

#Check if ANY of the homeologs of a given gene appear in a list of gene ids (e.g. DEGs)
#Returns TRUE or FALSE
any_homeologs_appear <- function(gene_id, homeologs, gene_vector, consider_singletons=F){
  gene_homeologs <- find_homeologs(gene_id, homeologs)
  if(consider_singletons){
    if(length(gene_homeologs)==1){
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
is_singleton <- function(gene_id, homeologs){
  if(find_homeologs(gene_id, homeologs) %>% length() == 1) return(T) else return(F)
}

#Given a gene extract the baseMean, log2FoldChange and padj of itself and the homeologs
#Returns a data frame with those values
retrieve_homeolog_results <- function(gene_id, homeologs, DESeq_results, merge_homeologs=F){
  gene_homeologs <- find_homeologs(gene_id, homeologs)
  gene_in_results <- gene_homeologs %in% rownames(DESeq_results)
  
  #Create a data frame of NAs for those homeologs that do not appear in the results
  #If all the homeologs appear in the results this will be empty and it will not affect the results
  NA_results <- data.frame(matrix(nrow=length(gene_homeologs[!gene_in_results]), ncol = 3))
  rownames(NA_results) <- gene_homeologs[!gene_in_results]
  colnames(NA_results) <- c("baseMean", "log2FoldChange", "padj")  
  homeolog_results <- rbind(as.data.frame(DESeq_results)[gene_homeologs[gene_in_results], c("baseMean", "log2FoldChange", "padj")], NA_results)
  if(merge_homeologs){
    homeolog_results %>%
      rownames_to_column() %>% 
      arrange(padj) %>% 
      format.data.frame() %>% #Make sure there won't be problems with decimals and scientific numbers
      mutate(group=1) %>% #Dummy group id to put the rows together
      mutate(merged_cols=paste(rowname, baseMean, log2FoldChange, padj, sep=";")) %>%
      group_by(group) %>%
      summarise(merged_cols_and_rows=paste0(merged_cols, collapse = ";")) %>%
      select(merged_cols_and_rows) %>%
      as.character() 
  }else{
    return(homeolog_results)
  }
}

homeolog_results_as_data_frame <- function(gene_vector, homeologs, DESeq_results){
  homeolog_results_list <- lapply(gene_vector,
                                  function(gene) retrieve_homeolog_results(gene, homeologs, DESeq_results, T)) %>%
                                  unlist() %>%
                                  strsplit(split=";")
  n.obs <- sapply(homeolog_results_list, length)
  seq.max <- seq_len(max(n.obs))
  homeolog_results_df <- t(sapply(homeolog_results_list, "[", i=seq.max)) %>% replace_na("") %>% as.data.frame()
  homeolog_results_df <- homeolog_results_df %>% distinct()
  colnames(homeolog_results_df) <- rep(c("Gene id", "baseMean", "log2FoldChange", "padj"), max(n.obs)/4)
  return(homeolog_results_df)
}

all_homeologs_same_FC_sign <- function(gene_id, homeologs, DESeq_results){
  homeolog_FC <- retrieve_homeolog_results(gene_id, homeologs, DESeq_results, F) %>% select(log2FoldChange)
  if(anyNA(homeolog_FC)){
    return(FALSE)
  }else{
    return(homeolog_FC %>% sign() %>% colMeans %>% abs() == 1)
  }
}