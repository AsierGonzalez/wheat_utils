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
find_homeologs <- function(gene_id, homeologs){
  homeologs_vector <- grep(gene_id, homeologs, value = T) %>% strsplit(",") %>% unlist()
  if(is.null(homeologs_vector)){
    homeologs_vector <- gene_id
  }
  return(homeologs_vector)
}

#Check if ALL the homeologs of a given gene appear in a list of gene ids (e.g. DEGs)
#Returns TRUE or FALSE
all_homeologs_appear <- function(gene_id, homeologs, gene_vector, consider_singletons){
  gene_homeologs <- find_homeologs(gene_id, homeologs)
  if(consider_singletons){
    return(length(gene_homeologs) == lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length())
  }else{ #Don't consider genes which do not have any homeologs, aka singletons
    if(length(gene_homeologs)>1){
      return(length(gene_homeologs) == lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length())
    }else{
      return(FALSE)
    }
  }
}

#Check if ANY of the homeologs of a given gene appear in a list of gene ids (e.g. DEGs)
#Returns TRUE or FALSE
any_homeologs_appear <- function(gene_id, homeologs, gene_vector, consider_singletons){
  gene_homeologs <- find_homeologs(gene_id, homeologs)
  if(consider_singletons){
    if(length(gene_homeologs<=1)){
      return(TRUE)
    }else{
      return(lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length() > 1)
    }
  }else{ #Don't consider genes which do not have any homeologs, aka singletons
    return(lapply(gene_homeologs, function(gene) grep(gene, gene_vector)) %>% unlist() %>% length() > 1)
  }
}

#Given a gene