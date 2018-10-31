homeologs <- read.table("homoe_W5-10degW1-2nodeg_up.tab", sep="\t", comment="", as.is=T, header = T)
homeologs %>% group_by(Gene_stable_ID) %>% summarise(homoeologs=paste(Triticum_aestivum_homoeologue_gene_stable_ID, collapse = ", ")) %>% write.xlsx("HIvLO_N_W10-5W1-2_comp_homoeo.xlsx", sheetName = "homeo W1-2degW5-10nodeg down", append = T)
