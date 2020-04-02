#####  Here I make a Tibble with the GO terms I'm interested in & grab all the Entrez gene names included in each set
GO_cellular_metabolic_process <- read_tsv(file="basket.tsv", col_names = FALSE)
additional_terms <- read_tsv(file="basket (1).tsv", col_names=FALSE)
GO_cellular_metabolic_process <- bind_rows(GO_cellular_metabolic_process, additional_terms)
GO_query_list <- c()
for (GO_term in 1:nrow(GO_cellular_metabolic_process)){ ###collect the GO terms into a list to send to biomaRt
  GO_query_list <- c(GO_query_list, as.character(GO_cellular_metabolic_process[GO_term,1]))
}

GO_metabolic_list <- c() ##empty vector that will collect all the GO terms & 
n=1 ## n represents the index for things put in the GO_metabolic_list vector
for (term in 1:length(GO_query_list)){
  GO_conversion <- getBM(attributes = c('hgnc_symbol',
                                        'entrezgene_id'),
                         filters = "go_parent_term", 
                         values = GO_query_list[term],
                         mart = mart)
  print(paste("completed:", GO_query_list[term])) ###just to let you know how its running (will be very slow)
  if (length(GO_conversion$entrezgene_id) > 50 && length(GO_conversion$entrezgene_id) < 300){ ##I want sets containing 50 to 300 genes
    GO_metabolic_list[[n]] <-  as.list(GO_conversion$entrezgene_id)
    names(GO_metabolic_list)[[n]] <- GO_query_list[term]
    n=n+1 
  }
}
assign("GO_metabolic_list", GO_metabolic_list, envir=.GlobalEnv)
save(GO_metabolic_list, file="GO_metabolic_list.Rdata")
