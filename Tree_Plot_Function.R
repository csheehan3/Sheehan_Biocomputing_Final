Tree_Plot_Function <- function(the_dataframe, the_gene, the_GO_term){
if(!exists("GO_metabolic_list", mode = "list")){
  print("Warning: need to load GO list for this function")
  load("GO_metabolic_list.Rdata")
}
  BACH1_row <- filter(the_dataframe, the_dataframe$entrez_dictionary==the_gene) %>% 
    dplyr::select(matches("TCGA")) %>% 
    as.numeric()
  filtered_frame <- filter(the_dataframe, the_dataframe$entrez_dictionary %in% GO_metabolic_list[[the_GO_term]]) ##filters frame to only genes part of this GO term
  GO_single_spearman_values <- c()
  GO_single_p_values <- c()
  for (n in 1:nrow(filtered_frame)){ ##runs correlations for each gene in set against the gene-of-interest
    single_gene_row <- filtered_frame[n,] %>% 
      dplyr::select(matches("TCGA")) %>% 
      as.numeric()
    correlation_value <- cor.test(BACH1_row, single_gene_row, method=c("spearman"))
    GO_single_spearman_values <- c(GO_single_spearman_values, correlation_value$estimate)
    GO_single_p_values <- c(GO_single_p_values, correlation_value$p.value)
  }
  hgnc_conversion <- getBM(attributes = c('entrezgene_id', ###get all the gene symbols for the genes of the GO set for plotting
                                          'hgnc_symbol'),
                           filters = 'entrezgene_id', 
                           values = filtered_frame$entrez_dictionary,
                           mart = mart)
  OXPHOS_tibble <- tibble("Gene_Name"=hgnc_conversion$hgnc_symbol, "Spearman_Coefficient"=GO_single_spearman_values, "P-Value"=GO_single_p_values)
  OXPHOS_tibble <- arrange(OXPHOS_tibble, OXPHOS_tibble$`P-Value`)
  ####Generate the Tree Plots for each of the interesting terms
  ggplot(OXPHOS_tibble, aes(x=reorder(OXPHOS_tibble$Gene_Name, OXPHOS_tibble$`P-Value`), y=OXPHOS_tibble$Spearman_Coefficient)) + 
    geom_bar(stat="identity" , width=0.3, col="red") +
    theme_light() +
    theme(axis.text=element_text(size=3)) +
    labs(title="BACH1 correlations with OXPHOS Gene Set", y="Spearman Coefficient", x="P-Value")
}
