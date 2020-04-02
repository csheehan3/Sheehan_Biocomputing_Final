fgsea_function <- function(the_dataframe, the_gene){
  global_mean_gene_expression <- filter(the_dataframe, the_dataframe$entrez_dictionary==the_gene) %>% 
    dplyr::select(matches("TCGA")) %>% ###get only the columns for tumor samples
    as.numeric() %>% 
    mean()
  gene_row <- filter(the_dataframe, the_dataframe$entrez_dictionary==the_gene) %>% dplyr::select(matches("TCGA"))
  high_gene_samples <- the_dataframe[1:nrow(the_dataframe), which(gene_row > global_mean_gene_expression)] %>% 
    dplyr::select(matches("TCGA"))
  low_gene_samples <- the_dataframe[1:nrow(the_dataframe), which(gene_row < global_mean_gene_expression)] %>% 
    dplyr::select(matches("TCGA"))
  average_high_list <- c(apply(high_gene_samples,1,mean))
  average_low_list <- c(apply(low_gene_samples,1,mean))
  difference_in_expression <- average_high_list - average_low_list #produces vector containing difference in mean expression for every gene of high and low gene samples
  entrez_expression_dataframe <- data.frame(entrez_dictionary, difference_in_expression) %>% #organizes everything into a dataframe, removes NA genes, and ranked by difference
    filter(., is.na(.$entrez_dictionary)==FALSE) %>% 
    .[order(-.$difference_in_expression),]
  #######Proper format for the fgsea package uses a ranked, named vector with the Entrez ID and difference in expression
  ranks_vector <- entrez_expression_dataframe$difference_in_expression
  names(ranks_vector) <- as.character(entrez_expression_dataframe$entrez_dictionary)
  assign("ranks_vector", ranks_vector, envir=.GlobalEnv) #This is used later so need to save to global environment
  library(fgsea)
  pathways <- reactomePathways(as.character(entrez_expression_dataframe$entrez_dictionary)) #grabs potential pathways to explore with gsea based on the genes being submitted
  fgseaRes <- fgsea(pathways, ranks_vector, nperm=1000, maxSize=500)
  #######Here we produce can produce a number of graphs from the output of the GSEA analysis.
  #######figure that displays the 15 most significant lists for positive and negative enrichments
  #######with the gene high tumor samples. Next are two enrichments plots for two different lists of interest. 
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks_vector, fgseaRes, 
                gseaParam = 0.5)
}