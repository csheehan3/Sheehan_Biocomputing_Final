simulation_function <- function(the_dataframe, set_size, total_simulations){
  set.seed(10)
  ##
  mean_spearman_values <- c()
  simulation_count=0
  while (simulation_count<total_simulations) {
    A_very_interesting_gene_index <- sample(1:nrow(the_dataframe), 1, replace=TRUE)
    A_very_interesting_gene <- as.numeric(the_dataframe[A_very_interesting_gene_index, 3:ncol(the_dataframe)])
    random_row_indices <- sample(1:nrow(the_dataframe), 200, replace=TRUE)
    random_rows <- the_dataframe[random_row_indices, ] %>% dplyr::select(matches("TCGA"))
    spearman_list <- c()
    for (A_gene in 1:set_size){
      A_correlation <- cor(as.numeric(random_rows[A_gene,]), A_very_interesting_gene, method = c("spearman"))
      spearman_list <- c(spearman_list, A_correlation) ##throw the correlation coefficient in the list for the gene set
      if(is.na(A_correlation)==TRUE){ ##This will skip a simulation, if value is NA for corr. coefficient for any of the genes
        print("warning! skipping simulation")
        print(head(as.numeric(random_rows[A_gene,])))
        print(head(A_very_interesting_gene))
        break()
      }
    }
    mean_spearman_values <- c(mean_spearman_values, mean(spearman_list)) #If no NAs present, update the list to include last simulation
    simulation_count <- length(mean_spearman_values)
    print(paste("simulation count is equal to:", simulation_count))
  }
  simulation_density <- kdensity(mean_spearman_values, kernel=c("gaussian"))
  assign("simulation_density_raw", mean_spearman_values, envir=.GlobalEnv)
  assign("simulation_density", simulation_density, envir=.GlobalEnv)
  ########produces a density plot of the simulation
  plot(simulation_density, col="purple", ylim=c(0,20), xlim=c(-0.15,0.15), xlab="Mean Spearman Coefficient", main=NULL) 
  legend(-0.1, 15, legend=c(paste(set_size, " genes")),
         col=c("purple"), lty=1:2, cex=1.1)
}