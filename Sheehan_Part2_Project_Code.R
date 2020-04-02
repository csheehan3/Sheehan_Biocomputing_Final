
######Get these packages or load them 
library("kdensity")
library(devtools)
library(tidyverse)
library("ggpubr")
library("biomaRt")
library("rlist")
######need to get rid of the genes with average value 0
BRCA_data <- read_tsv(file = "TCGA_example_data")
BRCA_data <- cbind(data.frame(average_values = apply(dplyr::select(BRCA_data, matches("TCGA")),1,mean)), BRCA_data)
BRCA_data <- filter(BRCA_data, BRCA_data$average_values!=0) ##these two steps are to eliminate genes filled with zeroes from the data
##### run 500 simulations of correlating one random gene with a set of 100 other random genes
##### this is done with simulation function. Input is simulation_function(the_dataframe, size of gene set, number of simulations)
source("simulation_function.R")
simulation_function(BRCA_data, 100, 500)


##### Here I add the Entrez gene name conversions to the BRCA dataset
##### copied from the fgsea code for translating Ensemble ID names
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast") ##servers can get overloaded, change mirror to circumvent this issue
Ensembl_ID_list <- c()
for (ensemblename in dplyr::select(BRCA_data, matches("Ensembl"))){ #gene names are collected into a list vector 
  Ensembl_ID_list <- c(ensemblename, Ensembl_ID_list)
}
cleaved_gene_names <- str_replace(Ensembl_ID_list, #Ensemble_ID gene names are modified to a more general ID name
                                  pattern = ".[0-9]+$",
                                  replacement = "")
Entrez_conversion <- getBM(attributes = c('ensembl_gene_id', #query sent to biomaRt and desired Entrez names return
                                          'entrezgene_id'),
                           filters = 'ensembl_gene_id', 
                           values = cleaved_gene_names,
                           mart = mart)
######## To deal with incomplete data and cases having multiple matched Entrez IDs, we iterate over gene rows 
######## in the BRCA_data table to match the Entrez IDs or fill in NAs for missing data
n=0
entrez_dictionary <- c()
for (translat in BRCA_data$Ensembl_ID[1: nrow(BRCA_data)]){
  cleaved_translat <- str_replace(translat, pattern = ".[0-9]+$", replacement = "")
  newline <- filter(Entrez_conversion, Entrez_conversion$ensembl_gene_id==cleaved_translat)
  n <- n + 1 #####this is just a counter to let you know how its running
  if (n %% 100 == 0) {
    print(paste("Working on line:",n," of ", nrow(BRCA_data))) 
  }
  if (nrow(newline) > 1){
    newline <- newline[1, 1:2] ###picks the first for cases with multiple gene transcripts
  }
  if (length(newline$entrezgene_id)==0 || is.na(newline$entrezgene_id)==TRUE) { ###puts in NA for cases that missed a conversion
    entrez_dictionary <- c(entrez_dictionary, NA)
  } else {
    entrez_dictionary <- c(entrez_dictionary, newline$entrezgene_id) ####grabs the Entrez conversion for all other cases
  }
}
BRCA_data <- cbind(as.data.frame(entrez_dictionary), BRCA_data)
###### Here I put together a list of metabolic GO terms that i was interested in exploring,
###### along with all the genes classified under that GO term. The GO terms were manually collected
###### online from the GO database; however, all the genes classified under that term had to be
###### collected from biomaRt. Because biomaRt is slow & the query is fairly large, I ran this code
###### already & provided the results as an R object to load. The code for making this object 
###### can be found in GO_term_collecting_code.R
load("GO_metabolic_list.Rdata")


#####calculate the mean of the multiple correlations tests
BACH1_row <- filter(BRCA_data, BRCA_data$entrez_dictionary==571) %>% ##571 is the Entrez ID for bach1
  dplyr::select(matches("TCGA")) %>% 
  as.numeric()
GO_spearman_values <- c()
mean_GO_spearman_values <- c()
for (term in 1:length(GO_metabolic_list)) { ## iterates over GO terms
  filtered_frame <- filter(BRCA_data, BRCA_data$entrez_dictionary %in% GO_metabolic_list[[term]]) ##grabs rows for genes included in the GO term bein worked on
  for (n in 1:nrow(filtered_frame)){ ## iterates over the genes for 1 GO term
    test_row <- filtered_frame[n,] %>%
      dplyr::select(matches("TCGA")) %>% 
      as.numeric()
    GO_spearman_values <- c(GO_spearman_values, cor(BACH1_row, test_row, method=c("spearman"))) ## add coefficient value to a list
  }
  mean_GO_spearman_values <- c(mean_GO_spearman_values, mean(GO_spearman_values)) ## add mean from last list of coefficients to a list of means
  GO_spearman_values <- c() ##reset the GO_spearman_values list for the next GO term to fill up
  print(paste("completed:", names(GO_metabolic_list)[term])) #reporter
}
GO_tibble <- tibble("Metabolic_GO_Term" = names(GO_metabolic_list), "Mean_Spearman_Value" = mean_GO_spearman_values)
arrange(GO_tibble, desc(GO_tibble$Mean_Spearman_Value))
######now for a probability value,
######first I need to find spearman value where integration to +Inf is = 0.5
starting_value=0.03
while (integrate(function(x) simulation_density(x), lower=starting_value, upper=Inf)$value > 0.5){ #increases value until integration is equal to 0.5
  starting_value = starting_value + 0.00001
}
Probability_function = function(X1){ #Function that passes value to integration function depending on whether it is above or below the midway point
  if (X1 > starting_value){
    integrate(function(x) simulation_density(x), lower=X1, upper=Inf)$value
  } else {
    integrate(function(x) simulation_density(x), lower=-Inf, upper=X1)$value
  }
}
GO_tibble <- sapply(GO_tibble$Mean_Spearman_Value, Probability_function) %>% ##Gets the probability from each mean and adds to tibble
  cbind(GO_tibble, "Probability_Value"=.) 
Bunch_of_zeros <- GO_tibble$Mean_Spearman_Value * 0 #used for plotting
GO_tibble <- cbind(GO_tibble, "density"=Bunch_of_zeros)
###### This piece generates a plot with the density from the simulations with random data along with the actual data from GO terms
###### Had issues with aesthetic mapping when trying to overlay density plot + scatter plot from different datasets.
###### To get around this, I just generated data from the density function & plotted the two sets of data as a
###### scatter plot & an area plot, so no conflicts with the aesthetics
density_data <- data.frame("density"=simulation_density(-0.3), "Spearman_Value"= -0.3)
n=1
stepsize = 0.6 / 300 ### 0.6 is the range of spearman coefficients I am interested in, 300 how many datapoints I want
while (n < 300){
  density_data <- density_data %>% add_row(density=simulation_density(-0.3 + (n*stepsize)), Spearman_Value= -0.3 + (n*stepsize))
  n = n + 1
}
cols <- c('red', 'blue')
ggplot(data=GO_tibble, aes(x=Mean_Spearman_Value, y=density)) +
  geom_point(aes(color="red"), size=1.5) +
  geom_area(data=density_data, aes(x=Spearman_Value, y=density, fill="blue"), alpha=0.3) +
  theme_minimal() +
  scale_colour_manual(name='GO term data', values=c("red")) +
  scale_fill_manual(name='Simulation Density', values=c("blue")) +
  xlab("Mean Spearman Coefficient") +
  ylab("Density")
####### I want to see who the top hit is
Top_hit <- arrange(GO_tibble, desc(GO_tibble$Probability_Value)) %>% 
  tail(n=1)
OXPHOS <- as.character(Top_hit[1])
####### Now I can look at the tree plot for that particular GO term. 
####### Too generate the tree plot for that GO term, can simply pass it on to the Tree_Plot_Function.
####### Input is Tree_Plot_Function(dataframe, the entrez ID for gene, the GO term to test against)
####### Output is the bar plot ranked by P-value, with each gene represented as a bar on the X axis
source("Tree_Plot_Function.R")
Tree_Plot_Function(BRCA_data, 571, OXPHOS)



###### This last part generates a tree plot of random data from simulation
set.seed(14)
##
spearman_values <- c()
P_values <- c()
A_very_interesting_gene_index <- sample(1:nrow(BRCA_data), 1, replace=TRUE)
A_very_interesting_gene <- BRCA_data[A_very_interesting_gene_index,] %>% 
  dplyr::select(matches("TCGA")) %>% 
  as.numeric()
row_indices <- sample(1:nrow(BRCA_data), 400, replace=TRUE)
random_rows <- BRCA_data[row_indices,] %>% 
  dplyr::select(matches("TCGA")) 
for (A_gene in 1:400){
  A_correlation <- cor.test(as.numeric(random_rows[A_gene,]), A_very_interesting_gene, method = c("spearman"))
  spearman_values <- c(spearman_values, A_correlation$estimate)
  P_values <- c(P_values, A_correlation$p.value)
}
random_rows <- BRCA_data[row_indices,]
Simulation_tibble <- tibble("Gene_Name"=random_rows$entrez_dictionary, "Spearman_Coefficient"=spearman_values, "P_Value"=P_values)
Simulation_tibble <- arrange(Simulation_tibble, Simulation_tibble$`P_Value`)
ggplot(data=subset(Simulation_tibble, !is.na(Simulation_tibble$Gene_Name)), aes(x=reorder(Gene_Name, P_Value), y=Spearman_Coefficient)) +
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="Simulation Tree Plot", y="Spearman Coefficient", x="P-Value")



