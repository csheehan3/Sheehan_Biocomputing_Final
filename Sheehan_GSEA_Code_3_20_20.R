######## You may need to initially download these bioinformatics packages for the analyses 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("fgsea")
BiocManager::install("reactome.db")
install.packages("devtools")
######## THIS WILL DOWNLOAD A LARGE FILE >300 MB TO YOUR DIRECTORY!!!
download.file("https://gdc.xenahubs.net/download/TCGA-BRCA.htseq_fpkm.tsv.gz", "TCGA-BRCA.htseq_fpkm.tsv", method = "auto", quiet=FALSE) 
library(devtools)
library(tidyverse)
library("ggpubr")
library("readr")
library("biomaRt")
library("rlist")
######## Ensemble_ID gene names are converted to Entrez ID names for the fgsea package 
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
BRCA_data <- read_tsv(file = "TCGA-BRCA.htseq_fpkm-uq.tsv")
row_total <- dim(BRCA_data)[1]
Ensembl_ID_list <- c()
for (ensemblename in BRCA_data[1:row_total, 1]){ #gene names are collected into a list vector 
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
for (translat in BRCA_data$Ensembl_ID[1:row_total]){
  cleaved_translat <- str_replace(translat, pattern = ".[0-9]+$", replacement = "")
  newline <- filter(Entrez_conversion, Entrez_conversion$ensembl_gene_id==cleaved_translat)
  n <- n + 1 #####this is just a counter to let you know how its running
  if (n %% 100 == 0) {
    print(paste("Working on line:",n," of ",row_total)) 
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
BRCA_genes_altered <- cbind(as.data.frame(entrez_dictionary), BRCA_genes)
########Now with all genes having their conversion name, we can manipulate the TCGA data for GSEA analysis
global_mean_gene_expression <- filter(BRCA_genes_altered, BRCA_genes_altered$entrez_dictionary==571)[1, 3:1219] %>% #gives you the global mean of the gene of interest, 571 is entrez ID of BACH1
  sum() / 1219
gene_row <- filter(BRCA_genes_altered, BRCA_genes_altered$entrez_dictionary==571)[1,3:1219] 
high_gene_samples <- BRCA_genes_altered[1:60483, c(2, 2 + which(gene_row > global_mean_gene_expression))] #gets the indices for samples with high gene expression
low_gene_samples <- BRCA_genes_altered[1:60483, c(2, 2 + which(gene_row < global_mean_gene_expression))] #same for low gene samples
high_gene_rows <- dim(high_gene_samples)[2]
low_gene_rows <- dim(low_gene_samples)[2]
average_high_list <- c(apply(high_gene_samples[1:60483, 2:high_gene_rows],1,mean))
average_low_list <- c(apply(low_gene_samples[1:60483, 2:low_gene_rows],1,mean))
difference_in_expression <- average_high_list - average_low_list #produces vector containing difference in mean expression for every gene of high and low gene samples
entrez_expression_dataframe <- data.frame(entrez_dictionary, difference_in_expression) %>% #organizes everything into a dataframe without NA genes, and ranked by difference
  filter(., is.na(.[1:60483,1])==FALSE) %>% 
  .[order(-.$difference_in_expression),]
#######Proper format for the fgsea package uses a ranked, named vector with the Entrez ID and difference in expression
ranks_vector <- entrez_expression_dataframe$difference_in_expression
names(ranks_vector) <- as.character(entrez_expression_dataframe$entrez_dictionary)
library(fgsea)
library(reactome.db)
pathways <- reactomePathways(as.character(entrez_expression_dataframe$entrez_dictionary)) #grabs potential gene set pathways based on your query
fgseaRes <- fgsea(pathways, ranks_vector, nperm=1000, maxSize=500) #runs the analysis, ignore the Warning message
#######Here we produce can produce a number of graphs from the output of the GSEA analysis.
#######The first is a figure that displays the 15 most significant lists for positive and negative enrichments
#######with the BACH1 high tumor samples. Next are two enrichments plots for two different lists of interest. 
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(pathways[topPathways], ranks_vector, fgseaRes, 
              gseaParam = 0.5)

plotEnrichment(pathways[["Cell-Cell communication"]], ranks_vector) + labs(title="Cell-Cell communication")

plotEnrichment(pathways[["Transcriptional regulation of pluripotent stem cells"]], ranks_vector) + labs(title="Stem Cell Regulators")

filter(fgseaRes, pathway=="Transcriptional regulation of pluripotent stem cells") %>% head()
########Here I create a new set of genes to run the GSEA test against for BACH1 high and BACH1 low tumors,
########using biomart to grab all the genes from the gene ontology (G0) term that I'm interested in.
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
stemness_GO_conversion <- getBM(attributes = c('hgnc_symbol',
                                               'entrezgene_id'),
                                filters = "go_parent_term", 
                                values = "GO:0019827",
                                mart = mart)
GO_stem_cell_maintenance_list <- list(c(as.character(stemness_GO_conversion$entrezgene_id)))
names(GO_stem_cell_maintenance_list) <- 'GO:0019827'
GO_stem_cell_fgsea <- fgsea(GO_stem_cell_maintenance_list, ranks_vector, nperm=1000, maxSize=500)
########Here we can produce an enrichment plot & look at the statistical output

plotEnrichment(GO_stem_cell_maintenance_list[['GO:0019827']], ranks_vector) + labs(title="Stem Cell Population Maintenance")

head(GO_stem_cell_fgsea)
########Thank you for a great quarter!

