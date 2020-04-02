
######## You may need to initially download these bioinformatics packages for the analyses 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if (!requireNamespace("fgsea", quietly = TRUE))
  install.packages("fgsea")
if (!requireNamespace("reactome.db", quietly = TRUE))
  install.packages("reactome.db")
if (!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("GEOquery")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(GEOquery)
######## THIS WILL DOWNLOAD A LARGE FILE >300 MB TO YOUR DIRECTORY!!!
if (file.exists("TCGA-BRCA_fpkm.tsv")==FALSE){
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA.htseq_fpkm.tsv.gz", "TCGA-BRCA.htseq_fpkm.tsv", method = "auto", quiet=FALSE) ##Takes around 2min
  gunzip("TCGA-BRCA.htseq_fpkm.tsv", destname = gsub("[.]gz$", "", "TCGA-BRCA_fpkm.tsv"), overwrite=TRUE) 
}
library(devtools)
library(tidyverse)
library("ggpubr")
library("readr")
library("biomaRt")
library("rlist")
######## Ensemble_ID gene names are converted to Entrez ID names for the fgsea package 
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast") ##servers can get overloaded, change mirror to circumvent this issue
BRCA_data <- read_tsv(file = "TCGA-BRCA_fpkm.tsv")
Ensembl_ID_list <- c()
for (ensemblename in BRCA_data[1:nrow(BRCA_data), 1]){ #gene names are collected into a list vector 
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
######  Now with all genes having their conversion name, we can manipulate the TCGA data for GSEA analysis.
######  This is now done by the fgsea_function.
######  Inputs for the fgsea_function(dataframe, gene-of-interest) will perform an fgsea analysis after inputing the 
######  cleaned-up dataframe & the entrez ID of the gene you are interested in looking at. Output is the top positive 
######  and top negative pathways associated with high-versus-low expression of the gene of interest.
######  Example: 571 is the ID number for the transcription factor, BACH1
source("fgsea_function.R")
fgsea_function(BRCA_data, 571)
########Here I create a new set of genes to run the GSEA test against for BACH1 high and BACH1 low tumors,
########using biomart to grab all the genes from the gene ontology (G0) term that I'm interested in.
stemness_GO_conversion <- getBM(attributes = c('hgnc_symbol',
                                               'entrezgene_id'),
                                filters = "go_parent_term", 
                                values = "GO:0019827",
                                mart = mart)
GO_stem_cell_maintenance_list <- list(c(as.character(stemness_GO_conversion$entrezgene_id)))
names(GO_stem_cell_maintenance_list) <- 'GO:0019827'
GO_stem_cell_fgsea <- fgsea(GO_stem_cell_maintenance_list, ranks_vector, nperm=1000, maxSize=500)
head(GO_stem_cell_fgsea)
########Here we can produce an enrichment plot & look at the statistical output
plotEnrichment(GO_stem_cell_maintenance_list[['GO:0019827']], ranks_vector) + labs(title="Stem Cell Population Maintenance")
########Thank you for a great quarter!

