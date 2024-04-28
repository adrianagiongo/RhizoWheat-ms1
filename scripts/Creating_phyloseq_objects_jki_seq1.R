#Transform dataset to phyloseq format

#Clear R's brain
rm(list = ls())

#load libraries
library("phyloseq")
library("ggplot2")
library("readxl")
library("dplyr")

#move excel files to R objects and metadata on .csv format with names on row 1. 
otu_mat<-read_excel("~/Documents/R_analysis/jki_seq1/data_jki_seq1/JKI_seq1_otu.xlsx")
taxo_mat<-read_excel("~/Documents/R_analysis/jki_seq1/data_jki_seq1/JKI_seq1_taxa.xlsx")
metadata<-read.csv("~/Documents/R_analysis/jki_seq1/data_jki_seq1/JKI_seq1_metadata.csv", row.names = 1)

#define the row names from the otu column under header "OTU"
#remove the column otu from otu, taxo and metadata tables since it is now used as a row name
row.names(otu_mat)<-otu_mat$OTU
otu_mat_unamed<-otu_mat %>% select (-OTU)

row.names(taxo_mat)<-taxo_mat$OTU
taxo_mat_unamed<-taxo_mat %>% select (-OTU)

#transform otu and taxo tables into matrix
otu_matrix<-as.matrix(otu_mat_unamed)
taxo_matrix<-as.matrix(taxo_mat_unamed)

#transform to phyloseq objects using functions otu_table, tax_table and sample_data
OTU=otu_table(otu_matrix, taxa_are_rows = TRUE)
TAX=tax_table(taxo_matrix)
SAMPLE=sample_data(metadata, errorIfNULL = T)

jki_seq1_data<-phyloseq(OTU, TAX, SAMPLE)

#visualize phyloseq-class
jki_seq1_data

#visualize data
sample_names(jki_seq1_data)
rank_names(jki_seq1_data)
sample_variables(jki_seq1_data)

#go through the next scripts ("Rename_NA_taxa.R")

