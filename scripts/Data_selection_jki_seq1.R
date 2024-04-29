###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

###Sub-setting samples from cleaned dataset to keep only samples that represent every case on the variable  
##Using subset_samples()
#Rotation

psO_jki_seq1_WR <- subset_samples(psO_jki_seq1, Rotation=="WR")
psO_jki_seq1_WR

psO_jki_seq1_WM <- subset_samples(psO_jki_seq1, Rotation=="WM")
psO_jki_seq1_WM

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

psO_jki_seq1_WR_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WR) > 0, psO_jki_seq1_WR)
psO_jki_seq1_WR_filt

psO_jki_seq1_WM_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WM) > 0, psO_jki_seq1_WM)
psO_jki_seq1_WM_filt

####Create tables

df_psO_jki_seq1_WR_filt <- data.frame(tax_table(psO_jki_seq1_WR_filt),otu_table(psO_jki_seq1_WR_filt))
write.csv(df_psO_jki_seq1_WR_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt.csv")

df_psO_jki_seq1_WM_filt <- data.frame(tax_table(psO_jki_seq1_WM_filt),otu_table(psO_jki_seq1_WM_filt))
write.csv(df_psO_jki_seq1_WM_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt.csv")

###Agglomerating taxa per sample at the appropriated taxonomic rank
##Using tax_glom()
#Annotation

colnames(tax_table(psO_jki_seq1_WR_filt))
psO_jki_seq1_WR_filt_annotation <- tax_glom(psO_jki_seq1_WR_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_WR_filt); ntaxa(psO_jki_seq1_WR_filt_annotation)

colnames(tax_table(psO_jki_seq1_WM_filt))
psO_jki_seq1_WM_filt_annotation <- tax_glom(psO_jki_seq1_WM_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_WM_filt); ntaxa(psO_jki_seq1_WM_filt_annotation)

### Transform to Relative abundance

psO_jki_seq1_WR_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_WR_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_WR_filt_annotation_rel
head(otu_table(psO_jki_seq1_WR_filt_annotation_rel))

psO_jki_seq1_WM_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_WM_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_WM_filt_annotation_rel
head(otu_table(psO_jki_seq1_WM_filt_annotation_rel))

###Create tables

df_psO_jki_seq1_WR_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_WR_filt_annotation_rel),otu_table(psO_jki_seq1_WR_filt_annotation_rel))
write.csv(df_psO_jki_seq1_WR_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_annotation_rel1.csv")

df_psO_jki_seq1_WM_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_WM_filt_annotation_rel),otu_table(psO_jki_seq1_WM_filt_annotation_rel))
write.csv(df_psO_jki_seq1_WM_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_annotation_rel1.csv")

################################################################
##Using subset_samples() subset Rotation data by Microhabitat to compare layers

psO_jki_seq1_WR_filt_RA <- subset_samples(psO_jki_seq1_WR_filt, Microhabitat=="RA")
psO_jki_seq1_WR_filt_RA
psO_jki_seq1_WR_filt_RH <- subset_samples(psO_jki_seq1_WR_filt, Microhabitat=="RH")
psO_jki_seq1_WR_filt_RH
psO_jki_seq1_WR_filt_RP <- subset_samples(psO_jki_seq1_WR_filt, Microhabitat=="RP")
psO_jki_seq1_WR_filt_RP

psO_jki_seq1_WM_filt_RA <- subset_samples(psO_jki_seq1_WM_filt, Microhabitat=="RA")
psO_jki_seq1_WM_filt_RA
psO_jki_seq1_WM_filt_RH <- subset_samples(psO_jki_seq1_WM_filt, Microhabitat=="RH")
psO_jki_seq1_WM_filt_RH
psO_jki_seq1_WM_filt_RP <- subset_samples(psO_jki_seq1_WM_filt, Microhabitat=="RP")
psO_jki_seq1_WM_filt_RP

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

psO_jki_seq1_WR_filt_RA_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WR_filt_RA) > 0, psO_jki_seq1_WR_filt_RA)
psO_jki_seq1_WR_filt_RA_filt
psO_jki_seq1_WR_filt_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WR_filt_RH) > 0, psO_jki_seq1_WR_filt_RH)
psO_jki_seq1_WR_filt_RH_filt
psO_jki_seq1_WR_filt_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WR_filt_RP) > 0, psO_jki_seq1_WR_filt_RP)
psO_jki_seq1_WR_filt_RP_filt

psO_jki_seq1_WM_filt_RA_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WM_filt_RA) > 0, psO_jki_seq1_WM_filt_RA)
psO_jki_seq1_WM_filt_RA_filt
psO_jki_seq1_WM_filt_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WM_filt_RH) > 0, psO_jki_seq1_WM_filt_RH)
psO_jki_seq1_WM_filt_RH_filt
psO_jki_seq1_WM_filt_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq1_WM_filt_RP) > 0, psO_jki_seq1_WM_filt_RP)
psO_jki_seq1_WM_filt_RP_filt

#### Create tables

df_psO_jki_seq1_WR_filt_RA_filt <- data.frame(tax_table(psO_jki_seq1_WR_filt_RA_filt),otu_table(psO_jki_seq1_WR_filt_RA_filt))
write.csv(df_psO_jki_seq1_WR_filt_RA_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_RA_filt.csv")
df_psO_jki_seq1_WR_filt_RH_filt <- data.frame(tax_table(psO_jki_seq1_WR_filt_RH_filt),otu_table(psO_jki_seq1_WR_filt_RH_filt))
write.csv(df_psO_jki_seq1_WR_filt_RH_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_RH_filt.csv")
df_psO_jki_seq1_WR_filt_RP_filt <- data.frame(tax_table(psO_jki_seq1_WR_filt_RP_filt),otu_table(psO_jki_seq1_WR_filt_RP_filt))
write.csv(df_psO_jki_seq1_WR_filt_RP_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_RP_filt.csv")

df_psO_jki_seq1_WM_filt_RA_filt <- data.frame(tax_table(psO_jki_seq1_WM_filt_RA_filt),otu_table(psO_jki_seq1_WM_filt_RA_filt))
write.csv(df_psO_jki_seq1_WM_filt_RA_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_RA_filt.csv")
df_psO_jki_seq1_WM_filt_RH_filt <- data.frame(tax_table(psO_jki_seq1_WM_filt_RH_filt),otu_table(psO_jki_seq1_WM_filt_RH_filt))
write.csv(df_psO_jki_seq1_WM_filt_RH_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_RH_filt.csv")
df_psO_jki_seq1_WM_filt_RP_filt <- data.frame(tax_table(psO_jki_seq1_WM_filt_RP_filt),otu_table(psO_jki_seq1_WM_filt_RP_filt))
write.csv(df_psO_jki_seq1_WM_filt_RP_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_RP_filt.csv")

###############################################################
###Sub-setting samples from cleaned dataset to keep only samples that represent every case on the variable  
##Using subset_samples()
#Microhabitat

sample_variables(psO_jki_seq1)
psO_jki_seq1_RA <- subset_samples(psO_jki_seq1, Microhabitat=="RA")
psO_jki_seq1_RA

psO_jki_seq1_RH <- subset_samples(psO_jki_seq1, Microhabitat=="RH")
psO_jki_seq1_RH

psO_jki_seq1_RP <- subset_samples(psO_jki_seq1, Microhabitat=="RP")
psO_jki_seq1_RP

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

psO_jki_seq1_RA_filt<-prune_taxa(taxa_sums(psO_jki_seq1_RA) > 0, psO_jki_seq1_RA)
psO_jki_seq1_RA_filt

psO_jki_seq1_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq1_RH) > 0, psO_jki_seq1_RH)
psO_jki_seq1_RH_filt

psO_jki_seq1_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq1_RP) > 0, psO_jki_seq1_RP)
psO_jki_seq1_RP_filt

#### Create tables 

df_psO_jki_seq1_RA_filt <- data.frame(tax_table(psO_jki_seq1_RA_filt),otu_table(psO_jki_seq1_RA_filt))
write.csv(df_psO_jki_seq1_RA_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RA_filt.csv")

df_psO_jki_seq1_RH_filt <- data.frame(tax_table(psO_jki_seq1_RH_filt),otu_table(psO_jki_seq1_RH_filt))
write.csv(df_psO_jki_seq1_RH_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RH_filt.csv")

df_psO_jki_seq1_RP_filt <- data.frame(tax_table(psO_jki_seq1_RP_filt),otu_table(psO_jki_seq1_RP_filt))
write.csv(df_psO_jki_seq1_RP_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RP_filt.csv")

###Agglomerating taxa per sample at the appropriated taxonomic rank
##Using tax_glom()
#Annotation

colnames(tax_table(psO_jki_seq1_RA_filt))
psO_jki_seq1_RA_filt_annotation <- tax_glom(psO_jki_seq1_RA_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_RA_filt); ntaxa(psO_jki_seq1_RA_filt_annotation)

colnames(tax_table(psO_jki_seq1_RH_filt))
psO_jki_seq1_RH_filt_annotation <- tax_glom(psO_jki_seq1_RH_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_RH_filt); ntaxa(psO_jki_seq1_RH_filt_annotation)

colnames(tax_table(psO_jki_seq1_RP_filt))
psO_jki_seq1_RP_filt_annotation <- tax_glom(psO_jki_seq1_RP_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_RP_filt); ntaxa(psO_jki_seq1_RP_filt_annotation)

###Transform to Relative abundance

psO_jki_seq1_RA_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_RA_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_RA_filt_annotation_rel
head(otu_table(psO_jki_seq1_RA_filt_annotation_rel))

psO_jki_seq1_RH_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_RH_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_RH_filt_annotation_rel
head(otu_table(psO_jki_seq1_RH_filt_annotation_rel))

psO_jki_seq1_RP_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_RP_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_RP_filt_annotation_rel
head(otu_table(psO_jki_seq1_RP_filt_annotation_rel))

####Create tables 

df_psO_jki_seq1_RA_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_RA_filt_annotation_rel),otu_table(psO_jki_seq1_RA_filt_annotation_rel))
write.csv(df_psO_jki_seq1_RA_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RA_filt_annotation_rel1.csv")

df_psO_jki_seq1_RH_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_RH_filt_annotation_rel),otu_table(psO_jki_seq1_RH_filt_annotation_rel))
write.csv(df_psO_jki_seq1_RH_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RH_filt_annotation_rel1.csv")

df_psO_jki_seq1_RP_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_RP_filt_annotation_rel),otu_table(psO_jki_seq1_RP_filt_annotation_rel))
write.csv(df_psO_jki_seq1_RP_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RP_filt_annotation_rel1.csv")

#############################################################################
###Sub-setting samples from cleaned dataset to keep only samples that represent every case on the variable  
##Using subset_samples()
#Layer
psO_jki_seq1_layer1 <- subset_samples(psO_jki_seq1, Layer=="L1")
psO_jki_seq1_layer1

psO_jki_seq1_layer2 <- subset_samples(psO_jki_seq1, Layer=="L2")
psO_jki_seq1_layer2

psO_jki_seq1_layer3 <- subset_samples(psO_jki_seq1, Layer=="L3")
psO_jki_seq1_layer3

psO_jki_seq1_layer4 <- subset_samples(psO_jki_seq1, Layer=="L4")
psO_jki_seq1_layer4

psO_jki_seq1_layer5 <- subset_samples(psO_jki_seq1, Layer=="L5")
psO_jki_seq1_layer5

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()
psO_jki_seq1_layer1_filt<-prune_taxa(taxa_sums(psO_jki_seq1_layer1) > 0, psO_jki_seq1_layer1)
psO_jki_seq1_layer1_filt

psO_jki_seq1_layer2_filt<-prune_taxa(taxa_sums(psO_jki_seq1_layer2) > 0, psO_jki_seq1_layer2)
psO_jki_seq1_layer2_filt

psO_jki_seq1_layer3_filt<-prune_taxa(taxa_sums(psO_jki_seq1_layer3) > 0, psO_jki_seq1_layer3)
psO_jki_seq1_layer3_filt

psO_jki_seq1_layer4_filt<-prune_taxa(taxa_sums(psO_jki_seq1_layer4) > 0, psO_jki_seq1_layer4)
psO_jki_seq1_layer4_filt

psO_jki_seq1_layer5_filt<-prune_taxa(taxa_sums(psO_jki_seq1_layer5) > 0, psO_jki_seq1_layer5)
psO_jki_seq1_layer5_filt

#### Create tables 
df_psO_jki_seq1_layer1_filt <- data.frame(tax_table(psO_jki_seq1_layer1_filt),otu_table(psO_jki_seq1_layer1_filt))
write.csv(df_psO_jki_seq1_layer1_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_layer1_filt.csv")

df_psO_jki_seq1_layer2_filt <- data.frame(tax_table(psO_jki_seq1_layer2_filt),otu_table(psO_jki_seq1_layer2_filt))
write.csv(df_psO_jki_seq1_layer2_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_layer2_filt.csv")

df_psO_jki_seq1_layer3_filt <- data.frame(tax_table(psO_jki_seq1_layer3_filt),otu_table(psO_jki_seq1_layer3_filt))
write.csv(df_psO_jki_seq1_layer3_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_layer3_filt.csv")

df_psO_jki_seq1_layer4_filt <- data.frame(tax_table(psO_jki_seq1_layer4_filt),otu_table(psO_jki_seq1_layer4_filt))
write.csv(df_psO_jki_seq1_layer4_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_layer4_filt.csv")

df_psO_jki_seq1_layer5_filt <- data.frame(tax_table(psO_jki_seq1_layer5_filt),otu_table(psO_jki_seq1_layer5_filt))
write.csv(df_psO_jki_seq1_layer5_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_layer5_filt.csv")




#go through the next scripts ("Alpha_diversity.R")






