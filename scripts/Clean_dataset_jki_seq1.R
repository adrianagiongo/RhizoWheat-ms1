##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

####################################### 
#Show available ranks in the dataset  --> original dataset = 36 Samples 
rank_names(jki_seq1_data)
sample_names(jki_seq1_data)
jki_seq1_data

####################################### 
#Show available ranks in the dataset
rank_names(jki_seq1_data)

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level (Keeping Archaea)
#Kingdom
table(tax_table(jki_seq1_data)[,"Kingdom"], exclude = NULL)
psO_jki_seq1 <- subset_taxa(jki_seq1_data, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota", "NA", ""))
table(tax_table(psO_jki_seq1)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_jki_seq1)[,"Phylum"], exclude = NULL)
psO_jki_seq1 <- subset_taxa(psO_jki_seq1, !Phylum %in% c("Bacteria", "Archaea"))
table(tax_table(psO_jki_seq1)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_jki_seq1)[,"Order"], exclude = NULL)
psO_jki_seq1 <- subset_taxa(psO_jki_seq1, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(psO_jki_seq1)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_jki_seq1)[,"Family"], exclude = NULL)
psO_jki_seq1 <- subset_taxa(psO_jki_seq1, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(psO_jki_seq1)[,"Family"], exclude = NULL)

#Genus
#table(tax_table(psO_jki_seq1)[,"Genus"], exclude = NULL)

#Annotation
#table(tax_table(psO_jki_seq1)[,"Annotation"], exclude = NULL)

##Observe psO after clean spurious taxa
jki_seq1_data
psO_jki_seq1

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(psO_jki_seq1),
               MARGIN = ifelse(taxa_are_rows(psO_jki_seq1), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(psO_jki_seq1), tax_table(psO_jki_seq1), otu_table(psO_jki_seq1))

#Write data frame in csv format
write.csv(prevdf, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/prevdf_psO_jki_seq1.csv")

#go through the next scripts ("Data_selection.R")

