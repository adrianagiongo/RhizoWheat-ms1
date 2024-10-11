##Select taxa on the dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

####################################### 
#Show available ranks in the dataset 
rank_names(psO_jki_seq1)
sample_names(psO_jki_seq1)
psO_jki_seq1

####################################### 
#Show available ranks in the dataset
rank_names(psO_jki_seq1)

##Create table with number of features for Myxococcota only
##Keeping only the specific taxonomic level
#Kingdom
table(tax_table(psO_jki_seq1)[,"Kingdom"], exclude = NULL)
psO_jki_seq1_Myxococcota <- subset_taxa(psO_jki_seq1, !is.na(Kingdom) & !Kingdom %in% c("Archaea"))
table(tax_table(psO_jki_seq1_Myxococcota)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_jki_seq1_Myxococcota)[,"Phylum"], exclude = NULL)
psO_jki_seq1_Myxococcota <- subset_taxa(psO_jki_seq1_Myxococcota, Phylum %in% c("Myxococcota"))
table(tax_table(psO_jki_seq1_Myxococcota)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_jki_seq1_Myxococcota)[,"Order"], exclude = NULL)
#psO_jki_seq1 <- subset_taxa(psO_jki_seq1, !is.na(Order) & !Order %in% c("Chloroplast"))
#table(tax_table(psO_jki_seq1)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_jki_seq1_Myxococcota)[,"Family"], exclude = NULL)
#psO_jki_seq1 <- subset_taxa(psO_jki_seq1, !is.na(Family) & !Family %in% c("Mitochondria"))
#table(tax_table(psO_jki_seq1)[,"Family"], exclude = NULL)

#Genus
table(tax_table(psO_jki_seq1_Myxococcota)[,"Genus"], exclude = NULL)

#Annotation
table(tax_table(psO_jki_seq1_Myxococcota)[,"Annotation"], exclude = NULL)

##Observe psO after select target taxa
psO_jki_seq1
psO_jki_seq1_Myxococcota

#Compute prevalence of each feature and store as data.frame
prevdf_myxococcota = apply(X = otu_table(psO_jki_seq1_Myxococcota),
               MARGIN = ifelse(taxa_are_rows(psO_jki_seq1_Myxococcota), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_myxococcota = data.frame(Prevalence=prevdf_myxococcota, TotalAbundance = taxa_sums(psO_jki_seq1_Myxococcota), tax_table(psO_jki_seq1_Myxococcota), otu_table(psO_jki_seq1_Myxococcota))
head(prevdf_myxococcota)

#Write data frame in csv format
write.csv(prevdf_myxococcota, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/prevdf_myxococcota_psO_jki_seq1_Myxococcota.csv")

###############################  Plot
#Plot Myxococcota

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_Myxococcota, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq1_Myxococcota_annotation <- tax_glom(psO_jki_seq1_Myxococcota, "Annotation", NArm= TRUE)
psO_jki_seq1_Myxococcota_annotation

###Transforming absolute abundance to relative abundance
##Using transform_sample_counts() - relative abundance (%)
psO_jki_seq1_Myxococcota_annotation_rel<-transform_sample_counts(psO_jki_seq1_Myxococcota_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_Myxococcota_annotation_rel
head(otu_table(psO_jki_seq1_Myxococcota_annotation_rel))

#### Create tables 
df_psO_jki_seq1_Myxococcota_annotation_rel <- data.frame(tax_table(psO_jki_seq1_Myxococcota_annotation_rel),otu_table(psO_jki_seq1_Myxococcota_annotation_rel))
write.csv(df_psO_jki_seq1_Myxococcota_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_Myxococcota_annotation_rel.csv")


#Load library
library(data.table)

# create dataframe from phyloseq object
dat_psO_jki_seq1_Myxococcota_annotation_rel <- data.table(psmelt(psO_jki_seq1_Myxococcota_annotation_rel))
# convert Annotation to a character vector from a factor because R
dat_psO_jki_seq1_Myxococcota_annotation_rel$Annotation <- as.character(dat_psO_jki_seq1_Myxococcota_annotation_rel$Annotation)
# group dataframe by Annotation, calculate median rel. abundance
dat_psO_jki_seq1_Myxococcota_annotation_rel[, median := median(Abundance, na.rm = TRUE), 
    by = "Annotation"]

# Plot and run statistics for subset sample
stat_psO_jki_seq1_Myxococcota_annotation_rel <- ggplot(dat_psO_jki_seq1_Myxococcota_annotation_rel[Abundance > 0], aes(x = reorder(Annotation, Abundance, Fun = median, .desc = FALSE), y = Abundance)) +
  geom_boxplot(fill = "#0099f8") +
  #geom_jitter(position=position_jitter(0.3)) +
  ylab("Relative abundance (%)") +
  xlab("Myxococcota") +
  theme_bw(base_size=14) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 10), label = c("0", "10", "20", "30", "40")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=14, color = "black", angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black"))+
  coord_flip()
  #scale_y_log10()
stat_psO_jki_seq1_Myxococcota_annotation_rel

ggsave("stat_psO_jki_seq1_Myxococcota_annotation_rel.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm", dpi = 300, device = "png")






#############################
#############################
###Sub-setting samples from cleaned dataset to keep only samples that represent every case on the variable  
##Using subset_samples()
#Rotation

psO_jki_seq1_Myxococcota_WR <- subset_samples(psO_jki_seq1_Myxococcota, Rotation=="WR")
psO_jki_seq1_Myxococcota_WR

psO_jki_seq1_Myxococcota_WM <- subset_samples(psO_jki_seq1_Myxococcota, Rotation=="WM")
psO_jki_seq1_Myxococcota_WM

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

psO_jki_seq1_Myxococcota_WR_filt<-prune_taxa(taxa_sums(psO_jki_seq1_Myxococcota_WR) > 0, psO_jki_seq1_Myxococcota_WR)
psO_jki_seq1_Myxococcota_WR_filt

psO_jki_seq1_Myxococcota_WM_filt<-prune_taxa(taxa_sums(psO_jki_seq1_Myxococcota_WM) > 0, psO_jki_seq1_Myxococcota_WM)
psO_jki_seq1_Myxococcota_WM_filt

####Create tables

df_psO_jki_seq1_Myxococcota_WR_filt <- data.frame(tax_table(psO_jki_seq1_Myxococcota_WR_filt),otu_table(psO_jki_seq1_Myxococcota_WR_filt))
write.csv(df_psO_jki_seq1_Myxococcota_WR_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_Myxococcota_WR_filt.csv")

df_psO_jki_seq1_Myxococcota_WM_filt <- data.frame(tax_table(psO_jki_seq1_Myxococcota_WM_filt),otu_table(psO_jki_seq1_Myxococcota_WM_filt))
write.csv(df_psO_jki_seq1_Myxococcota_WM_filt, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_Myxococcota_WM_filt.csv")

###Agglomerating taxa per sample at the appropriated taxonomic rank
##Using tax_glom()
#Annotation

colnames(tax_table(psO_jki_seq1_Myxococcota_WR_filt))
psO_jki_seq1_Myxococcota_WR_filt_annotation <- tax_glom(psO_jki_seq1_Myxococcota_WR_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_Myxococcota_WR_filt); ntaxa(psO_jki_seq1_Myxococcota_WR_filt_annotation)

colnames(tax_table(psO_jki_seq1_Myxococcota_WM_filt))
psO_jki_seq1_Myxococcota_WM_filt_annotation <- tax_glom(psO_jki_seq1_Myxococcota_WM_filt, taxrank = "Annotation")
ntaxa(psO_jki_seq1_Myxococcota_WM_filt); ntaxa(psO_jki_seq1_Myxococcota_WM_filt_annotation)

### Transform to Relative abundance

psO_jki_seq1_Myxococcota_WR_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_Myxococcota_WR_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_Myxococcota_WR_filt_annotation_rel
head(otu_table(psO_jki_seq1_Myxococcota_WR_filt_annotation_rel))

psO_jki_seq1_Myxococcota_WM_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_Myxococcota_WM_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_Myxococcota_WM_filt_annotation_rel
head(otu_table(psO_jki_seq1_Myxococcota_WM_filt_annotation_rel))

###Create tables

df_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_Myxococcota_WR_filt_annotation_rel),otu_table(psO_jki_seq1_Myxococcota_WR_filt_annotation_rel))
write.csv(df_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel1.csv")

df_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_Myxococcota_WM_filt_annotation_rel),otu_table(psO_jki_seq1_Myxococcota_WM_filt_annotation_rel))
write.csv(df_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel1.csv")

#######

###############################  Plot
#Plot Myxococcota WR

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_Myxococcota_WR_filt_annotation_rel, taxonomic.rank = "Annotation"))

#Load library
library(data.table)

# create dataframe from phyloseq object
dat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel <- data.table(psmelt(psO_jki_seq1_Myxococcota_WR_filt_annotation_rel))
# convert Annotation to a character vector from a factor because R
dat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel$Annotation <- as.character(dat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel$Annotation)
# group dataframe by Annotation, calculate median rel. abundance
dat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel[, median := median(Abundance, na.rm = TRUE), 
                                                    by = "Annotation"]

# Plot and run statistics for subset sample
stat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel <- ggplot(dat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel[Abundance > 0], aes(x = reorder(Annotation, Abundance, Fun = median, .desc = FALSE), y = Abundance)) +
  geom_boxplot(fill = "#0099f8") +
  #geom_jitter(position=position_jitter(0.3)) +
  ylab("Relative abundance (%)") +
  xlab("Myxococcota") +
  theme_bw(base_size=14) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 10), label = c("0", "10", "20", "30", "40")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=14, color = "black", angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black"))+
  coord_flip()
#scale_y_log10()
stat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel

ggsave("stat_psO_jki_seq1_Myxococcota_WR_filt_annotation_rel.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm", dpi = 300, device = "png")



#Plot Myxococcota WM

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_Myxococcota_WM_filt_annotation_rel, taxonomic.rank = "Annotation"))

#Load library
library(data.table)

# create dataframe from phyloseq object
dat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel <- data.table(psmelt(psO_jki_seq1_Myxococcota_WM_filt_annotation_rel))
# convert Annotation to a character vector from a factor because R
dat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel$Annotation <- as.character(dat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel$Annotation)
# group dataframe by Annotation, calculate median rel. abundance
dat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel[, median := median(Abundance, na.rm = TRUE), 
                                                    by = "Annotation"]

# Plot and run statistics for subset sample
stat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel <- ggplot(dat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel[Abundance > 0], aes(x = reorder(Annotation, Abundance, Fun = median, .desc = FALSE), y = Abundance)) +
  geom_boxplot(fill = "#0099f8") +
  #geom_jitter(position=position_jitter(0.3)) +
  ylab("Relative abundance (%)") +
  xlab("Myxococcota") +
  theme_bw(base_size=14) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 10), label = c("0", "10", "20", "30", "40")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=14, color = "black", angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20)) +
  theme(legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black"))+
  coord_flip()
#scale_y_log10()
stat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel

ggsave("stat_psO_jki_seq1_Myxococcota_WM_filt_annotation_rel.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm", dpi = 300, device = "png")

