##BoxPlot comparison 
#Palette for Rotation  WR = "#DED93E", WM = "#196418"
#Palette for Microhabitat  RA = "#b05644", RH = "#d9b967", RP = "#57896A"
#Palette for Layer 1 = "#dcc8ba", Layer 2 = "#DEB3AD", Layer 3 = "#DE847B", Layer 4 = "#B95C50", Layer 5 = "#3B0404"

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")
library("gridExtra")
library("ggpubr")
library("RColorBrewer")
library("agricolae")   

########################################## Compare phylum - Rotation 
#by rotation in every microhabitat
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1, taxonomic.rank = "Phylum"))


#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq1_phylum <- tax_glom(psO_jki_seq1, "Phylum", NArm= TRUE)
psO_jki_seq1_phylum
ntaxa(psO_jki_seq1); ntaxa(psO_jki_seq1_phylum)

#### Create tables 
df_psO_jki_seq1_phylum <- data.frame(tax_table(psO_jki_seq1_phylum),otu_table(psO_jki_seq1_phylum))
write.csv(df_psO_jki_seq1_phylum, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_phylum.csv")

###Transforming absolute abundance to relative abundance
##Using transform_sample_counts()
#relative abundance (%)
#For psO_jki_seq1_phylum
psO_jki_seq1_phylum_rel<-transform_sample_counts(psO_jki_seq1_phylum, function(x) (x*100)/sum(x))
psO_jki_seq1_phylum_rel
head(otu_table(psO_jki_seq1_phylum_rel))

#### Create tables 
df_psO_jki_seq1_phylum_rel <- data.frame(tax_table(psO_jki_seq1_phylum_rel),otu_table(psO_jki_seq1_phylum_rel))
write.csv(df_psO_jki_seq1_phylum_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_phylum_rel.csv")

##Subset by phylum
ps_phylum_high= subset_taxa(psO_jki_seq1_phylum_rel, Phylum == "Proteobacteria" | Phylum == "Actinobacteriota" )
ps_phylum_low= subset_taxa(psO_jki_seq1_phylum_rel, Phylum == "Crenarchaeota" | Phylum == "Acidobacteriota" | Phylum == "Bacteroidota" | Phylum == "Gemmatimonadota" | Phylum == "Verrucomicrobiota")
ps_phylum_signif= subset_taxa(psO_jki_seq1_phylum_rel, Phylum == "Nitrospirota" | Phylum == "Firmicutes")

##Create a function to plot abundances for a subset of samples (Subset Phylum and facet Phylum)
plot_abundance_stat_phylum = function(phyloseq, title = "",
                                      Facet ="Phylum", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    #scale_x_discrete(limits=c("WR","WM")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") +
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
jki_seq1_rotation_colors <- c("WM" = "#196418", "WR" = "#DED93E")



# Plot and run statistics for subset sample
stat_plot_ps_phylum_high <- plot_abundance_stat_phylum(ps_phylum_high, Facet = "Phylum", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20), label = c("0", "20", "40", "60", "80")) 
 # stat_compare_means(method = "wilcox", label= "p", label.y = 78, label.x = 1, size=5)
stat_plot_ps_phylum_high

ggsave("stat_plot_ps_phylum_high.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 10, units = "cm", dpi = 300, device = "png")

stat_plot_ps_phylum_low <- plot_abundance_stat_phylum(ps_phylum_low, Facet = "Phylum", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), label = c("0", "5", "10", "15", "20")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 19.5, label.x = 1, size=5)
stat_plot_ps_phylum_low

ggsave("stat_plot_ps_phylum_low.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 20, height = 10, units = "cm", dpi = 300, device = "png")


stat_plot_ps_phylum_signif <- plot_abundance_stat_phylum(ps_phylum_signif, Facet = "Phylum", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), label = c("0", "2", "4", "6", "8", "10")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 9.5, label.x = 1, size=5)
stat_plot_ps_phylum_signif

ggsave("stat_plot_ps_phylum_signif.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 10, units = "cm", dpi = 300, device = "png")



##Create a function to plot abundances for a subset of samples (Subset Phylum and facet Phylum)
plot_abundance_stat_archaea_phylum = function(phyloseq, title = "",
                                      Facet ="Phylum", Color = "Rotation") {
  p1k_arch = subset_taxa(phyloseq, Kingdom %in% c("Archaea"))
  mkk_arch = psmelt(p1k_arch)
  ggplot(data = mkk_arch, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    #scale_x_discrete(limits=c("WR","WM")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") +
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Plot and run statistics for subset sample
stat_plot_ps_archaea_phylum <- plot_abundance_stat_archaea_phylum(ps_phylum_low, Facet = "Phylum", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
 # theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), label = c("0", "2", "4", "6", "8", "10")) 
  #stat_compare_means(method = "wilcox", label= "p", label.y = 9.5, label.x = 1, size=5)
stat_plot_ps_archaea_phylum

ggsave("stat_plot_ps_archaea_phylum.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 5.4, height = 10, units = "cm", dpi = 300, device = "png")


#####@@@@@@####@@#@#@@#@##@@#@##@ GENERAL PGPR
########################################## Compare Annotation - Rotation by microhabitat
#by microhabitat overall
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_RA_filt, taxonomic.rank = "Annotation"))
length(get_taxa_unique(psO_jki_seq1_RH_filt, taxonomic.rank = "Annotation"))
length(get_taxa_unique(psO_jki_seq1_RP_filt, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq1_RA_filt_annotation <- tax_glom(psO_jki_seq1_RA_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_RA_filt_annotation
ntaxa(psO_jki_seq1_RA_filt); ntaxa(psO_jki_seq1_RA_filt_annotation)

psO_jki_seq1_RH_filt_annotation <- tax_glom(psO_jki_seq1_RH_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_RH_filt_annotation
ntaxa(psO_jki_seq1_RH_filt); ntaxa(psO_jki_seq1_RH_filt_annotation)

psO_jki_seq1_RP_filt_annotation <- tax_glom(psO_jki_seq1_RP_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_RP_filt_annotation
ntaxa(psO_jki_seq1_RP_filt); ntaxa(psO_jki_seq1_RP_filt_annotation)

#### Create tables 
df_psO_jki_seq1_RA_filt_annotation <- data.frame(tax_table(psO_jki_seq1_RA_filt_annotation),otu_table(psO_jki_seq1_RA_filt_annotation))
write.csv(df_psO_jki_seq1_RA_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RA_filt_annotation.csv")

df_psO_jki_seq1_RH_filt_annotation <- data.frame(tax_table(psO_jki_seq1_RH_filt_annotation),otu_table(psO_jki_seq1_RH_filt_annotation))
write.csv(df_psO_jki_seq1_RH_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RH_filt_annotation.csv")

df_psO_jki_seq1_RP_filt_annotation <- data.frame(tax_table(psO_jki_seq1_RP_filt_annotation),otu_table(psO_jki_seq1_RP_filt_annotation))
write.csv(df_psO_jki_seq1_RP_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RP_filt_annotation.csv")

###Transforming absolute abundance to relative abundance
##Using transform_sample_counts()
#relative abundance (%)
psO_jki_seq1_RA_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_RA_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_RA_filt_annotation_rel
head(otu_table(psO_jki_seq1_RA_filt_annotation_rel))

psO_jki_seq1_RH_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_RH_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_RH_filt_annotation_rel
head(otu_table(psO_jki_seq1_RH_filt_annotation_rel))

psO_jki_seq1_RP_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_RP_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_RP_filt_annotation_rel
head(otu_table(psO_jki_seq1_RP_filt_annotation_rel))

#### Create tables 
df_psO_jki_seq1_RA_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_RA_filt_annotation_rel),otu_table(psO_jki_seq1_RA_filt_annotation_rel))
write.csv(df_psO_jki_seq1_RA_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RA_filt_annotation_rel.csv")

df_psO_jki_seq1_RH_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_RH_filt_annotation_rel),otu_table(psO_jki_seq1_RH_filt_annotation_rel))
write.csv(df_psO_jki_seq1_RH_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RH_filt_annotation_rel.csv")

df_psO_jki_seq1_RP_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_RP_filt_annotation_rel),otu_table(psO_jki_seq1_RP_filt_annotation_rel))
write.csv(df_psO_jki_seq1_RP_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_RP_filt_annotation_rel.csv")

####Select taxa 
ps_pgpr_annotation_RA1= subset_taxa(psO_jki_seq1_RA_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Devosia"| Annotation == "Gaiella" | Annotation == "MB-A2-108" | Annotation == "Nitrospira" | Annotation == "ANPR" | Annotation == "Pedobacter")
ps_pgpr_annotation_RA2= subset_taxa(psO_jki_seq1_RA_filt_annotation_rel, Annotation == "Pseudoxanthomonas" | Annotation == "Phyllobacterium" | Annotation == "Bradyrhizobium" | Annotation == "Bacillus" | Annotation == "Pseudolabrys" | Annotation == "Acinetobacter" | Annotation == "Shinella")

ps_pgpr_annotation_RH1= subset_taxa(psO_jki_seq1_RH_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Devosia"| Annotation == "Gaiella" | Annotation == "MB-A2-108" | Annotation == "Nitrospira" | Annotation == "ANPR" | Annotation == "Pedobacter")
ps_pgpr_annotation_RH2= subset_taxa(psO_jki_seq1_RH_filt_annotation_rel, Annotation == "Pseudoxanthomonas" | Annotation == "Phyllobacterium" | Annotation == "Bradyrhizobium" | Annotation == "Bacillus" | Annotation == "Pseudolabrys" | Annotation == "Acinetobacter" | Annotation == "Shinella")

ps_pgpr_annotation_RP1= subset_taxa(psO_jki_seq1_RP_filt_annotation_rel, Annotation == "Devosia" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas" | Annotation == "Nitrospira" | Annotation == "Pedobacter" | Annotation == "Nocardioides" | Annotation == "Streptomyces" | Annotation == "Paenibacillus")
ps_pgpr_annotation_RP2= subset_taxa(psO_jki_seq1_RP_filt_annotation_rel, Annotation == "Shinella" | Annotation == "Sphingopyxis" | Annotation == "Massilia" |  Annotation == "Variovorax" | Annotation == "Agromyces" | Annotation == "Brevundimonas" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium")

# ps_pgpr_WM1= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Devosia" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas" | Annotation == "Nitrospira" | Annotation == "Pedobacter" | Annotation == "Nocardioides" | Annotation == "Streptomyces" | Annotation == "Paenibacillus")
# ps_pgpr_WM2= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Shinella" | Annotation == "Sphingopyxis" | Annotation == "Massilia" |  Annotation == "Variovorax" | Annotation == "Agromyces" | Annotation == "Brevundimonas" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium")


##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)
#for Layer
plot_abundance_stat_rotation_microhabitat = function(phyloseq, title = "",
                                      Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") + 
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
jki_seq1_rotation_colors <- c("WM" = "#196418", "WR" = "#DED93E")

# Plot for Rotation and run statistics for subset sample
#RA
stat_plot_ps_pgpr_annotation_RA1 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RA1, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), label = c("0", "2", "4", "6", "8", "10")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 9.5, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RA1

ggsave("stat_plot_ps_pgpr_annotation_RA1.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 30, height = 10, units = "cm", dpi = 200, device = "png")

stat_plot_ps_pgpr_annotation_RA2 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RA2, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0, 4, by = 1), label = c("0", "1", "2", "3", "4")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 3.8, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RA2

ggsave("stat_plot_ps_pgpr_annotation_RA2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 30, height = 10, units = "cm", dpi = 200, device = "png")


#RH
stat_plot_ps_pgpr_annotation_RH1 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RH1, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), label = c("0", "2", "4", "6", "8", "10")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 9.5, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RH1

ggsave("stat_plot_ps_pgpr_annotation_RH1.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 30, height = 10, units = "cm", dpi = 200, device = "png")

stat_plot_ps_pgpr_annotation_RH2 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RH2, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0, 4, by = 1), label = c("0", "1", "2", "3", "4")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 3.8, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RH2

ggsave("stat_plot_ps_pgpr_annotation_RH2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 30, height = 10, units = "cm", dpi = 200, device = "png")


#RP
stat_plot_ps_pgpr_annotation_RP1a <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RP1, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), label = c("0", "2", "4", "6", "8", "10")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 9.5, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RP1a

ggsave("stat_plot_ps_pgpr_annotation_RP1a.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 30, height = 10, units = "cm", dpi = 200, device = "png")

stat_plot_ps_pgpr_annotation_RP2a <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RP2, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0, 4, by = 1), label = c("0", "1", "2", "3", "4")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 3.8, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RP2a

ggsave("stat_plot_ps_pgpr_annotation_RP2a.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 30, height = 10, units = "cm", dpi = 200, device = "png")

#################################################

#####@@@@@@####@@#@#@@#@##@@#@##@ RESPONDERS DESeq2 ---> >1% relative abundance 

########################################## Compare Annotation - Rotation by microhabitat
#by microhabitat overall
####Select taxa 
ps_pgpr_annotation_RA_deseq2= subset_taxa(psO_jki_seq1_RA_filt_annotation_rel, Annotation == "MB-A2-108" | Annotation == "wb1-A12" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium" | Annotation == "Devosia")

ps_pgpr_annotation_RH_deseq2= subset_taxa(psO_jki_seq1_RH_filt_annotation_rel, Annotation == "MB-A2-108" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium" | Annotation == "Pseudomonas" | Annotation == "Paenibacillus")

ps_pgpr_annotation_RP_deseq2= subset_taxa(psO_jki_seq1_RP_filt_annotation_rel, Annotation == "MB-A2-108" | Annotation == "Gemmatimonadaceae" | Annotation == "Shinella" | Annotation == "Nitrospira" | Annotation == "Pedobacter" | Annotation == "Phyllobacterium"| Annotation == "Gaiellales" | Annotation == "Agromyces")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)
#for Layer
plot_abundance_stat_rotation_microhabitat = function(phyloseq, title = "",
                                                     Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") + 
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
#jki_seq1_rotation_colors <- c("WM" = "#196418", "WR" = "#DED93E")

# Plot for Rotation and run statistics for subset sample
#RA
stat_plot_ps_pgpr_annotation_RA_deseq2 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RA_deseq2, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), label = c("0", "2", "4", "6", "8")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 7.5, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RA_deseq2

ggsave("stat_plot_ps_pgpr_annotation_RA_deseq2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 25, height = 8, units = "cm", dpi = 200, device = "png")

#RH
stat_plot_ps_pgpr_annotation_RH_deseq2 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RH_deseq2, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), label = c("0", "2", "4", "6", "8")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 7.5, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RH_deseq2

ggsave("stat_plot_ps_pgpr_annotation_RH_deseq2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 25, height = 8, units = "cm", dpi = 200, device = "png")



#RP
stat_plot_ps_pgpr_annotation_RP_deseq2 <- plot_abundance_stat_rotation_microhabitat(ps_pgpr_annotation_RP_deseq2, Facet = "Annotation", Color = "Rotation") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1), label = c("0", "1", "2", "3", "4")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 3.7, label.x = 1, size=5)
stat_plot_ps_pgpr_annotation_RP_deseq2

ggsave("stat_plot_ps_pgpr_annotation_RP_deseq2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 40, height = 8, units = "cm", dpi = 200, device = "png")


#####
#####
#####
#####
#####
#####
#####
#####


########## Microhabitat - by rotation
#by rotation overall
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_WM_filt, taxonomic.rank = "Annotation"))
length(get_taxa_unique(psO_jki_seq1_WR_filt, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq1_WM_filt_annotation <- tax_glom(psO_jki_seq1_WM_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_WM_filt_annotation
ntaxa(psO_jki_seq1_WM_filt); ntaxa(psO_jki_seq1_WM_filt_annotation)

psO_jki_seq1_WR_filt_annotation <- tax_glom(psO_jki_seq1_WR_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_WR_filt_annotation
ntaxa(psO_jki_seq1_WR_filt); ntaxa(psO_jki_seq1_WR_filt_annotation)

#### Create tables 
df_psO_jki_seq1_WM_filt_annotation <- data.frame(tax_table(psO_jki_seq1_WM_filt_annotation),otu_table(psO_jki_seq1_WM_filt_annotation))
write.csv(df_psO_jki_seq1_WM_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_annotation.csv")

df_psO_jki_seq1_WR_filt_annotation <- data.frame(tax_table(psO_jki_seq1_WR_filt_annotation),otu_table(psO_jki_seq1_WR_filt_annotation))
write.csv(df_psO_jki_seq1_WR_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_annotation.csv")

###Transforming absolute abundance to relative abundance
##Using transform_sample_counts()
#relative abundance (%)
#For psO_jki_seq1_WM_filt_annotation
psO_jki_seq1_WM_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_WM_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_WM_filt_annotation_rel
head(otu_table(psO_jki_seq1_WM_filt_annotation_rel))

#For psO_jki_seq1_WM_filt_annotation
psO_jki_seq1_WR_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_WR_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_WR_filt_annotation_rel
head(otu_table(psO_jki_seq1_WR_filt_annotation_rel))

#### Create tables 
df_psO_jki_seq1_WM_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_WM_filt_annotation_rel),otu_table(psO_jki_seq1_WM_filt_annotation_rel))
write.csv(df_psO_jki_seq1_WM_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_annotation_rel.csv")

df_psO_jki_seq1_WR_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_WR_filt_annotation_rel),otu_table(psO_jki_seq1_WR_filt_annotation_rel))
write.csv(df_psO_jki_seq1_WR_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_annotation_rel.csv")

####Select taxa
ps_pgpr_WM1= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Devosia" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas" | Annotation == "Nitrospira" | Annotation == "Pedobacter" | Annotation == "Nocardioides" | Annotation == "Streptomyces" | Annotation == "Paenibacillus")
ps_pgpr_WM2= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Shinella" | Annotation == "Sphingopyxis" | Annotation == "Massilia" |  Annotation == "Variovorax" | Annotation == "Agromyces" | Annotation == "Brevundimonas" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium")

ps_pgpr_WR1= subset_taxa(psO_jki_seq1_WR_filt_annotation_rel, Annotation == "Devosia" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas" |Annotation == "Nitrospira" | Annotation == "Pedobacter"| Annotation == "Nocardioides" | Annotation == "Streptomyces" | Annotation == "Paenibacillus" )
ps_pgpr_WR2= subset_taxa(psO_jki_seq1_WR_filt_annotation_rel, Annotation == "Shinella" | Annotation == "Sphingopyxis" | Annotation == "Massilia" |  Annotation == "Variovorax" | Annotation == "Agromyces"  | Annotation == "Brevundimonas" | Annotation == "Bacillus"| Annotation == "Bradyrhizobium")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)
#for Microhabitat
plot_abundance_stat_microhabitat_rotation = function(phyloseq, title = "",
                                         Facet ="Annotation", Color = "Microhabitat") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Microhabitat", y="Abundance", color = Color, fill = Color)) +
    theme_bw() +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") +
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
jki_seq1_microhabitat_colors <- c("RA" = "#b05644", "RH" = "#d9b967", "RP" = "#57896A")

#Select variable of comparison
comparison_microhabitat <- list(c("RA","RH"),c("RA","RP"),c("RH","RP"))
                                
# Plot for Microhabitat
p_abund_stat_microhabitat_pgpr_WR1 <- plot_abundance_stat_microhabitat_rotation(ps_pgpr_WR1, Facet = "Annotation", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 16, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  scale_color_manual(values = jki_seq1_microhabitat_colors) +
  scale_fill_manual(values = jki_seq1_microhabitat_colors) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 3), label = c("0", "3", "6", "9")) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(6.1,6.6,7.1)) 
  #stat_compare_means(method = "kruskal", label= "p", label.y = 8.9, label.x = 1, size=4)
p_abund_stat_microhabitat_pgpr_WR1

ggsave("p_abund_stat_microhabitat_pgpr_WR1.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 35, height = 12, units = "cm", dpi = 300, device = "png")

p_abund_stat_microhabitat_pgpr_WM1 <- plot_abundance_stat_microhabitat_rotation(ps_pgpr_WM1, Facet = "Annotation", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 16, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  scale_color_manual(values = jki_seq1_microhabitat_colors) +
  scale_fill_manual(values = jki_seq1_microhabitat_colors) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 3), label = c("0", "3", "6", "9")) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(6.1,6.6,7.1)) 
  #stat_compare_means(method = "kruskal", label= "p", label.y = 8.9, label.x = 1, size=4)
p_abund_stat_microhabitat_pgpr_WM1

ggsave("p_abund_stat_microhabitat_pgpr_WM1.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 35, height = 12, units = "cm", dpi = 300, device = "png")


p_abund_stat_microhabitat_pgpr_WR2 <- plot_abundance_stat_microhabitat_rotation(ps_pgpr_WR2, Facet = "Annotation", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 16, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  scale_color_manual(values = jki_seq1_microhabitat_colors) +
  scale_fill_manual(values = jki_seq1_microhabitat_colors) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1), label = c("0", "1", "2", "3")) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(2.3,2.5,2.7)) 
  #stat_compare_means(method = "kruskal", label= "p", label.y = 3, label.x = 1, size=4)
p_abund_stat_microhabitat_pgpr_WR2

ggsave("p_abund_stat_microhabitat_pgpr_WR2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 35, height = 12, units = "cm", dpi = 300, device = "png")

p_abund_stat_microhabitat_pgpr_WM2 <- plot_abundance_stat_microhabitat_rotation(ps_pgpr_WM2, Facet = "Annotation", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 16, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  scale_color_manual(values = jki_seq1_microhabitat_colors) +
  scale_fill_manual(values = jki_seq1_microhabitat_colors) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1), label = c("0", "1", "2", "3")) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(2.3,2.5,2.7)) 
  #stat_compare_means(method = "kruskal", label= "p", label.y = 3, label.x = 1, size=4)
p_abund_stat_microhabitat_pgpr_WM2

ggsave("p_abund_stat_microhabitat_pgpr_WM2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 35, height = 12, units = "cm", dpi = 300, device = "png")






########## Layer - by rotation
#by rotation overall
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_WM_filt, taxonomic.rank = "Annotation"))
length(get_taxa_unique(psO_jki_seq1_WR_filt, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq1_WM_filt_annotation <- tax_glom(psO_jki_seq1_WM_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_WM_filt_annotation
ntaxa(psO_jki_seq1_WM_filt); ntaxa(psO_jki_seq1_WM_filt_annotation)

psO_jki_seq1_WR_filt_annotation <- tax_glom(psO_jki_seq1_WR_filt, "Annotation", NArm= TRUE)
psO_jki_seq1_WR_filt_annotation
ntaxa(psO_jki_seq1_WR_filt); ntaxa(psO_jki_seq1_WR_filt_annotation)

#### Create tables 
df_psO_jki_seq1_WM_filt_annotation <- data.frame(tax_table(psO_jki_seq1_WM_filt_annotation),otu_table(psO_jki_seq1_WM_filt_annotation))
write.csv(df_psO_jki_seq1_WM_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_annotation.csv")

df_psO_jki_seq1_WR_filt_annotation <- data.frame(tax_table(psO_jki_seq1_WR_filt_annotation),otu_table(psO_jki_seq1_WR_filt_annotation))
write.csv(df_psO_jki_seq1_WR_filt_annotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_annotation.csv")

###Transforming absolute abundance to relative abundance
##Using transform_sample_counts()
#relative abundance (%)
#For psO_jki_seq1_WM_filt_annotation
psO_jki_seq1_WM_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_WM_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_WM_filt_annotation_rel
head(otu_table(psO_jki_seq1_WM_filt_annotation_rel))

#For psO_jki_seq1_WM_filt_annotation
psO_jki_seq1_WR_filt_annotation_rel<-transform_sample_counts(psO_jki_seq1_WR_filt_annotation, function(x) (x*100)/sum(x))
psO_jki_seq1_WR_filt_annotation_rel
head(otu_table(psO_jki_seq1_WR_filt_annotation_rel))

#### Create tables 
df_psO_jki_seq1_WM_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_WM_filt_annotation_rel),otu_table(psO_jki_seq1_WM_filt_annotation_rel))
write.csv(df_psO_jki_seq1_WM_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WM_filt_annotation_rel.csv")

df_psO_jki_seq1_WR_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq1_WR_filt_annotation_rel),otu_table(psO_jki_seq1_WR_filt_annotation_rel))
write.csv(df_psO_jki_seq1_WR_filt_annotation_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_WR_filt_annotation_rel.csv")

####Select taxa
ps_pgpr_WR1_deseq2= subset_taxa(psO_jki_seq1_WR_filt_annotation_rel, Annotation == "Acinetobacter" | Annotation == "Agromyces" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium"| Annotation == "Gaiellales")
ps_pgpr_WR2_deseq2= subset_taxa(psO_jki_seq1_WR_filt_annotation_rel, Annotation == "Gemmatimonadaceae" | Annotation == "MB-A2-108"| Annotation == "Nitrospira" | Annotation == "Paenibacillus" | Annotation == "Pedobacter")
ps_pgpr_WR3_deseq2= subset_taxa(psO_jki_seq1_WR_filt_annotation_rel, Annotation == "Phyllobacterium" | Annotation == "Pseudomonas" | Annotation == "Rhodococcus" | Annotation == "Shinella" | Annotation == "wb1-A12")

ps_pgpr_WR4_deseq2= subset_taxa(psO_jki_seq1_WR_filt_annotation_rel, Annotation == "Agromyces" | Annotation == "Bacillus" | Annotation == "Gaiellales" | Annotation == "MB-A2-108" | Annotation == "Nitrospira" | Annotation == "Paenibacillus" | Annotation == "Pedobacter" | Annotation == "Phyllobacterium" |Annotation == "Rhodococcus" | Annotation == "Shinella")

ps_pgpr_WM1_deseq2= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Acinetobacter" | Annotation == "Agromyces" | Annotation == "Bacillus" | Annotation == "Bradyrhizobium"| Annotation == "Gaiellales")
ps_pgpr_WM2_deseq2= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Gemmatimonadaceae" | Annotation == "MB-A2-108"| Annotation == "Nitrospira" | Annotation == "Paenibacillus" | Annotation == "Pedobacter")
ps_pgpr_WM3_deseq2= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Phyllobacterium" | Annotation == "Pseudomonas" | Annotation == "Rhodococcus" | Annotation == "Shinella" | Annotation == "wb1-A12")

ps_pgpr_WM4_deseq2= subset_taxa(psO_jki_seq1_WM_filt_annotation_rel, Annotation == "Agromyces" | Annotation == "Bacillus" | Annotation == "Gaiellales" | Annotation == "MB-A2108" | Annotation == "Nitrospira" | Annotation == "Paenibacillus" | Annotation == "Pedobacter" | Annotation == "Phyllobacterium" |Annotation == "Rhodococcus" | Annotation == "Shinella")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)
#for Layer
plot2_abundance_stat_depth_WR = function(phyloseq, title = "",
                                         Facet ="Annotation", Color = "Layer") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Microhabitat", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, ncol = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#3B0404","#B95C50","#DE847B", "#DEB3AD","#dcc8ba")) +
    scale_fill_manual(values =c("#3B0404","#B95C50","#DE847B", "#DEB3AD","#dcc8ba"))
  
}

# plot2_abundance_stat_depth_WM = function(phyloseq, title = "",
#                                          Facet ="Annotation", Color = "Layer") {
#   p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
#   mkk_bac = psmelt(p1k_bac)
#   ggplot(data = mkk_bac, mapping = aes_string(x ="Layer", y="Abundance", color = Color, fill = Color)) +
#     geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
#     facet_wrap(facets = Facet, ncol = 1) +
#     theme(legend.position = "right") +
#     scale_color_manual(values =c("#b05644", "#d9b967", "#57896A")) +
# #     scale_fill_manual(values =c("#b05644", "#d9b967", "#57896A"))
#   
# }
# Plot for Layer
p_abund_stat_depth_pgpr_WR1 <- plot2_abundance_stat_depth_WR(ps_pgpr_WR1_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WR1

ggsave("p_abund_stat_depth_pgpr_WR1.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm",dpi = 300)

p_abund_stat_depth_pgpr_WR2 <- plot2_abundance_stat_depth_WR(ps_pgpr_WR2_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WR2

ggsave("p_abund_stat_depth_pgpr_WR2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm",dpi = 300)

p_abund_stat_depth_pgpr_WR3 <- plot2_abundance_stat_depth_WR(ps_pgpr_WR3_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 5.5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WR3

ggsave("p_abund_stat_depth_pgpr_WR3.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm",dpi = 300)


p_abund_stat_depth_pgpr_WR4 <- plot2_abundance_stat_depth_WR(ps_pgpr_WR4_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 5.5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WR4

ggsave("p_abund_stat_depth_pgpr_WR4.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 32, units = "cm",dpi = 300)

# Plot for Layer
p_abund_stat_depth_pgpr_WM1 <- plot2_abundance_stat_depth_WR(ps_pgpr_WM1_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WM1

ggsave("p_abund_stat_depth_pgpr_WM1.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm",dpi = 300)

p_abund_stat_depth_pgpr_WM2 <- plot2_abundance_stat_depth_WR(ps_pgpr_WM2_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WM2

ggsave("p_abund_stat_depth_pgpr_WM2.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm",dpi = 300)

p_abund_stat_depth_pgpr_WM3 <- plot2_abundance_stat_depth_WR(ps_pgpr_WM3_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 5.5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WM3

ggsave("p_abund_stat_depth_pgpr_WM3.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 28, units = "cm",dpi = 300)


p_abund_stat_depth_pgpr_WM4 <- plot2_abundance_stat_depth_WR(ps_pgpr_WM4_deseq2, Facet = "Annotation", Color = "Layer") +
  theme_bw() + 
  ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(size=20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 5.5, size=3) +
  theme(strip.text.x = element_text(size = 14, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abund_stat_depth_pgpr_WM4

ggsave("p_abund_stat_depth_pgpr_WM4.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 32, units = "cm",dpi = 300)




########################################## Compare phylum - Microhabitat

#by rotation in every microhabitat   ------> WM
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_WM_filt_phylum, taxonomic.rank = "Phylum"))

###Transforming absolute abundance to relative abundance  
##Using transform_sample_counts()
#relative abundance (%)
#For psO_jki_seq1_phylum
psO_jki_seq1_WM_phylum_rel<-transform_sample_counts(psO_jki_seq1_WM_filt_phylum, function(x) (x*100)/sum(x))
psO_jki_seq1_WM_phylum_rel
head(otu_table(psO_jki_seq1_WM_phylum_rel))

#### Create tables 
df_psO_jki_seq1_phylum_rel <- data.frame(tax_table(psO_jki_seq1_WM_phylum_rel),otu_table (psO_jki_seq1_WM_phylum_rel))
write.csv(df_psO_jki_seq1_phylum_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_phylum_rel.csv")

##Subset by phylum
ps_phylum_WM_microh_high= subset_taxa(psO_jki_seq1_WM_phylum_rel, Phylum == "Proteobacteria" | Phylum == "Actinobacteriota" )
ps_phylum_WM_microh_low= subset_taxa(psO_jki_seq1_WM_phylum_rel, Phylum == "Acidobacteriota" | Phylum == "Bacteroidota" | Phylum == "Chloroflexi" | Phylum == "Gemmatimonadota" | Phylum == "Verrucomicrobiota" | Phylum == "Nitrospirota" | Phylum == "Firmicutes" | Phylum == "Myxococcota")

##Create a function to plot abundances for a subset of samples (Subset Phylum and facet Phylum)
plot_abundance_stat_phylum = function(phyloseq, title = "",
                                      Facet ="Phylum", Color = "Microhabitat") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Microhabitat", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    #scale_x_discrete(limits=c("WR","WM")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") +
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
jki_seq1_rotation_colors <- c("RA" = "#b05644", "RH" = "#d9b967", "RP" = "#57896A")

#Select variable of comparison
comparison_microhabitat <- list(c("RA","RH"),c("RA","RP"),c("RH","RP"))

# Plot and run statistics for subset sample
stat_plot_ps_WM_phylum_microh_high <- plot_abundance_stat_phylum(ps_phylum_WM_microh_high, Facet = "Phylum", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(60,65,70)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20), label = c("0", "20", "40", "60", "80")) 
  #stat_compare_means(method = "wilcox", label= "p", label.y = 78, label.x = 1, size=5)
stat_plot_ps_WM_phylum_microh_high

ggsave("stat_plot_ps_WM_phylum_microh_high.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 10, units = "cm", dpi = 300, device = "png")


stat_plot_ps_WM_phylum_microh_low <- plot_abundance_stat_phylum(ps_phylum_WM_microh_low, Facet = "Phylum", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(15,16,17)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), label = c("0", "5", "10", "15", "20")) 
#stat_compare_means(method = "wilcox", label= "p", label.y = 78, label.x = 1, size=5)
stat_plot_ps_WM_phylum_microh_low

ggsave("stat_plot_ps_WM_phylum_microh_low.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 35, height = 10, units = "cm", dpi = 300, device = "png")


#by rotation in every microhabitat   ------> WR
#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq1_WR_filt_phylum, taxonomic.rank = "Phylum"))

###Transforming absolute abundance to relative abundance  
##Using transform_sample_counts()
#relative abundance (%)
#For psO_jki_seq1_phylum
psO_jki_seq1_WR_phylum_rel<-transform_sample_counts(psO_jki_seq1_WR_filt_phylum, function(x) (x*100)/sum(x))
psO_jki_seq1_WR_phylum_rel
head(otu_table(psO_jki_seq1_WR_phylum_rel))

#### Create tables 
df_psO_jki_seq1_phylum_rel <- data.frame(tax_table(psO_jki_seq1_WR_phylum_rel),otu_table (psO_jki_seq1_WR_phylum_rel))
write.csv(df_psO_jki_seq1_phylum_rel, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_phylum_rel.csv")

##Subset by phylum
ps_phylum_WR_microh_high= subset_taxa(psO_jki_seq1_WR_phylum_rel, Phylum == "Proteobacteria" | Phylum == "Actinobacteriota" )
ps_phylum_WR_microh_low= subset_taxa(psO_jki_seq1_WR_phylum_rel, Phylum == "Acidobacteriota" | Phylum == "Bacteroidota" | Phylum == "Chloroflexi" | Phylum == "Gemmatimonadota" | Phylum == "Verrucomicrobiota" | Phylum == "Nitrospirota" | Phylum == "Firmicutes" | Phylum == "Myxococcota")


##Create a function to plot abundances for a subset of samples (Subset Phylum and facet Phylum)
plot_abundance_stat_phylum = function(phyloseq, title = "",
                                      Facet ="Phylum", Color = "Microhabitat") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Microhabitat", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    #scale_x_discrete(limits=c("WR","WR")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") +
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
jki_seq1_rotation_colors <- c("RA" = "#b05644", "RH" = "#d9b967", "RP" = "#57896A")

#Select variable of comparison
comparison_microhabitat <- list(c("RA","RH"),c("RA","RP"),c("RH","RP"))

# Plot and run statistics for subset sample
stat_plot_ps_WR_phylum_microh_high <- plot_abundance_stat_phylum(ps_phylum_WR_microh_high, Facet = "Phylum", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(60,65,70)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20), label = c("0", "20", "40", "60", "80")) 
#stat_compare_means(method = "wilcox", label= "p", label.y = 78, label.x = 1, size=5)
stat_plot_ps_WR_phylum_microh_high

ggsave("stat_plot_ps_WR_phylum_microh_high.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 10, height = 10, units = "cm", dpi = 300, device = "png")


stat_plot_ps_WR_phylum_microh_low <- plot_abundance_stat_phylum(ps_phylum_WR_microh_low, Facet = "Phylum", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 12, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16, colour = "black")) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(15,16,17)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), label = c("0", "5", "10", "15", "20")) 
#stat_compare_means(method = "wilcox", label= "p", label.y = 78, label.x = 1, size=5)
stat_plot_ps_WR_phylum_microh_low

ggsave("stat_plot_ps_WR_phylum_microh_low.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 35, height = 10, units = "cm", dpi = 300, device = "png")



##Subset by phylum
ps_phylum_WM_archaea_microh= subset_taxa(psO_jki_seq1_WM_phylum_rel, Phylum == "Crenarchaeota" )
ps_phylum_WR_archaea_microh= subset_taxa(psO_jki_seq1_WR_phylum_rel, Phylum == "Crenarchaeota" )


##Create a function to plot abundances for a subset of samples (Subset Phylum and facet Phylum)
plot_abundance_stat_archaea_phylum_microh = function(phyloseq, title = "",
                                              Facet ="Phylum", Color = "Microhabitat") {
  p1k_arch = subset_taxa(phyloseq, Kingdom %in% c("Archaea"))
  mkk_arch = psmelt(p1k_arch)
  ggplot(data = mkk_arch, mapping = aes_string(x ="Microhabitat", y="Abundance", color = Color, fill = Color)) +
    theme_bw() + 
    #scale_x_discrete(limits=c("WR","WM")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") +
    geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Plot and run statistics for subset sample
stat_plot_ps_WM_archaea_phylum_microh <- plot_abundance_stat_archaea_phylum_microh(ps_phylum_WM_archaea_microh, Facet = "Phylum", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  # theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), label = c("0", "2", "4", "6", "8")) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(6,6.5,7))
stat_plot_ps_WM_archaea_phylum_microh

ggsave("stat_plot_ps_WM_archaea_phylum_microh.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 5, height = 10, units = "cm", dpi = 300, device = "png")



stat_plot_ps_WR_archaea_phylum_microh <- plot_abundance_stat_archaea_phylum_microh(ps_phylum_WR_archaea_microh, Facet = "Phylum", Color = "Microhabitat") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_text(size=14, color = "black"), axis.text.y = element_text(size=14, color = "black")) +
  # theme(legend.text = element_text(size = 20), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 14, face="italic", color = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_color_manual(values = jki_seq1_rotation_colors) +
  scale_fill_manual(values = jki_seq1_rotation_colors) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), label = c("0", "2", "4", "6", "8")) +
  stat_compare_means(comparisons = comparison_microhabitat, method = "wilcox.test", label= "p", bracket.size = .1, size=4, label.y = c(6,6.5,7))
stat_plot_ps_WR_archaea_phylum_microh

ggsave("stat_plot_ps_WR_archaea_phylum_microh.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Stats_abundance_jki_seq1/", width = 5, height = 10, units = "cm", dpi = 300, device = "png")


#The end! Enjoy your day!  : )
