##Creating heatmap for selected dataset

#Loading package
library("microbiome")
library("phyloseq")
library("dplyr")
library("ggpubr")
library("vegan")


#Define colors
jki_seq2_seq4_phy_colors <- c("Actinobacteriota" = "#ff8b60", 
                         "Proteobacteria" = "#363781",
                    "Acidobacteriota" = "#ffc331",
                    "Firmicutes" = "#979191",
                    "Chloroflexi" = "#2b8f22",
                    "Nitrospirota" = "#1361cf",
                    "Crenarchaeota" = "#90c541",
                    "Gemmatimonadota" = "#ca250f",
                    "Verrucomicrobiota" = "#ced2ce",
                    "Bacteroidota"  = "#Bff4be",
                    "Myxococcota" = "#FFC0CB", 
                    "Others" = "#3a3845")




################################################################################
# Plots geral

#count number of samples on dataset using "nsamples" function, and check the number of OTUs on the dataset
nsamples(psO_jki_seq2_seq4)
psO_jki_seq2_seq4

#create function aggregate_rare_taxa
aggregate_rare_taxa <- function (x, level, detection, prevalence, include.lowest=FALSE, ...) {
  x <- aggregate_taxa(x, level)
  rare <- rare_members(x, detection, prevalence, include.lowest)
  tax <- tax_table(x)
  inds <- which(rownames(tax) %in% rare)
  tax[inds, level] <- "Others"
  tax_table(x) <- tax
  tt <- tax_table(x)[, level]
  tax_table(x) <- tax_table(tt)
  aggregate_taxa(x, level)
}

##Aggregate rare taxa
#define total sums
total_samples_psO_jki_seq2_seq4 <- phyloseq::nsamples(psO_jki_seq2_seq4)
total_sum_psO_psO_jki_seq2_seq4 = sample_sums(psO_jki_seq2_seq4)

#aggregate based on phylum (detection higher than 40 copies and present in more than 85% of samples)
psO_jki_seq2_seq4_aggreg_rare_phy <- aggregate_rare_taxa(psO_jki_seq2_seq4, level="Phylum", detection = 1000, prevalence =60/total_samples_psO_jki_seq2_seq4)
psO_jki_seq2_seq4_aggreg_rare_phy

##transform absolute to relative abundance using "transform" function
#for phylum
psO_jki_seq2_seq4_aggreg_rare_phy_rel <- microbiome::transform(psO_jki_seq2_seq4_aggreg_rare_phy, "compositional")
psO_jki_seq2_seq4_aggreg_rare_phy_rel

##plot using "plot_composition" function from microbiome package (plot,type = "barplot" or "heatmap")
###################

#### AVERAGE Rotation
plot_bar_jki_seq2_seq4_rotation <- plot_composition(psO_jki_seq2_seq4_aggreg_rare_phy_rel, otu.sort = "abundance", x.label = "Site", plot.type = "barplot", average_by = "Site", verbose = FALSE) +
  theme_bw() +
  labs(x="", y= "Relative abundance (%)") +
  scale_fill_manual(values = jki_seq2_seq4_phy_colors, breaks=c('Proteobacteria', 'Actinobacteriota', 'Acidobacteriota', 'Bacteroidota', 'Chloroflexi', 
                                                           'Firmicutes', 'Verrucomicrobiota', 'Gemmatimonadota', 'Nitrospirota', 'Crenarchaeota', 'Myxococcota', 'Others')) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), label = c("0", "25", "50", "75", "100")) +
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 16, colour = "black"), axis.title.y = element_text(size = 16, colour = "black")) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 16, colour = "black"))
plot_bar_jki_seq2_seq4_rotation

ggsave("plot_bar_jki_seq2_seq4_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Compositional_barplot_jki_seq2_seq4/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")


## The end! Enjoy!  :)

