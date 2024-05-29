#Create rarefied data from phyloseq object using vegan package or microbiome package
#Palette for Rotation  WR = "#B1AD31", WM = "#196418" 
#Palette for Microhabitat  RA = "#b05644", RH = "#d9b967", RP = "#57896A"
#Palette for Layer 1 = "#3B0404", Layer 2 = "#B95C50", Layer 3 = "#DE847B", Layer 4 = "#DEB3AD", Layer 5 = "#dcc8ba"

#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")

#Keep theme_bw for the rest of the R session
theme_set(theme_bw())

##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1)), step=20, cex=0.5, col = "blue")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq1_rarefied<-rarefy_even_depth(psO_jki_seq1, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1)),trimOTUs=TRUE)
psO_jki_seq1_rarefied
sample_sums(psO_jki_seq1_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alfa diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq1_rarefied_1 <- richness(psO_jki_seq1_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_rarefied_2 <- evenness(psO_jki_seq1_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_rarefied_3 <- microbiome::diversity(psO_jki_seq1_rarefied, index = 'shannon', zeroes = TRUE)

##Boxplot alpha diversity using ggviolin function from ggpubr
#Prepare file using meta function
psO_jki_seq1_rarefied.meta <- meta(psO_jki_seq1_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for all
psO_jki_seq1_rarefied.meta$observed <- alfa_div_psO_jki_seq1_rarefied_1$observed
psO_jki_seq1_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_rarefied_1$chao1
psO_jki_seq1_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_rarefied_2$pielou
psO_jki_seq1_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_rarefied_3$shannon

head(psO_jki_seq1_rarefied.meta)
write.csv(psO_jki_seq1_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_rarefied_meta.csv")

##check for the distribution of the diversity using "hist" function to plot, 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05) 
#then the null hypotesis (the population is normally distributed) is rejected
#for data with only two variable

hist(psO_jki_seq1_rarefied.meta$shannon)
shapiro.test(psO_jki_seq1_rarefied.meta$shannon)
qqnorm(psO_jki_seq1_rarefied.meta$shannon)

hist(psO_jki_seq1_rarefied.meta$chao1)
shapiro.test(psO_jki_seq1_rarefied.meta$chao1)
qqnorm(psO_jki_seq1_rarefied.meta$chao1)

hist(psO_jki_seq1_rarefied.meta$pielou)
shapiro.test(psO_jki_seq1_rarefied.meta$pielou)
qqnorm(psO_jki_seq1_rarefied.meta$pielou)

####################################### Overall comparison WR vs WM 

# #Select variable of comparison
# alpha_div_comparison_rotation <- list(c("WR", "WM"))

###plot using ggviolin
##test kruskall_wallis
#to compare rotations Shannon
#no color, colorful dots by rotation
violin_plot_psO_jki_seq1_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_rarefied.meta, x="Rotation", y= "shannon", ylab = "Shannon index", xlab = "") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font("legend.text", size = 20) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_x_discrete(limits=c("WM","WR")) +
  scale_y_continuous(limits = c(2.0, 8.0), breaks = seq(2.0, 8.0, by = 2.0), label = c("2.0", "4.0", "6.0", "8.0"))
  #stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_rarefied_diversity_rotation.png", path = "~/path/", width = 10, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_rarefied.meta, x="Rotation", y= "chao1", ylab = "Chao1 index", xlab = "") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font("legend.text", size = 20) +
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_x_discrete(limits=c("WM","WR")) +
  scale_y_continuous(limits = c(400, 4400), breaks = seq(400, 4400, by = 1000), label = c("400", "1400", "2400", "3400", "4400"))
  #stat_compare_means(method = "wilcox", label= "p", label.y = 4300, size=6)
violin_plot_psO_jki_seq1_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_rarefied_richness_rotation.png", path = "~/path/", width = 10, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_rarefied.meta, x="Rotation", y= "pielou", ylab = "Pielou index", xlab = "") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font("legend.text", size = 22) +
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_x_discrete(limits=c("WM","WR")) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, by = 0.25), label = c("0.5", "0.75", "1.0"))
  #stat_compare_means(method = "wilcox", label= "p", label.y = 0.97, size=6)
violin_plot_psO_jki_seq1_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_rarefied_evenness_rotation.png", path = "~/path/", width = 10, height = 12, units = "cm", dpi = 300, device = "png")


#######################################     For each Rotation groups

##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_WR_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq1_WM_filt)), step=20, cex=0.5, col = "red")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq1_WR_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_WR_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_WR_filt)),trimOTUs=TRUE)
psO_jki_seq1_WR_filt_rarefied
sample_sums(psO_jki_seq1_WR_filt_rarefied)

psO_jki_seq1_WM_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_WM_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_WM_filt)),trimOTUs=TRUE)
psO_jki_seq1_WM_filt_rarefied
sample_sums(psO_jki_seq1_WM_filt_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_WR_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_WM_filt_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq1_WR_filt_rarefied_1 <- richness(psO_jki_seq1_WR_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_WM_filt_rarefied_1 <- richness(psO_jki_seq1_WM_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq1_WR_filt_rarefied_2 <- evenness(psO_jki_seq1_WR_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_WM_filt_rarefied_2 <- evenness(psO_jki_seq1_WM_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq1_WR_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_WR_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_WM_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_WM_filt_rarefied, index = 'shannon', zeroes = TRUE)

##Boxplot alpha diversity using ggviolin function from ggpubr
#Prepare file using meta function
psO_jki_seq1_WR_filt_rarefied.meta <- meta(psO_jki_seq1_WR_filt_rarefied)
psO_jki_seq1_WM_filt_rarefied.meta <- meta(psO_jki_seq1_WM_filt_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for WR
psO_jki_seq1_WR_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_WR_filt_rarefied_1$observed
psO_jki_seq1_WR_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_WR_filt_rarefied_1$chao1
psO_jki_seq1_WR_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_WR_filt_rarefied_2$pielou
psO_jki_seq1_WR_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_WR_filt_rarefied_3$shannon

head(psO_jki_seq1_WR_filt_rarefied.meta)
write.csv(psO_jki_seq1_WR_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_WR_filt_rarefied_meta.csv")

#for WM
psO_jki_seq1_WM_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_WM_filt_rarefied_1$observed
psO_jki_seq1_WM_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_WM_filt_rarefied_1$chao1
psO_jki_seq1_WM_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_WM_filt_rarefied_2$pielou
psO_jki_seq1_WM_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_WM_filt_rarefied_3$shannon

head(psO_jki_seq1_WM_filt_rarefied.meta)
write.csv(psO_jki_seq1_WM_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_WM_filt_rarefied.csv")

########## Compare layers in rotations
#Select variable of comparison
alfa_div_comparison_layer <- list(c("L1","L2"),c("L1","L3"),c("L1","L4"),c("L1","L5"),c("L2","L3"),c("L2","L4"),c("L2","L5"),c("L3","L4"),c("L3","L5"),c("L4","L5"))
alfa_div_comparison_layer_WR_shannon <- list(c("L1","L4"),c("L1","L5"),c("L2","L5"))
alfa_div_comparison_layer_WR_chao1 <- list(c("L1","L4"),c("L1","L5"))
alfa_div_comparison_layer_WR_pielou <- list(c("L1","L4"), c("L1","L5"), c("L2","L5"))
alfa_div_comparison_layer_WM_shannon <- list(c("L1","L3"),c("L1","L4"),c("L1","L5"),c("L2","L5"),c("L3","L5"))
alfa_div_comparison_layer_WM_chao1 <- list(c("L1","L2"),c("L1","L3"),c("L1","L4"),c("L1","L5"),c("L2","L5"),c("L3","L5"))
alfa_div_comparison_layer_WM_pielou <- list(c("L1","L4"),c("L1","L5"),c("L2","L5"),c("L3","L5"),c("L4","L5"))


###plot using ggviolin   -----> FLIP (coord_flip)
##test kruskall_wallis and wilcox
#to compare layers from WR
violin_plot_psO_jki_seq1_WR_rarefied_diversity_layer <- ggviolin(psO_jki_seq1_WR_filt_rarefied.meta, x="Layer", y= "shannon", ylab = "Shannon index", xlab = "", title = "WR") +
  geom_jitter(aes(color=Layer), width=0.1, size = 3.5) +
  theme_bw() +
  scale_x_discrete(limits=c("L5","L4", "L3", "L2", "L1")) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(1.0, 9.5), breaks = seq(1, 9, by = 4.0), label = c("1.0", "5.0", "9.0")) +
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font(legend.title = element_text (size = 20)) +
  #font(legend.text = element_text (size = 22)) +
  #font("y.text", size = 22) +
  font("xlab", size = 18) +
  coord_flip()+
  scale_color_manual(values =c("#3B0404", "#B95C50", "#DE847B", "#DEB3AD", "#DCC8BA")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 2.5, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_layer_WR_shannon, method = "wilcox", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=5, label.y = c(8.5,8.8,9.1), vjust = 0.8)
violin_plot_psO_jki_seq1_WR_rarefied_diversity_layer

ggsave("violin_plot_psO_jki_seq1_WR_rarefied_diversity_layer.png", path = "~/path/", width = 8, height = 12, units = "cm", dpi = 300, device = "png")

#to compare layers from WM
violin_plot_psO_jki_seq1_WM_rarefied_diversity_layer <- ggviolin(psO_jki_seq1_WM_filt_rarefied.meta, x="Layer", y= "shannon", ylab = "Shannon index", xlab = "", title = "WM") +
  geom_jitter(aes(color=Layer), width=0.1, size = 3.5) +
  theme_bw() +
  scale_x_discrete(limits=c("L5","L4", "L3", "L2", "L1")) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(1.0, 9.5), breaks = seq(1, 9, by = 4.0), label = c("1.0", "5.0", "9.0")) +
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font(legend.title = element_text (size = 20)) +
  #font(legend.text = element_text (size = 22)) +
  #font("y.text", size = 22) +
  font("xlab", size = 18) +
  coord_flip()+
  scale_color_manual(values =c("#3B0404", "#B95C50", "#DE847B", "#DEB3AD", "#DCC8BA")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 2.5, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_layer_WM_shannon, method = "wilcox", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=5, label.y = c(8.5,8.8,9.1), vjust = 0.8)
violin_plot_psO_jki_seq1_WM_rarefied_diversity_layer

ggsave("violin_plot_psO_jki_seq1_WM_rarefied_diversity_layer.png", path = "~/path/", width = 8, height = 12, units = "cm", dpi = 300, device = "png")


violin_plot_psO_jki_seq1_WR_rarefied_richness_layer <- ggviolin(psO_jki_seq1_WR_filt_rarefied.meta, x="Layer", y= "chao1", ylab = "Chao1 index", xlab = "", title = "WR") +
  geom_jitter(aes(color=Layer), width=0.1, size = 3.5) +
  theme_bw() +
  scale_x_discrete(limits=c("L5","L4", "L3", "L2", "L1")) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 4000, by = 2000), label = c("0", "2000", "4000")) +
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font(legend.title = element_text (size = 20)) +
  #font(legend.text = element_text (size = 22)) +
  #font("y.text", size = 22) +
  font("xlab", size = 18) +
  coord_flip()+
  scale_color_manual(values =c("#3B0404", "#B95C50", "#DE847B", "#DEB3AD", "#DCC8BA")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 9.9, size=6) +
  stat_compare_means(comparisons = alfa_div_comparison_layer_WR_chao1, method = "wilcox", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=5, label.y = c(4400, 4550), vjust = 0.8)
violin_plot_psO_jki_seq1_WR_rarefied_richness_layer

ggsave("violin_plot_psO_jki_seq1_WR_rarefied_richness_layer.png", path = "~/path/", width = 8, height = 12, units = "cm", dpi = 300, device = "png")

#to compare layers from WM
violin_plot_psO_jki_seq1_WM_rarefied_richness_layer <- ggviolin(psO_jki_seq1_WM_filt_rarefied.meta, x="Layer", y= "chao1", ylab = "Chao1 index", xlab = "", title = "WM") +
  geom_jitter(aes(color=Layer), width=0.1, size = 3.5) +
  theme_bw() +
  scale_x_discrete(limits=c("L5","L4", "L3", "L2", "L1")) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 4000, by = 2000), label = c("0", "2000", "4000")) +
  theme(axis.text.x = element_text(size=20, angle = 0, vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=20)) +
  #font(legend.title = element_text (size = 20)) +
  #font(legend.text = element_text (size = 22)) +
  #font("y.text", size = 22) +
  font("xlab", size = 18) +
  coord_flip()+
  scale_color_manual(values =c("#3B0404", "#B95C50", "#DE847B", "#DEB3AD", "#DCC8BA")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 9.9, size=6) +
  stat_compare_means(comparisons = alfa_div_comparison_layer_WM_chao1, method = "wilcox", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=5, label.y = c(4000,4150, 4300, 4450, 4600,4750,4900, 5050), vjust = 0.8)
violin_plot_psO_jki_seq1_WM_rarefied_richness_layer 

ggsave("violin_plot_psO_jki_seq1_WM_rarefied_richness_layer.png", path = "~/path/", width = 8, height = 12, units = "cm", dpi = 300, device = "png")







########## Compare microhabitats in rotations
#Select variable of comparison
alfa_div_comparison_microhabitat <- list(c("RA","RP"),c("RA","RH"), c("RH","RP"))
alfa_div_comparison_microhabitat_WR <- list(c("RA","RP"),c("RH","RP"))
alfa_div_comparison_microhabitat_WM <- list(c("RA","RP"),c("RH","RP"))

#to compare microhabitat from WR
violin_plot_psO_jki_seq1_WR_rarefied_diversity_microhabitat <- ggviolin(psO_jki_seq1_WR_filt_rarefied.meta, x="Microhabitat", y= "shannon", ylab = "Shannon index", xlab = "", title = "WR") +
  geom_jitter(aes(color=Microhabitat), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  #theme(axis.text.x = element_text()) +
  #font(legend.title = element_text (size = 20)) +
  #font(legend.text = element_text (size = 22)) +
  font("xy.text", size = 22) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#b05644", "#d9b967", "#57896A")) +
  scale_y_continuous(limits = c(3, 9), breaks = seq(3, 9, by = 3.0), label = c("3.0", "6.0", "9.0")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 8.9, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat_WR, method = "wilcox", label= "p", bracket.size = .6, size=7, label.y = c(8,8.5, 9), vjust = 0.7) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat_WR, method = "wilcox", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=5, label.y = c(8,8.5), vjust = 0.8) 
violin_plot_psO_jki_seq1_WR_rarefied_diversity_microhabitat

ggsave("violin_plot_psO_jki_seq1_WR_rarefied_diversity_microhabitat.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")


violin_plot_psO_jki_seq1_WM_rarefied_diversity_microhabitat <- ggviolin(psO_jki_seq1_WM_filt_rarefied.meta, x="Microhabitat", y= "shannon", ylab = "Shannon index", xlab = "", title = "WM") +
  geom_jitter(aes(color=Microhabitat), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  #theme(axis.text.x = element_text()) +
  #font(legend.title = element_text (size = 20)) +
  #font(legend.text = element_text (size = 22)) +
  font("xy.text", size = 22) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#b05644", "#d9b967", "#57896A")) +
  scale_y_continuous(limits = c(3, 9), breaks = seq(3, 9, by = 3.0), label = c("3.0", "6.0", "9.0")) +
  stat_compare_means(method = "kruskal.test", label= "p", label.y = 8.9, size=5) +
  #stat_compare_means(comparisons = alfa_div_comparison_microhabitat_WM, method = "wilcox", label= "p", bracket.size = .6, size=7, label.y = c(8,8.5, 9), vjust = 0.7) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat_WM, method = "wilcox", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=5, label.y = c(8,8.5), vjust = 0.8) 
violin_plot_psO_jki_seq1_WM_rarefied_diversity_microhabitat

ggsave("violin_plot_psO_jki_seq1_WM_rarefied_diversity_microhabitat.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")



###########################################     For each Microhabitat group

##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_RA_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq1_RH_filt)), step=20, cex=0.5, col = "red")
rarecurve(t(otu_table(psO_jki_seq1_RP_filt)), step=20, cex=0.5, col = "red")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq1_RA_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_RA_filt, rngseed=2022, sample.size = min(sample_sums(psO_jki_seq1_RA_filt)),trimOTUs=TRUE)
psO_jki_seq1_RA_filt_rarefied
sample_sums(psO_jki_seq1_RA_filt_rarefied)

psO_jki_seq1_RH_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_RH_filt, rngseed=2022, sample.size = min(sample_sums(psO_jki_seq1_RH_filt)),trimOTUs=TRUE)
psO_jki_seq1_RH_filt_rarefied
sample_sums(psO_jki_seq1_RH_filt_rarefied)

psO_jki_seq1_RP_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_RP_filt, rngseed=2022, sample.size = min(sample_sums(psO_jki_seq1_RP_filt)),trimOTUs=TRUE)
psO_jki_seq1_RP_filt_rarefied
sample_sums(psO_jki_seq1_RP_filt_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_RA_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_RH_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_RP_filt_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq1_RA_filt_rarefied_1 <- richness(psO_jki_seq1_RA_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_RH_filt_rarefied_1 <- richness(psO_jki_seq1_RH_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_RP_filt_rarefied_1 <- richness(psO_jki_seq1_RP_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq1_RA_filt_rarefied_2 <- evenness(psO_jki_seq1_RA_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_RH_filt_rarefied_2 <- evenness(psO_jki_seq1_RH_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_RP_filt_rarefied_2 <- evenness(psO_jki_seq1_RP_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq1_RA_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_RA_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_RH_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_RH_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_RP_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_RP_filt_rarefied, index = 'shannon', zeroes = TRUE)

##Boxplot alpha diversity using ggviolin function from ggpubr
#Prepare file using meta function
psO_jki_seq1_RA_filt_rarefied.meta <- meta(psO_jki_seq1_RA_filt_rarefied)
psO_jki_seq1_RH_filt_rarefied.meta <- meta(psO_jki_seq1_RH_filt_rarefied)
psO_jki_seq1_RP_filt_rarefied.meta <- meta(psO_jki_seq1_RP_filt_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for RA
psO_jki_seq1_RA_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_RA_filt_rarefied_1$observed
psO_jki_seq1_RA_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_RA_filt_rarefied_1$chao1
psO_jki_seq1_RA_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_RA_filt_rarefied_2$pielou
psO_jki_seq1_RA_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_RA_filt_rarefied_3$shannon

head(psO_jki_seq1_RA_filt_rarefied.meta)
write.csv(psO_jki_seq1_RA_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_RA_filt_rarefied_meta.csv")

#for RH
psO_jki_seq1_RH_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_RH_filt_rarefied_1$observed
psO_jki_seq1_RH_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_RH_filt_rarefied_1$chao1
psO_jki_seq1_RH_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_RH_filt_rarefied_2$pielou
psO_jki_seq1_RH_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_RH_filt_rarefied_3$shannon

head(psO_jki_seq1_RH_filt_rarefied.meta)
write.csv(psO_jki_seq1_RH_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_RH_filt_rarefied_meta.csv")

#for RP
psO_jki_seq1_RP_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_RP_filt_rarefied_1$observed
psO_jki_seq1_RP_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_RP_filt_rarefied_1$chao1
psO_jki_seq1_RP_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_RP_filt_rarefied_2$pielou
psO_jki_seq1_RP_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_RP_filt_rarefied_3$shannon

head(psO_jki_seq1_RP_filt_rarefied.meta)
write.csv(psO_jki_seq1_RP_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_RP_filt_rarefied_meta.csv")

###plot using ggviolin
##test Wilcoxon
#to compare rotations Shannon
#no color, colorful dots by rotation
violin_plot_psO_jki_seq1_RA_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_RA_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "Shannon index", xlab = "", title = "RA") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(2.0, 8.0), breaks = seq(2, 8, by = 2.0), label = c("2.0", "4.0", "6.0", "8.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 7.9, size=6)
violin_plot_psO_jki_seq1_RA_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_RA_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RA_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_RA_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "Chao1 index", xlab = "", title = "RA") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(300, 4300), breaks = seq(300, 4300, by = 1000), label = c("300", "1300", "2300", "3300", "4300")) +
stat_compare_means(method = "wilcox", label= "p", label.y = 4000, size=6)
violin_plot_psO_jki_seq1_RA_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_RA_filt_rarefied_richness_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RA_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_RA_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "Pielou index", xlab = "", title = "RA") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(0.25, 1.00), breaks = seq(0.25, 1.0, by = 0.25), label = c("0.25", "0.50", "0.75", "1.00")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 0.97, size=6)
violin_plot_psO_jki_seq1_RA_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_RA_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RH_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_RH_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "RH") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(2.0, 8.0), breaks = seq(2, 8, by = 2.0), label = c("2.0", "4.0", "6.0", "8.0"))
#stat_compare_means(method = "wilcox", label= "p", label.y = 7.9, size=6)
violin_plot_psO_jki_seq1_RH_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_RH_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RH_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_RH_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "Chao1 index", xlab = "", title = "RH") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(300, 4300), breaks = seq(300, 4300, by = 1000), label = c("300", "1300", "2300", "3300", "4300"))
#stat_compare_means(method = "wilcox", label= "p", label.y = 4000, size=6)
violin_plot_psO_jki_seq1_RH_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_RH_filt_rarefied_richness_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RH_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_RH_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "Pielou index", xlab = "", title = "RH") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(0.25, 1.0), breaks = seq(0.25, 1.0, by = 0.25), label = c("0.25", "0.50", "0.75", "1.00"))
#stat_compare_means(method = "wilcox", label= "p", label.y = 0.97, size=6)
violin_plot_psO_jki_seq1_RH_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_RH_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RP_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_RP_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "RP") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(2.0, 8.0), breaks = seq(2, 8, by = 2.0), label = c("2.0", "4.0", "6.0", "8.0"))
#stat_compare_means(method = "wilcox", label= "p", label.y = 7.9, size=6)
violin_plot_psO_jki_seq1_RP_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_RP_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RP_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_RP_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "Chao1 index", xlab = "", title = "RP") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(300, 4300), breaks = seq(300, 4300, by = 1000), label = c("300", "1300", "2300", "3300", "4300"))
#stat_compare_means(method = "wilcox", label= "p", label.y = 4000, size=6)
violin_plot_psO_jki_seq1_RP_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_RP_filt_rarefied_richness_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_RP_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_RP_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "Pielou index", xlab = "", title = "RP") +
  geom_jitter(aes(color=Rotation), width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  font("y.text", size = 20) +
  font("ylab", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  scale_y_continuous(limits = c(0.25, 1.0), breaks = seq(0.25, 1.0, by = 0.25), label = c("0.25", "0.50", "0.75", "1.00"))
#stat_compare_means(method = "wilcox", label= "p", label.y = 0.97, size=6)
violin_plot_psO_jki_seq1_RP_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_RP_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")




###########################################     For each Layer group

##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_layer1_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq1_layer2_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq1_layer3_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq1_layer4_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq1_layer5_filt)), step=20, cex=0.5, col = "blue")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq1_layer1_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_layer1_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_layer1_filt)),trimOTUs=TRUE)
psO_jki_seq1_layer1_filt_rarefied
sample_sums(psO_jki_seq1_layer1_filt_rarefied)

psO_jki_seq1_layer2_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_layer2_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_layer2_filt)),trimOTUs=TRUE)
psO_jki_seq1_layer2_filt_rarefied
sample_sums(psO_jki_seq1_layer2_filt_rarefied)

psO_jki_seq1_layer3_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_layer3_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_layer3_filt)),trimOTUs=TRUE)
psO_jki_seq1_layer3_filt_rarefied
sample_sums(psO_jki_seq1_layer3_filt_rarefied)

psO_jki_seq1_layer4_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_layer4_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_layer4_filt)),trimOTUs=TRUE)
psO_jki_seq1_layer4_filt_rarefied
sample_sums(psO_jki_seq1_layer4_filt_rarefied)

psO_jki_seq1_layer5_filt_rarefied<-rarefy_even_depth(psO_jki_seq1_layer5_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq1_layer5_filt)),trimOTUs=TRUE)
psO_jki_seq1_layer5_filt_rarefied
sample_sums(psO_jki_seq1_layer5_filt_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq1_layer1_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_layer2_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_layer3_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_layer4_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq1_layer5_filt_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq1_layer1_filt_rarefied_1 <- richness(psO_jki_seq1_layer1_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_layer2_filt_rarefied_1 <- richness(psO_jki_seq1_layer2_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_layer3_filt_rarefied_1 <- richness(psO_jki_seq1_layer3_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_layer4_filt_rarefied_1 <- richness(psO_jki_seq1_layer4_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq1_layer5_filt_rarefied_1 <- richness(psO_jki_seq1_layer5_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq1_layer1_filt_rarefied_2 <- evenness(psO_jki_seq1_layer1_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_layer2_filt_rarefied_2 <- evenness(psO_jki_seq1_layer2_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_layer3_filt_rarefied_2 <- evenness(psO_jki_seq1_layer3_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_layer4_filt_rarefied_2 <- evenness(psO_jki_seq1_layer4_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq1_layer5_filt_rarefied_2 <- evenness(psO_jki_seq1_layer5_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq1_layer1_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_layer1_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_layer2_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_layer2_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_layer3_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_layer3_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_layer4_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_layer4_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq1_layer5_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq1_layer5_filt_rarefied, index = 'shannon', zeroes = TRUE)

##Boxplot alpha diversity using ggviolin function from ggpubr
#Prepare file using meta function
psO_jki_seq1_layer1_filt_rarefied.meta <- meta(psO_jki_seq1_layer1_filt_rarefied)
psO_jki_seq1_layer2_filt_rarefied.meta <- meta(psO_jki_seq1_layer2_filt_rarefied)
psO_jki_seq1_layer3_filt_rarefied.meta <- meta(psO_jki_seq1_layer3_filt_rarefied)
psO_jki_seq1_layer4_filt_rarefied.meta <- meta(psO_jki_seq1_layer4_filt_rarefied)
psO_jki_seq1_layer5_filt_rarefied.meta <- meta(psO_jki_seq1_layer5_filt_rarefied)

#select column from output file to be plotted and tested (only one by time)
psO_jki_seq1_layer1_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_layer1_filt_rarefied_1$observed
psO_jki_seq1_layer1_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_layer1_filt_rarefied_1$chao1
psO_jki_seq1_layer1_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_layer1_filt_rarefied_2$pielou
psO_jki_seq1_layer1_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_layer1_filt_rarefied_3$shannon

head(psO_jki_seq1_layer1_filt_rarefied.meta)
write.csv(psO_jki_seq1_layer1_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_layer1_filt_rarefied_meta.csv")

psO_jki_seq1_layer2_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_layer2_filt_rarefied_1$observed
psO_jki_seq1_layer2_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_layer2_filt_rarefied_1$chao1
psO_jki_seq1_layer2_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_layer2_filt_rarefied_2$pielou
psO_jki_seq1_layer2_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_layer2_filt_rarefied_3$shannon

head(psO_jki_seq1_layer2_filt_rarefied.meta)
write.csv(psO_jki_seq1_layer2_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_layer2_filt_rarefied_meta.csv")

psO_jki_seq1_layer3_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_layer3_filt_rarefied_1$observed
psO_jki_seq1_layer3_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_layer3_filt_rarefied_1$chao1
psO_jki_seq1_layer3_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_layer3_filt_rarefied_2$pielou
psO_jki_seq1_layer3_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_layer3_filt_rarefied_3$shannon

head(psO_jki_seq1_layer3_filt_rarefied.meta)
write.csv(psO_jki_seq1_layer3_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_layer3_filt_rarefied_meta.csv")

psO_jki_seq1_layer4_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_layer4_filt_rarefied_1$observed
psO_jki_seq1_layer4_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_layer4_filt_rarefied_1$chao1
psO_jki_seq1_layer4_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_layer4_filt_rarefied_2$pielou
psO_jki_seq1_layer4_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_layer4_filt_rarefied_3$shannon

head(psO_jki_seq1_layer4_filt_rarefied.meta)
write.csv(psO_jki_seq1_layer4_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_layer4_filt_rarefied_meta.csv")

psO_jki_seq1_layer5_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq1_layer5_filt_rarefied_1$observed
psO_jki_seq1_layer5_filt_rarefied.meta$chao1 <- alfa_div_psO_jki_seq1_layer5_filt_rarefied_1$chao1
psO_jki_seq1_layer5_filt_rarefied.meta$pielou <- alfa_div_psO_jki_seq1_layer5_filt_rarefied_2$pielou
psO_jki_seq1_layer5_filt_rarefied.meta$shannon <- alfa_div_psO_jki_seq1_layer5_filt_rarefied_3$shannon

head(psO_jki_seq1_layer5_filt_rarefied.meta)
write.csv(psO_jki_seq1_layer5_filt_rarefied.meta, "~/path/alpha_div_psO_jki_seq1_layer5_filt_rarefied_meta.csv")


###plot using ggviolin
##test Wilcoxon test
#to compare rotations Shannon
#no color, colorful dots by rotation
violin_plot_psO_jki_seq1_layer1_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_layer1_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "L1") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer1_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_layer1_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer2_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_layer2_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "L2") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer2_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_layer2_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer3_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_layer3_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "L3") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer3_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_layer3_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer4_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_layer4_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "L4") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer4_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_layer4_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png" )

violin_plot_psO_jki_seq1_layer5_filt_rarefied_diversity_rotation <- ggviolin(psO_jki_seq1_layer5_filt_rarefied.meta, x="Rotation", y= "shannon", ylab = "", xlab = "", title = "L5") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer5_filt_rarefied_diversity_rotation

ggsave("violin_plot_psO_jki_seq1_layer5_filt_rarefied_diversity_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")



#####Richness

violin_plot_psO_jki_seq1_layer1_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_layer1_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "", xlab = "", title = "L1") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer1_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_layer1_filt_rarefied_richness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer2_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_layer2_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "", xlab = "", title = "L2") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer2_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_layer2_filt_rarefied_richness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer3_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_layer3_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "", xlab = "", title = "L3") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer3_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_layer3_filt_rarefied_richness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer4_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_layer4_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "", xlab = "", title = "L4") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer4_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_layer4_filt_rarefied_richness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png" )

violin_plot_psO_jki_seq1_layer5_filt_rarefied_richness_rotation <- ggviolin(psO_jki_seq1_layer5_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "", xlab = "", title = "L5") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer5_filt_rarefied_richness_rotation

ggsave("violin_plot_psO_jki_seq1_layer5_filt_rarefied_richness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")



#####evenness

violin_plot_psO_jki_seq1_layer1_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_layer1_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "", xlab = "", title = "L1") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer1_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_layer1_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer2_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_layer2_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "", xlab = "", title = "L2") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer2_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_layer2_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer3_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_layer3_filt_rarefied.meta, x="Rotation", y= "chao1", ylab = "", xlab = "", title = "L3") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer3_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_layer3_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

violin_plot_psO_jki_seq1_layer4_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_layer4_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "", xlab = "", title = "L4") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer4_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_layer4_filt_rarefied_evenness_rotation.png", path = "~/path/", width = 11, height = 12, units = "cm", dpi = 300, device = "png" )

violin_plot_psO_jki_seq1_layer5_filt_rarefied_evenness_rotation <- ggviolin(psO_jki_seq1_layer5_filt_rarefied.meta, x="Rotation", y= "pielou", ylab = "", xlab = "", title = "L5") +
  geom_jitter(aes(color=Rotation),width=0.1, size = 3.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  #font("legend.text", size = 22) + 
  #font("legend.title", size = 20) +
  #font("ylab", size = 20) +
  #font("y.text", size = 20) +
  scale_color_manual(values =c("#B1AD31", "#196418")) +
  #scale_y_continuous(limits = c(4.0, 7.0), breaks = seq(4, 7, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0")) +
  stat_compare_means(method = "wilcox", label= "p", label.y = 6.9, size=6)
violin_plot_psO_jki_seq1_layer5_filt_rarefied_evenness_rotation

ggsave("violin_plot_psO_jki_seq1_layer5_filt_rarefied_evenness_rotation.png", path = "~/", width = 11, height = 12, units = "cm", dpi = 300, device = "png")

