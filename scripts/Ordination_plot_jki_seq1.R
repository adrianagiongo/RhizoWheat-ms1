##Creating MDS plot for selected dataset
#Palette for Rotation  WR = "#B1AD31", WM = "#196418"
#Palette for Microhabitat  RA = "#b05644", RH = "#d9b967", RP = "#57896A"
#Palette for Layer 1 = "#dcc8ba", Layer 2 = "#DEB3AD", Layer 3 = "#DE847B", Layer 4 = "#B95C50", Layer 5 = "#3B0404"

#load package
library("phyloseq")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")

#set seed
set.seed(2022)

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
MDS_bray_psO_jki_seq1_sqr<-ordinate(psO_jki_seq1, "MDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_jki_seq1_sqr)

##Create a MDS plot (all ASV)
#Rotation
plot_MDS_bray_psO_jki_seq1_sqr_rotation<-plot_ordination(psO_jki_seq1,MDS_bray_psO_jki_seq1_sqr, type="sample",color="Rotation", shape = "Microhabitat") +
  theme_bw() +
  geom_point(size=5) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=16)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#196418", "#B1AD31"))
plot_MDS_bray_psO_jki_seq1_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq1_sqr_rotation.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 18, height = 12, units = "cm", dpi = 300, device = "png")

####################################### STATISTICS
#All
##Multivariate Anova for ordinate files
set.seed(2022)

#Transform abundances to square root
psO_jki_seq1_sqr <- microbiome::transform(psO_jki_seq1, "hellinger")

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq1_sqr_abundances <- abundances(psO_jki_seq1_sqr)
psO_jki_seq1_meta <- meta(psO_jki_seq1)

##Permanova using adonis function from vegan and print p value
#for ASV
#by null (to test the whole Model - overall)
set.seed(2022)
Permanova_psO_jki_seq1_rotation_model <- adonis2(t(psO_jki_seq1_sqr_abundances) ~Rotation + Microhabitat + Layer, data = psO_jki_seq1_meta, permutations = 10000, method = "bray", by = NULL)
Permanova_psO_jki_seq1_rotation_model
write.csv(Permanova_psO_jki_seq1_rotation_model, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_rotation_model.csv")

#by terms (it is default)
set.seed(2022)
Permanova_psO_jki_seq1_rotation <- adonis2(t(psO_jki_seq1_sqr_abundances) ~ Rotation, data = psO_jki_seq1_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq1_rotation
write.csv(Permanova_psO_jki_seq1_rotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_rotation.csv")

Permanova_psO_jki_seq1_rotation_byterm_all <- adonis2(t(psO_jki_seq1_sqr_abundances) ~Rotation + Microhabitat + Layer, data = psO_jki_seq1_meta, permutations = 10000, method = "bray", by =  "terms")
Permanova_psO_jki_seq1_rotation_byterm_all
write.csv(Permanova_psO_jki_seq1_rotation_byterm_all, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_rotation_byterm_all.csv")

#by margin
set.seed(2022)
Permanova_psO_jki_seq1_rotation_bymargin <- adonis2(t(psO_jki_seq1_sqr_abundances) ~Rotation:Microhabitat:Layer, data = psO_jki_seq1_meta, permutations = 10000, method = "bray", by =  "margin")
Permanova_psO_jki_seq1_rotation_bymargin
write.csv(Permanova_psO_jki_seq1_rotation_bymargin, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_rotation_bymargin.csv")


####################################### For each Rotation group

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
MDS_bray_psO_jki_seq1_WR_filt_sqr<-ordinate(psO_jki_seq1_WR_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq1_WM_filt_sqr<-ordinate(psO_jki_seq1_WM_filt, "MDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_jki_seq1_WR_filt_sqr)
head(MDS_bray_psO_jki_seq1_WM_filt_sqr)

##Create a MDS plot
#Layer
plot_MDS_bray_psO_jki_seq1_WR_filt_sqr_layer<-plot_ordination(psO_jki_seq1_WR_filt,MDS_bray_psO_jki_seq1_WR_filt_sqr, type="sample",color="Layer") +
  theme_bw() +
  geom_point(size=5) +
  stat_ellipse(type = "norm", linetype = 2, geom = "polygon", alpha=0.2) +
  #stat_ellipse(type = "t") +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=16)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#3B0404", "#B95C50", "#DE847B", "#DEB3AD", "#dcc8ba"))
plot_MDS_bray_psO_jki_seq1_WR_filt_sqr_layer

ggsave("plot_MDS_bray_psO_jki_seq1_WR_filt_sqr_layer.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 13, height = 10, units = "cm",dpi = 300)

#Microhabitat
plot_MDS_bray_psO_jki_seq1_WR_filt_sqr_microhabitat<-plot_ordination(psO_jki_seq1_WR_filt,MDS_bray_psO_jki_seq1_WR_filt_sqr, type="sample",color="Microhabitat") +
  theme_bw() +
  geom_point(size=5) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  #theme(legend.position = c(0.9, 0.875), legend.background = element_rect(fill = "white", color = "gray")) +
  theme(legend.text = element_text(size=14)) +
  scale_color_manual(values = c("#b05644","#d9b967","#57896A"))
plot_MDS_bray_psO_jki_seq1_WR_filt_sqr_microhabitat

ggsave("plot_MDS_bray_psO_jki_seq1_WR_filt_sqr_microhabitat.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 14, height = 11, units = "cm",dpi = 300)

#Layer
plot_MDS_bray_psO_jki_seq1_WM_filt_sqr_layer<-plot_ordination(psO_jki_seq1_WM_filt,MDS_bray_psO_jki_seq1_WM_filt_sqr, type="sample",color="Layer") +
  theme_bw() +
  geom_point(size=5) +
  stat_ellipse(type = "norm", linetype = 2, geom = "polygon", alpha=0.2) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=16)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#3B0404", "#B95C50", "#DE847B", "#DEB3AD", "#dcc8ba"))
plot_MDS_bray_psO_jki_seq1_WM_filt_sqr_layer

ggsave("plot_MDS_bray_psO_jki_seq1_WM_filt_sqr_layer.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 13, height = 10, units = "cm",dpi = 300)

#Microhabitat
plot_MDS_bray_psO_jki_seq1_WM_filt_sqr_microhabitat<-plot_ordination(psO_jki_seq1_WM_filt,MDS_bray_psO_jki_seq1_WM_filt_sqr, type="sample",color="Microhabitat") +
  theme_bw() +
  geom_point(size=5) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=16)) +
  theme(legend.title = element_blank()) +
  #theme(legend.position = c(0.9, 0.875), legend.background = element_rect(fill = "white", color = "gray")) +
  theme(legend.text = element_text(size=14)) +
  scale_color_manual(values = c("#b05644","#d9b967","#57896A"))
plot_MDS_bray_psO_jki_seq1_WM_filt_sqr_microhabitat

ggsave("plot_MDS_bray_psO_jki_seq1_WM_filt_sqr_microhabitat.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 14, height = 11, units = "cm",dpi = 300)

####################################### STATISTICS
#Rotation
##Multivariate Anova for ordinate files
set.seed(2022)

#Transform abundances to square root
psO_jki_seq1_WR_filt_sqr <- microbiome::transform(psO_jki_seq1_WR_filt, "hellinger")
psO_jki_seq1_WM_filt_sqr <- microbiome::transform(psO_jki_seq1_WM_filt, "hellinger")

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq1_WR_filt_sqr_abundances <- abundances(psO_jki_seq1_WR_filt_sqr)
psO_jki_seq1_WR_filt_meta <- meta(psO_jki_seq1_WR_filt)

psO_jki_seq1_WM_filt_sqr_abundances <- abundances(psO_jki_seq1_WM_filt_sqr)
psO_jki_seq1_WM_filt_meta <- meta(psO_jki_seq1_WM_filt)

##Permanova using adonis function from vegan and print p value
#WR data
Permanova_psO_jki_seq1_WR_filt_microhabitat <- adonis2(t(psO_jki_seq1_WR_filt_sqr_abundances) ~ Microhabitat + Layer, data = psO_jki_seq1_WR_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq1_WR_filt_microhabitat
write.csv(Permanova_psO_jki_seq1_WR_filt_microhabitat, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_WR_filt_microhabitat.csv")

##WM data
Permanova_psO_jki_seq1_WM_filt_microhabitat <- adonis2(t(psO_jki_seq1_WM_filt_sqr_abundances) ~ Microhabitat + Layer, data = psO_jki_seq1_WM_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq1_WM_filt_microhabitat
write.csv(Permanova_psO_jki_seq1_WM_filt_microhabitat, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_WM_filt_microhabitat.csv")


####################################### For each Microhabitat group


#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
MDS_bray_psO_jki_seq1_RA_filt_sqr<-ordinate(psO_jki_seq1_RA_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq1_RH_filt_sqr<-ordinate(psO_jki_seq1_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq1_RP_filt_sqr<-ordinate(psO_jki_seq1_RP_filt, "MDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_jki_seq1_RA_filt_sqr)
head(MDS_bray_psO_jki_seq1_RH_filt_sqr)
head(MDS_bray_psO_jki_seq1_RP_filt_sqr)

##Create a MDS plot
#Rotation
plot_MDS_bray_psO_jki_seq1_RA_filt_sqr_rotation<-plot_ordination(psO_jki_seq1_RA_filt,MDS_bray_psO_jki_seq1_RA_filt_sqr, type="sample", color="Rotation") +
  theme_bw() +
  geom_point(size=5) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.position = "none") +
  #theme(legend.title = element_text(size=20)) +
  #theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#196418", "#B1AD31"))
plot_MDS_bray_psO_jki_seq1_RA_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq1_RA_filt_sqr_rotation.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 12, height = 12, units = "cm", dpi = 300, device = "png")

plot_MDS_bray_psO_jki_seq1_RH_filt_sqr_rotation<-plot_ordination(psO_jki_seq1_RH_filt,MDS_bray_psO_jki_seq1_RH_filt_sqr, type="sample",color="Rotation") +
  theme_bw() +
  geom_point(size=5) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.position = "none") +
  #theme(legend.title = element_text(size=20)) +
  #theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#196418", "#B1AD31"))
plot_MDS_bray_psO_jki_seq1_RH_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq1_RH_filt_sqr_rotation.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 10, height = 10, units = "cm", dpi = 300, device = "png")

plot_MDS_bray_psO_jki_seq1_RP_filt_sqr_rotation<-plot_ordination(psO_jki_seq1_RP_filt,MDS_bray_psO_jki_seq1_RP_filt_sqr, type="sample",color="Rotation") +
  theme_bw() +
  geom_point(size=5) +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.position = "none") +
  #theme(legend.title = element_text(size=20)) +
  #theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#196418", "#B1AD31"))
plot_MDS_bray_psO_jki_seq1_RP_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq1_RP_filt_sqr_rotation.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Ordination_jki_seq1/", width = 10, height = 10, units = "cm", dpi = 300, device = "png")

####################################### STATISTICS
#Microhabitat
##Multivariate Anova for ordinate files
set.seed(2022)

#Transform abundances to square root
psO_jki_seq1_RA_filt_sqr <- microbiome::transform(psO_jki_seq1_RA_filt, "hellinger")
psO_jki_seq1_RH_filt_sqr <- microbiome::transform(psO_jki_seq1_RH_filt, "hellinger")
psO_jki_seq1_RP_filt_sqr <- microbiome::transform(psO_jki_seq1_RP_filt, "hellinger")

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq1_RA_filt_sqr_abundances <- abundances(psO_jki_seq1_RA_filt_sqr)
psO_jki_seq1_RA_filt_meta <- meta(psO_jki_seq1_RA_filt)

psO_jki_seq1_RH_filt_sqr_abundances <- abundances(psO_jki_seq1_RH_filt_sqr)
psO_jki_seq1_RH_filt_meta <- meta(psO_jki_seq1_RH_filt)

psO_jki_seq1_RP_filt_sqr_abundances <- abundances(psO_jki_seq1_RP_filt_sqr)
psO_jki_seq1_RP_filt_meta <- meta(psO_jki_seq1_RP_filt)

##Permanova using adonis function from vegan and print p value
Permanova_psO_jki_seq1_RA_filt_rotation <- adonis2(t(psO_jki_seq1_RA_filt_sqr_abundances) ~ Rotation + Layer, data = psO_jki_seq1_RA_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq1_RA_filt_rotation
write.csv(Permanova_psO_jki_seq1_RA_filt_rotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RA_filt_rotation.csv")

Permanova_psO_jki_seq1_RH_filt_rotation <- adonis2(t(psO_jki_seq1_RH_filt_sqr_abundances) ~ Rotation + Layer, data = psO_jki_seq1_RH_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq1_RH_filt_rotation
write.csv(Permanova_psO_jki_seq1_RH_filt_rotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RH_filt_rotation.csv")

Permanova_psO_jki_seq1_RP_filt_rotation <- adonis2(t(psO_jki_seq1_RP_filt_sqr_abundances) ~ Rotation + Layer, data = psO_jki_seq1_RP_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq1_RP_filt_rotation
write.csv(Permanova_psO_jki_seq1_RP_filt_rotation, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RP_filt_rotation.csv")


##Permanova using adonis function from vegan and print p value
#for ASV
#by null (to test the whole Model - overall)
set.seed(2022)
Permanova_psO_jki_seq1_RA_filt_model <- adonis2(t(psO_jki_seq1_RA_filt_sqr_abundances) ~Rotation + Layer, data = psO_jki_seq1_RA_filt_meta, permutations = 10000, method = "bray", by = NULL)
Permanova_psO_jki_seq1_RA_filt_model
write.csv(Permanova_psO_jki_seq1_RA_filt_model, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RA_filt_model.csv")

set.seed(2022)
Permanova_psO_jki_seq1_RH_filt_model <- adonis2(t(psO_jki_seq1_RH_filt_sqr_abundances) ~Rotation + Layer, data = psO_jki_seq1_RH_filt_meta, permutations = 10000, method = "bray", by = NULL)
Permanova_psO_jki_seq1_RH_filt_model
write.csv(Permanova_psO_jki_seq1_RH_filt_model, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RH_filt_model.csv")

set.seed(2022)
Permanova_psO_jki_seq1_RP_filt_model <- adonis2(t(psO_jki_seq1_RP_filt_sqr_abundances) ~Rotation + Layer, data = psO_jki_seq1_RP_filt_meta, permutations = 10000, method = "bray", by = NULL)
Permanova_psO_jki_seq1_RP_filt_model
write.csv(Permanova_psO_jki_seq1_RP_filt_model, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RP_filt_model.csv")

#by terms (it is default)
set.seed(2022)
Permanova_psO_jki_seq1_RA_filt_byterm <- adonis2(t(psO_jki_seq1_RA_filt_sqr_abundances) ~Rotation + Layer, data = psO_jki_seq1_RA_filt_meta, permutations = 10000, method = "bray", by = "terms")
Permanova_psO_jki_seq1_RA_filt_byterm
write.csv(Permanova_psO_jki_seq1_RA_filt_byterm, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RA_filt_byterm.csv")

set.seed(2022)
Permanova_psO_jki_seq1_RH_filt_byterm <- adonis2(t(psO_jki_seq1_RH_filt_sqr_abundances) ~Rotation + Layer, data = psO_jki_seq1_RH_filt_meta, permutations = 10000, method = "bray", by = "terms")
Permanova_psO_jki_seq1_RH_filt_byterm
write.csv(Permanova_psO_jki_seq1_RH_filt_byterm, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RH_filt_byterm.csv")

set.seed(2022)
Permanova_psO_jki_seq1_RP_filt_byterm <- adonis2(t(psO_jki_seq1_RP_filt_sqr_abundances) ~Rotation + Layer, data = psO_jki_seq1_RP_filt_meta, permutations = 10000, method = "bray", by = "terms")
Permanova_psO_jki_seq1_RP_filt_byterm
write.csv(Permanova_psO_jki_seq1_RP_filt_byterm, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RP_filt_byterm.csv")

#by margin
set.seed(2022)
Permanova_psO_jki_seq1_RA_filt_bymargin <- adonis2(t(psO_jki_seq1_RA_filt_sqr_abundances) ~Rotation : Layer, data = psO_jki_seq1_RA_filt_meta, permutations = 10000, method = "bray", by = "margin")
Permanova_psO_jki_seq1_RA_filt_bymargin
write.csv(Permanova_psO_jki_seq1_RA_filt_bymargin, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RA_filt_bymargin.csv")

set.seed(2022)
Permanova_psO_jki_seq1_RH_filt_bymargin <- adonis2(t(psO_jki_seq1_RH_filt_sqr_abundances) ~Rotation : Layer, data = psO_jki_seq1_RH_filt_meta, permutations = 10000, method = "bray",  by = "margin")
Permanova_psO_jki_seq1_RH_filt_bymargin
write.csv(Permanova_psO_jki_seq1_RH_filt_bymargin, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RH_filt_bymargin.csv")

set.seed(2022)
Permanova_psO_jki_seq1_RP_filt_bymargin <- adonis2(t(psO_jki_seq1_RP_filt_sqr_abundances) ~Rotation : Layer, data = psO_jki_seq1_RP_filt_meta, permutations = 10000, method = "bray",  by = "margin")
Permanova_psO_jki_seq1_RP_filt_bymargin
write.csv(Permanova_psO_jki_seq1_RP_filt_bymargin, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/Permanova_psO_jki_seq1_RP_filt_bymargin.csv")


#### ANOSIM
## testing of significance for the bray-curtis dissimilarity using ANOSIM
#calculate distance values between samples
dist_psO_jki_seq1 = phyloseq::distance(psO_jki_seq1, method = "bray")

#create a dataframe of the metadata
metadata_psO_jki_seq1 <- data.frame(sample_data(psO_jki_seq1))

#run Anosim using 10000 permutations
anosim_dist_rotation_psO_jki_seq1 <- anosim(dist_psO_jki_seq1, metadata_psO_jki_seq1$Rotation, permutations = 10000)
anosim_dist_microhabitat_psO_jki_seq1 <- anosim(dist_psO_jki_seq1, metadata_psO_jki_seq1$Microhabitat, permutations = 10000)
anosim_dist_layer_psO_jki_seq1 <- anosim(dist_psO_jki_seq1, metadata_psO_jki_seq1$Layer, permutations = 10000)

#print results
print(anosim_dist_rotation_psO_jki_seq1)
print(anosim_dist_microhabitat_psO_jki_seq1)
print(anosim_dist_layer_psO_jki_seq1)

####################################### PAIRWISE ANOSIM   ---> Microhabitat

## testing of pairwise significance for the bray-curtis dissimilarity using ANOSIM
#Create metadata with an unique variable and run Anosim using 10000 permutations
metadata_microhabitat_psO_jki_seq1 <- combn(x=unique(metadata_psO_jki_seq1$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq1 <- c()

for (i in 1:ncol(metadata_microhabitat_psO_jki_seq1)){
  ps_subs_psO_jki_seq1 <- subset_samples(psO_jki_seq1, Microhabitat %in% metadata_microhabitat_psO_jki_seq1 [,i])
  metadata_subs_psO_jki_seq1 <- data.frame(sample_data(ps_subs_psO_jki_seq1))
  pairwise_anosim_bray_dist_psO_jki_seq1 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq1, method= "bray"), metadata_subs_psO_jki_seq1$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq1 <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq1, pairwise_anosim_bray_dist_psO_jki_seq1$signif[1])
}

#Adjust statistics values using BH
p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq1 <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq1, method = "BH")

#Create a table with statistics results
p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq1 <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq1), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq1, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq1)

#print results
#print(pairwise_anosim_dist_psO_jki_seq1)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq1)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq1)


