#Create a venndiagram 
#Palette for Rotation  W1 = "#DED93E", W2 = "#8BCD50", WM = "#196418"
#Palette for Microhabitat  RA = "#b05644", RH = "#d9b967", RP = "#57896A"
#Palette for Depth 0-15 = "#dcc8ba", 15-30 = "#DEB3AD", 30-60 = "#DE847B", 60-90 = "#B95C50", 90-120 = "#3B0404"

#loading packages
library("phyloseq")
library("limma")
library("VennDiagram")
library("ggplot2")

###############
#Venndiagram Microhabitat

#verify name of variables
sample_variables(psO_jki_seq1_WR_filt_annotation)
sample_variables(psO_jki_seq1_WM_filt_annotation)

#merge samples by time
venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat <- merge_samples(psO_jki_seq1_WR_filt_annotation, "Microhabitat")
sample_data(venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat)
venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat

venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat <- merge_samples(psO_jki_seq1_WM_filt_annotation, "Microhabitat")
sample_data(venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat)
venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat

#create the object to calculate the variable intersections
t_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat <- t(otu_table(venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat))
t_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat <- t(otu_table(venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat))

#formatting results on table and save
df_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat <- data.frame(tax_table(venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat), otu_table(t_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat))
head(df_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat)
dim(df_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat)
write.csv(df_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat.csv")

df_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat <- data.frame(tax_table(venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat), otu_table(t_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat))
head(df_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat)
dim(df_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat)
write.csv(df_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat.csv")

#calculate the variable intersections
venn_counts_jki_seq1_WR_filt_annotation_microhabitat <- vennCounts(t_venndiagram_psO_jki_seq1_WR_filt_annotation_microhabitat)
venn_counts_jki_seq1_WR_filt_annotation_microhabitat

venn_counts_jki_seq1_WM_filt_annotation_microhabitat <- vennCounts(t_venndiagram_psO_jki_seq1_WM_filt_annotation_microhabitat)
venn_counts_jki_seq1_WM_filt_annotation_microhabitat


#plot Venndiagram WR
venn_plot_jki_seq1_WR_filt_annotation_microhabitat <- draw.triple.venn(
  area1 = 1390,
  area2 = 1363,
  area3 = 1300,
  n12 = 1128,
  n13 = 1075,
  n23 = 1043,
  n123 = 966,
  category = c("RA", "RH", "RP"),
  fill = c("#DC7726", "#d9b967", "#57896A"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#DC7726", "#d9b967", "#57896A")
);

tiff(filename = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Venndiagram_jki_seq1/venn_plot_jki_seq1_WR_filt_annotation_microhabitat.png", width = 15, height = 15, units = "cm", res = 300);
grid.draw(venn_plot_jki_seq1_WR_filt_annotation_microhabitat);
dev.off()


#plot Venndiagram WM
venn_plot_jki_seq1_WM_filt_annotation_microhabitat <- draw.triple.venn(
  area1 = 1265,
  area2 = 1390,
  area3 = 1244,
  n12 = 1078,
  n13 = 1019,
  n23 = 1055,
  n123 = 948,
  category = c("RA", "RH", "RP"),
  fill = c("#DC7726", "#d9b967", "#57896A"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#DC7726", "#d9b967", "#57896A")
);

tiff(filename = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Venndiagram_jki_seq1/venn_plot_jki_seq1_WM_filt_annotation_microhabitat.png", width = 15, height = 15, units = "cm", res = 300);
grid.draw(venn_plot_jki_seq1_WM_filt_annotation_microhabitat);
dev.off()

