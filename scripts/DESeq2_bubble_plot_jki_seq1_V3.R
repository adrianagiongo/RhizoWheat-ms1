## DESeq2 - bubble plot - jki_seq1
#### Adriana Giongo
#### (10.07.2023) 

#Load packages
library("ggplot2")
library("reshape2")

############
## DESeq results - by microhabitat and year 
############

## RA
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
deseq2_jki_seq1_RA = read.csv("~/Documents/R_analysis/jki_seq1/data_jki_seq1/Bubble_plot_data/input_deseq2_jki_seq1_RA_layer_average_01_V3.csv", header = TRUE)
deseq2_jki_seq1_RA

#convert data frame from a "wide" format to a "long" format
deseq2_jki_seq1_RA_melt <- melt(deseq2_jki_seq1_RA, id = c("Sample", "Field"))
deseq2_jki_seq1_RA_melt

colours = c("#B1AD31", "#196418")

deseq2_jki_seq1_RA_melt$Sample <- factor(deseq2_jki_seq1_RA_melt$Sample,levels=unique(deseq2_jki_seq1_RA$Sample))

bubbleplot_deseq2_jki_seq1_RA_melt = ggplot(deseq2_jki_seq1_RA_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Field), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.00001, 50), range = c(1,50), breaks = c(0.1, 1, 3, 6)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Field")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 16, colour ="black"), 
        legend.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.8), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)), breaks=c('WR', 'WM'))
bubbleplot_deseq2_jki_seq1_RA_melt

ggsave("bubbleplot_deseq2_jki_seq1_RA_melt.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Bubble_plot_jki_seq1/", width = 26, height = 22, units = "cm", dpi = 300, device = "png")


## RH
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
deseq2_jki_seq1_RH = read.csv("~/Documents/R_analysis/jki_seq1/data_jki_seq1/Bubble_plot_data/input_deseq2_jki_seq1_RH_layer_average_01_V3.csv", header = TRUE)
deseq2_jki_seq1_RH

#convert data frame from a "wide" format to a "long" format
deseq2_jki_seq1_RH_melt <- melt(deseq2_jki_seq1_RH, id = c("Sample", "Field"))
deseq2_jki_seq1_RH_melt

colours = c("#B1AD31", "#196418")

deseq2_jki_seq1_RH_melt$Sample <- factor(deseq2_jki_seq1_RH_melt$Sample,levels=unique(deseq2_jki_seq1_RH$Sample))

bubbleplot_deseq2_jki_seq1_RH_melt = ggplot(deseq2_jki_seq1_RH_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Field), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.00001, 50), range = c(1,50), breaks = c(0.1, 1, 3, 6)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Field")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 0, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 16, colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.6), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)), breaks=c('WR', 'WM'))
bubbleplot_deseq2_jki_seq1_RH_melt

ggsave("bubbleplot_deseq2_jki_seq1_RH_melt.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Bubble_plot_jki_seq1/", width = 26, height = 22, units = "cm", dpi = 300, device = "png")


## RP
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
deseq2_jki_seq1_RP = read.csv("~/Documents/R_analysis/jki_seq1/data_jki_seq1/Bubble_plot_data/input_deseq2_jki_seq1_RP_layer_average_01_V3.csv", header = TRUE)
deseq2_jki_seq1_RP

#convert data frame from a "wide" format to a "long" format
deseq2_jki_seq1_RP_melt <- melt(deseq2_jki_seq1_RP, id = c("Sample", "Field"))
deseq2_jki_seq1_RP_melt

colours = c("#B1AD31", "#196418")

deseq2_jki_seq1_RP_melt$Sample <- factor(deseq2_jki_seq1_RP_melt$Sample,levels=unique(deseq2_jki_seq1_RP$Sample))

bubbleplot_deseq2_jki_seq1_RP_melt = ggplot(deseq2_jki_seq1_RP_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Field), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.00001, 50), range = c(1,50), breaks = c(0.1, 1, 3, 6)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Field")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.9), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)), breaks=c('WR', 'WM'))
bubbleplot_deseq2_jki_seq1_RP_melt

ggsave("bubbleplot_deseq2_jki_seq1_RP_melt.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Bubble_plot_jki_seq1/", width = 25, height = 50, units = "cm", dpi = 300, device = "png")


### The end! Have fun : )