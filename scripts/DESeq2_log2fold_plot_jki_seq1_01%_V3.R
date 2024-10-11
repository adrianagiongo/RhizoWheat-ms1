##Define taxa deferentially abundant using DESeq2 on phyloseq pipeline

#Loading package
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")

##### This script uses a file from the output that was modified to contain only what presented more than 0.1% (relative abundance) 
sigtab_RA <- read.csv ("~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/sigtab_RA_filt_annotation_01_V3.csv", row.names = 1)
sigtab_RH <- read.csv ("~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/sigtab_RH_filt_annotation_01_V3.csv", row.names = 1)
sigtab_RP <- read.csv ("~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/sigtab_RP_filt_annotation_01_V3.csv", row.names = 1)

# colors_RP <- list(
#   Rotation = c("WR" = "#B1AD31", "WM" = "#196418"),
#   Layer = c("L1" = "#3B0404", "L2" = "#B95C50", "L3" = "#DE847B", "L4" = "#DEB3AD", "L5" = "#E9DDD4"),
#   WRvsWM = c("down" = "#DB1F48", "up" = "#0074B7"),
#   Phylum = c("Acidobacteriota" = "#ffc331", "Actinobacteriota" = "#ff8b60", "Bacteroidota"  = "#Bff4be", "Bdellovibrionota" = "#185113", 
#              "Chloroflexi" = "#2b8f22", "Entotheonellaeota" = "#d4ff31", "Firmicutes" = "#979191", "Gemmatimonadota" = "#363781",
#              "Latescibacterota" = "#4a2500", "Methylomirabilota" = "#964B00", "Myxococcota" = "#FFC0CB", "NB1-j" = "#b5814c", "Nitrospirota" = "#1361cf",
#              "Patescibacteria" = "#4b0096", "Planctomycetota" = "#9fc5e8", "Proteobacteria" = "#ca250f", "RCP2-54" = "#800080",
#              "Verrucomicrobiota" = "#ced2ce", "Zixibacteria" = "#e89fc5"))


##Convert variables characters in factors ---> Done with DESeq2_jki_seq1_like24
# Phylum order
x = tapply(sigtab_RA$log2FoldChange, sigtab_RA$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RA$Phylum = factor(as.character(sigtab_RA$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RA$log2FoldChange, sigtab_RA$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RA$Annotation = factor(as.character(sigtab_RA$Annotation), levels=names(x))

# Define color palette
# Define color palette
color.phylum_RA <- c("Actinobacteriota" = "#ff8b60", 
                        "Bacteroidota"  = "#Bff4be", 
                        "Calditrichota" = "#979191", 
                        "Chloroflexi" = "#2b8f22", 
                        "Firmicutes" = "#b300b3", 
                        "Methylomirabilota" = "#964B00", 
                        "Proteobacteria" = "#363781")

plot_sigtab_RA <- ggplot(sigtab_RA, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=6) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 14, color = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  scale_color_manual(values = color.phylum_RA) +
  scale_x_continuous(limits = c(-8, 4), breaks = seq(-8, 4, by = 4), label = c("-8", "-4", "0", "4")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RA

ggsave("plot_sigtab_RA.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/DESeq2_jki_seq1/", width = 24, height = 16, units = "cm", dpi = 200, device = "png")


##Convert variables characters in factors   ---> Done with DESeq2_jki_seq1_like24
# Phylum order
x = tapply(sigtab_RH$log2FoldChange, sigtab_RH$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH$Phylum = factor(as.character(sigtab_RH$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RH$log2FoldChange, sigtab_RH$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH$Annotation = factor(as.character(sigtab_RH$Annotation), levels=names(x))

# Define color palette
color.phylum_RH <- c("Acidobacteriota" = "#ffc331",
                     "Actinobacteriota" = "#ff8b60", 
                     "Bacteroidota"  = "#Bff4be", 
                     "Chloroflexi" = "#2b8f22", 
                     "Firmicutes" = "#b300b3", 
                     "Proteobacteria" = "#363781")

plot_sigtab_RH <- ggplot(sigtab_RH, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=6) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 14, color = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  scale_color_manual(values = color.phylum_RH) +
  scale_x_continuous(limits = c(-6, 3), breaks = seq(-6, 3, by = 3), label = c("-6", "-3", "0", "3")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RH

ggsave("plot_sigtab_RH.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/DESeq2_jki_seq1/", width = 24.5, height = 16, units = "cm", dpi = 200, device = "png")


##Convert variables characters in factors   ---> Done with DESeq2_jki_seq1_like24
# Phylum order
x = tapply(sigtab_RP$log2FoldChange, sigtab_RP$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP$Phylum = factor(as.character(sigtab_RP$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RP$log2FoldChange, sigtab_RP$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP$Annotation = factor(as.character(sigtab_RP$Annotation), levels=names(x))

# Define color palette
color.phylum_RP <- c("Acidobacteriota" = "#ffc331",
                        "Actinobacteriota" = "#ff8b60", 
                        "Bacteroidota"  = "#Bff4be", 
                        "Chloroflexi" = "#2b8f22", 
                        "Firmicutes" = "#b300b3", 
                        "Gemmatimonadota" = "#ca250f",
                        "Latescibacterota" = "#4a2500", 
                        "Methylomirabilota" = "#964B00", 
                        "NB1-j" = "#b5814c", 
                        "Nitrospirota" = "#1361cf",
                        "Patescibacteria" = "#4b0096", 
                        "Proteobacteria" = "#363781", 
                        "RCP2-54" = "#800080",
                        "Verrucomicrobiota" = "#ced2ce") 
                        

plot_sigtab_RP <- ggplot(sigtab_RP, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=6) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 14, color = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  scale_color_manual(values = color.phylum_RP) +
  scale_x_continuous(limits = c(-6, 3), breaks = seq(-6, 3, by = 3), label = c("-6", "-3", "0", "3")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RP

ggsave("plot_sigtab_RP.png", path = "~/Documents/R_analysis/jki_seq1/output_jki_seq1/DESeq2_jki_seq1/", width = 23, height = 37, units = "cm", dpi = 200, device = "png")



## The end! : )

