#Create rarefied data from phyloseq object using vegan package or microbiome package
#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("ggplot2")
library("tidyr")

##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
tab <- otu_table(psO_jki_seq1)
class(tab) <- "matrix" # as.matrix() will do nothing ## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)
#rarecurve(t(otu_table(psO_jki_seq1)), step=20, cex=0.5, col = "blue")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq1_rarefied<-rarefy_even_depth(psO_jki_seq1, rngseed=2022, sample.size = min(sample_sums(psO_jki_seq1)),trimOTUs=TRUE)
psO_jki_seq1_rarefied
sample_sums(psO_jki_seq1_rarefied)

#Prepare file with meta data using meta function
psO_jki_seq1_rarefied.meta <- meta(psO_jki_seq1_rarefied)
head(psO_jki_seq1_rarefied.meta)
psO_jki_seq1_rarefied.meta.rotation <- subset(psO_jki_seq1_rarefied.meta, select = -c(Sample_name, Sample_id, Core, Microhabitat, Replicates, Layer, Rotlayer, Depth_cm, Group, C, N, Ngkg, SoilMoisture))
head(psO_jki_seq1_rarefied.meta.rotation)

#Create and transpose matrix to make x axis for samples
df_psO_jki_seq1_rarefied <- data.frame(otu_table(psO_jki_seq1_rarefied))
#write.csv(df_psO_jki_seq1_rarefied, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1_rarefied.csv")

df_psO_jki_seq1_rarefied_t2 <- t(df_psO_jki_seq1_rarefied)

#Combine matrix in one
ASVs_rarefied_ed<-cbind(psO_jki_seq1_rarefied.meta.rotation, df_psO_jki_seq1_rarefied_t2)
ASVs_rarefied_ed[1:120,1]

##check for the distribution of any ASV counts in the sample group using "hist" function to plot, 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05) 
#then the null hypotesis (the population is normally distributed) is rejected
hist(ASVs_rarefied_ed$sp1)
shapiro.test(ASVs_rarefied_ed$sp1)
qqnorm(ASVs_rarefied_ed$sp1)

#Define colours
colour = rep(NA, length=length(ASVs_rarefied_ed[,1]))
colour[which(ASVs_rarefied_ed$Rotation=="WR")] = "#B1AD31"
colour[which(ASVs_rarefied_ed$Rotation=="WM")] = "#196418"


dim(ASVs_rarefied_ed)

#plot
tiff("~/Documents/R_analysis/jki_seq1/output_jki_seq1/Alpha_div_jki_seq1/rarecurve_jki_seq1_rarefied.tiff", units="cm", width=12, height=12, res=300)
rarecurve(ASVs_rarefied_ed [,2:27093], step=100, label=FALSE, col=colour, main="", xlab = "Reads", ylab = "ASVs")     #Package VEGAN
legend(legend=c("WR", "WM"),
       "bottomright", 
       bty="n",
       col=c("#B1AD31", "#196418"),
       pch=15,
       pt.cex=1.5,
       cex=0.75,
       ncol=3)
arrows(x0=41961, y0=300, y1=1800, angle=90, length=0, lty = 2)
dev.off()



###################  Unrarefied 

# #Prepare file with meta data using meta function
# psO_jki_seq1.meta <- meta(psO_jki_seq1)
# head(psO_jki_seq1.meta)
# psO_jki_seq1.meta.rotation <- subset(psO_jki_seq1.meta, select = -c(Sample_name, Sample_id, Core, Microhabitat, Replicates, Layer, Rotlayer, Depth_cm, Group, C, N, Ngkg, SoilMoisture))
# head(psO_jki_seq1.meta.rotation)
# 
# #Create and transpose matrix to make x axis for samples
# df_psO_jki_seq1 <- data.frame(otu_table(psO_jki_seq1))
# #write.csv(df_psO_jki_seq1, "~/Documents/R_analysis/jki_seq1/output_jki_seq1/Tables_jki_seq1/df_psO_jki_seq1.csv")
# 
# df_psO_jki_seq1_t2 <- t(df_psO_jki_seq1)
# 
# #Combine matrix in one
# ASVs_ed<-cbind(psO_jki_seq1.meta.rotation, df_psO_jki_seq1_t2)
# ASVs_ed[1:120,1]
# 
# ##check for the distribution of any ASV counts in the sample group using "hist" function to plot, 
# #"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# # If alpha result is lower than the alpha value chosen (p<0.05) 
# #then the null hypotesis (the population is normally distributed) is rejected
# hist(ASVs_ed$sp1)
# shapiro.test(ASVs_ed$sp1)
# qqnorm(ASVs_ed$sp1)
# 
# #Define colours
# colour = rep(NA, length=length(ASVs_ed[,1]))
# colour[which(ASVs_ed$Replicate=="WR")] = "#B1AD31"
# colour[which(ASVs_ed$Replicate=="WM")] = "#196418"
# 
# dim(ASVs_ed)
# 
# #plot
# tiff("~/Documents/R_analysis/jki_seq1/output_jki_seq1/Alpha_div_jki_seq1/rarecurve_jki_seq1_unrarified.tiff", units="cm", width=12, height=12, res=300)
# rarecurve(ASVs_ed [,2:27093], step=100, label=FALSE, col=colour, main="", xlab = "Reads", ylab = "ASVs")     #Package VEGAN
# legend(legend=c("WR", "WM"),
#        "bottomright", 
#        bty="n",
#        col=c("#B1AD31", "#196418"),
#        pch=15,
#        pt.cex=1.5,
#        cex=0.75,
#        ncol=3)
# arrows(x0=41961, y0=350, y1=1710, angle=90, length=0, lty = 2)
# dev.off()
# 
