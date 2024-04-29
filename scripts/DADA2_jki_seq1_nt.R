library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)

path <- "~/jki_seq1/jki_seq1/" # change to the right location
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern = "_R1_001.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

#plotQualityProfile(system.file, n = 5e+05, aggregate = FALSE)
pdf("~/jki_seq1/jki_seq1_dada2_output/plot_Quality_Forward.pdf")
plotQualityProfile(fnFs[1:6])
dev.off()

pdf("~/jki_seq1/jki_seq1_dada2_output/plot_Quality_Reverse.pdf")
plotQualityProfile(fnRs[1:6])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))

# Filtering on Windows,set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress = TRUE, truncQ = 2, trimLeft=c(17,20), maxN=0, maxEE=c(2,2), rm.phix=TRUE, matchIDs = TRUE, multithread=TRUE)
head(out)

#plot quality after filter
pdf("~/jki_seq1/jki_seq1_dada2_output/Quality_Forward2.pdf")
plotQualityProfile(filtFs[1:6])
dev.off()

pdf("~/jki_seq1/jki_seq1_dada2_output/Quality_Reverse2.pdf")
plotQualityProfile(filtRs[1:6])
dev.off()

#set run
set.seed(2022)

#predict error and plot
errF <- learnErrors(filtFs, multithread=TRUE)
pdf("~/jki_seq1/jki_seq1_dada2_output/plot_error_Forward.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

errR <- learnErrors(filtRs, multithread=TRUE)
pdf("~/jki_seq1/jki_seq1_dada2_output/plot_error_Reverse.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#dereplicate sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#run dada
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#save.image(file = "~/jki_seq3/jki_seq3_16S/jki_seq3_16S_ggt/jki_seq3_16S_ggt_dada2_output/JKI_seq3_16S_ggt.RData")
#load("~/jki_seq3/jki_seq3_16S/jki_seq3_16S_ggt/jki_seq3_16S_ggt_dada2_output/JKI_seq3_16S_ggt.RData")

#merge paired sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#get sequences, remove chimera, and write tables
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
seqtab.nochim.t<-t(seqtab.nochim)

#Export results
write.csv(seqtab.nochim.t, "~/jki_seq1/jki_seq1_dada2_output/jki_seq1_Phy_Table.csv")

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

##If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "~/jki_seq1/jki_seq1_dada2_output/jki_seq1_Quality_Control.csv")

saveRDS(seqtab.nochim, file = "~/jki_seq1/seqtab.nochim_jki_seq1.rds")
save.image(file = "~/jki_seq1/jki_seq1_dada2.RData")

#Load workspace or object
#load("~/jki_seq3/jki_seq3_16S/jki_seq3_16S_ggt_dada2.RData")
#seqtab.nochim <- readRDS("~/jki_seq3/jki_seq3_16S/seqtab.nochim_jki_seq3_16S_ggt.rds")

#Assign Taxonomy at Genus level based on the minBoot of bootstrap confidence - minBoot=50 is default
taxa <- assignTaxonomy(seqtab.nochim, "~/database/silva_nr_v138_train_set.fa", minBoot = 0, outputBootstraps = TRUE, multithread=TRUE)
write.csv(taxa, "~/jki_seq1/jki_seq1_dada2_output/jki_seq1_taxa_silva138_table.csv")

##Assign Taxonomy at Genus level and add Species based on the exact maching
#silva_138
#taxa <- assignTaxonomy(seqtab.nochim, "~/database/silva_nr_v138_train_set.fa", multithread=TRUE)
#taxa <- addSpecies(taxa, "~/database/silva_species_assignment_v138.fa")
#write.csv(taxa, "~/jki_seq1/jki_seq1/jki_seq1/jki_seq1_dada2_output/jki_seq1_taxa_species_silva138_table.csv")

