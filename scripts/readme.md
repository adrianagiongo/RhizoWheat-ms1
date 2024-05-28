### Github repository for 

## Dataset jki_seq1 (16S rRNA gene amplicon sequencing) pipeline analyses
Performed in R v.4.1.3 [R Core Team](https://www.r-project.org)

#### Color code
- Rotations
  - ![#B1AD31](https://placehold.co/15x15/B1AD31/B1AD31.png) `#B1AD31` (WR)
  - ![#196418](https://placehold.co/15x15/196418/196418.png) `#196418` (WM)
- Microhabitats
  - ![#b05644](https://placehold.co/15x15/b05644/b05644.png) `#b05644` (RA)
  - ![#d9b967](https://placehold.co/15x15/d9b967/d9b967.png) `#d9b967` (RH)
  - ![#57896A](https://placehold.co/15x15/57896A/57896A.png) `#57896A` (RP)
- Soil layers
  - ![#3B0404](https://placehold.co/15x15/3B0404/3B0404.png) `#3B0404` (L1)
  - ![#B95C50](https://placehold.co/15x15/B95C50/B95C50.png) `#B95C50` (L2)
  - ![#DE847B](https://placehold.co/15x15/DE847B/DE847B.png) `#DE847B` (L3)
  - ![#DEB3AD](https://placehold.co/15x15/DEB3AD/DEB3AD.png) `#DEB3AD` (L4)
  - ![#DCC8BA](https://placehold.co/15x15/DCC8BA/DCC8BA.png) `#DCC8BA` (L5)

### Packages
library("dada2")\
library("ShortRead")\
library("Biostrings")\
library("phyloseq")\
library("microbiome")\
library("vegan")\
library("DESeq2") \
library("dplyr")\
library("stringr")\
library("ggpubr")\
library("tidyr")\
library("ggplot2")\
library("readxl")\
library("RColorBrewer")
  
### 1. Dada2
This script uses the raw data obtained from the BioProject SRA.\
A R server is required. \
Database used: [SILVA 138 SSU](https://www.arb-silva.de/documentation/release-138/) 
- Non-truncate parameter
```
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress = TRUE, truncQ = 2, trimLeft=c(17,20), maxN=0, maxEE=c(2,2), rm.phix=TRUE, matchIDs = TRUE, multithread=TRUE)
```

### 2. Creating phyloseq object
This script creates a phyloseq object based on these files:

- jki_seq1_metadata.csv
- jki_seq1_otu.xlsx
- jki_seq1_seqs.xlsx
- jki_seq1_taxa.xlsx

### 3. Rename NA
This script replaces NA for the latest taxonomy found for an ASV.

### 4. Clean dataset
This script removes unwanted taxonomic groups from the dataset.
- Root NAs (Domain, Phylum)
- Eukaryotes (Domain)
- Chloroplasts (Order)
- Mitochondria (Family)

### 5. Data selection
This script selects a group of samples to be analyzed separately.

### 6. Rarefaction
This script performs rarefaction based on the minimum sequences. See Schloss (2024).

### 7. Alpha diversity
This script calculates alpha diversity based on the rarefied data.

### 8. Ordination 
This script creates MDS plots and calculates PERMANOVA and ANOSIM based on the rarefied data.

### 9. DESeq2
Based on the rarefied data, this script performs differential abundance (DA) between two groups.

### References
Callahan B, McMurdie P, Rosen M et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat. Met. 13:581–583. [DOI](http://10.1101/024034)

Schloss PD (2024) Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. mSphere 9:e00354-23. [DOI](http://10.1128/msphere.00354-23)

