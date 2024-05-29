### Github repository for 

## Dataset jki_seq1 (16S rRNA gene amplicon sequencing) pipeline analyses
Performed in R v.4.1.3 [R Core Team](https://www.r-project.org)

#### Files
- ASV (former OTU) information
  - ASVs names
  - Number of sequences for each OTU in each sample (original name of samples)
- Taxonomic information (jki_seq1_taxa) using Silva Database
  - ASVs names
  - Six taxonomic levels (Kingdom, Phylum, Class, Order, Family, Genus) 
  - One taxonomic level (Annotation/Taxa) created based on the Genus level (containing NAs)
- Sequences
  - ASVs names
  - Sequences of each ASV
- Metadata
  - Original name of samples (Sample)
  - Relevant metadata parameters (Sample_name,	Sample_id,	Core,	Microhabitat,	Rotation,	Replicates,	Layer,	Rotlayer,	Depth_cm,	Group,	C,	N,	Ngkg,	SoilMoisture
