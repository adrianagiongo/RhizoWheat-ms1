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
  - Relevant metadata parameters (example below)

| Sample     | Sample_name | Sample_id | Core | Microhabitat | Rotation | Replicates | Layer | Rotlayer | Depth_cm | Group  | C     | N     | Ngkg | SoilMoisture |
| ---------- | ----------- | --------- | ---- | ------------ | -------- | ---------- | ----- | -------- | -------- | ------ | ----- | ----- | ---- | ------------ |
| S105W1a1RA | WR1L1RA     | WR1L1     | C1   | RA           | WR       | WR1        | L1    | WRL1     | 0-15     | WRL1RA | 1.344 | 0.161 | 1.61 | 0.2          |
| S106W1a1RH | WR1L1RH     | WR1L1     | C1   | RH           | WR       | WR1        | L1    | WRL1     | 0-15     | WRL1RH | 1.365 | 0.161 | 1.61 | 0.21         |
| S107W1a1RP | WR1L1RP     | WR1L1     | C1   | RP           | WR       | WR1        | L1    | WRL1     | 0-15     | WRL1RP | 1.36  | 0.164 | 1.64 | 0.2          |
|            |
