## Core - jki_seq1
#### Adriana Giongo
#### (21.05.2023) 

#Load packages
library("phyloseq")
library("microbiome")
library("dplyr")

# Calculate compositional version of the data (relative abundances)
#collapse the taxa into "Annotation" level
rank_names(jki_seq1)
aggreg_jki_seq1_annotation <- aggregate_taxa(jki_seq1, "Annotation")

##Compute prevalence of each feature
prev_jki_seq1 = apply(X = otu_table(aggreg_jki_seq1_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_jki_seq1_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to a data frame
prev_jki_seq1 = data.frame(Prevalence=prev_jki_seq1, TotalAbundance = taxa_sums(aggreg_jki_seq1_annotation), tax_table(aggreg_jki_seq1_annotation), otu_table(aggreg_jki_seq1_annotation))

#Write data frame in csv format
write.csv(prev_jki_seq1, "~/path/prev_aggreg_jki_seq1_annotation.csv")

#taxa that exceed the given prevalence and detection thresholds
core_taxa_jki_seq1_annotation <- core_members(aggreg_jki_seq1_annotation, detection = 0, prevalence = 99/100)

#write core taxa as csv table
write.csv(core_taxa_jki_seq1_annotation, "~/path/core_taxa_jki_seq1_annotation.csv")

## The end! : )


