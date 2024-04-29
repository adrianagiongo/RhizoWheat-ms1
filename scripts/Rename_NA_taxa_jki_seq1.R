#Source = https://github.com/joey711/phyloseq/issues/850

#Load packages
library("phyloseq")
library("stringr")

#Export tax table
tax <- data.frame(tax_table(jki_seq1_data))

#Change NA to a empty string (changing the script to use is.na() is also an option)
tax.clean <- data.frame(row.names = row.names(tax), Kingdom = str_replace(tax[,1],"NA",""),
                        Phylum = str_replace(tax[,2], "NA",""),
                        Class = str_replace(tax[,3], "NA",""),
                        Order = str_replace(tax[,4], "NA",""),
                        Family = str_replace(tax[,5],"NA", ""),
                        Genus = str_replace(tax[,6], "NA",""),
                        Annotation = str_replace(tax[,7],"NA", ""),
                        stringsAsFactors = FALSE)
                        
#Change all columns to characters (otherwise everything becomes NA)
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

#Fill missing taxonomy
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste(tax.clean[i,1])
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste(tax.clean[i,2])
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste(tax.clean[i,3])
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste(tax.clean[i,4])
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste(tax.clean[i,5])
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Annotation[i] <- paste(tax.clean$Genus[i])
  }
}

#Return data.frame to a phyloseq object
tax_table(jki_seq1_data) <- as.matrix(tax.clean)
head(tax_table(jki_seq1_data))
tail(tax_table(jki_seq1_data))

## The end : )
#go through the next scripts ("Clean_dataset.R")

