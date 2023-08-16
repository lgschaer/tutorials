



### Create a phyloseq object with the new inocula sequencing data from Zymo to look at the decrease in diversity during the enrichment process.

library(phyloseq)
library(tidyverse)
library(csv)



toy_ps <- read_rds("/home/lgschaer/old/Tutorials/toy_phyloseq.rds")  
toy_ps


#This is how I made the inputs for the tutorial using toy_ps
#sdata <- as_tibble(sample_data(toy_ps))
#head(sdata)
#write_csv(sdata, "/home/lgschaer/old/Tutorials/phyloseq/sdata.csv")
#seqtab <- otu_table(toy_ps)
#colnames(seqtab) <- colnames(sequence_table)
#seqtab[1:5,1:5]
#dim(seqtab)
#write_rds(seqtab, "/home/lgschaer/old/Tutorials/phyloseq/seqtab.rds")
#taxtab <- tax_table(toy_ps)
#taxtab[1:5,1:5]
#write_rds(taxtab, "/home/lgschaer/old/Tutorials/phyloseq/taxtab.rds")

###create a phyloseq object with the new sequencing data from Zymo
#load metadata

sdata <- as.csv("/home/lgschaer/old/Tutorials/phyloseq/sdata.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/Tutorials/phyloseq/seqtab.rds")
sequence_table[1:3,1:3]

colnames(sequence_table) <- NULL
#sequence_table <- as.data.frame(sequence_table)
sequence_table[1:5,1:5]


sequence_table <- as.matrix(sequence_table)                                #change to matrix format
m <- (colSums(sequence_table, na.rm=TRUE) != 0)                        #T if colSum is not 0, F otherwise
nonzero <- sequence_table[, m]   



#check dimensions
dim(sdata)
dim(nonzero)


#load taxa table
taxa <- readRDS("/home/lgschaer/old/Tutorials/phyloseq/taxtab.rds")
taxa[1:5,1:5]

dim(nonzero)
dim(taxa)

#make phyloseq object
samdata = sample_data(sdata)
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)
taxtab = tax_table(taxa)

sample_names(samdata)
sample_names(seqtab)


taxa_names(seqtab)
taxa_names(taxtab)

#combine all components into a phyloseq object
toy_ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
toy_ps


write_rds(toy_ps, "/home/lgschaer/old/Tutorials/phyloseq/toy_ps_out.rds")



#Filter out eukaryotes and mitochondria
toy_ps2 <- toy_ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
toy_ps2

### normalize

sample_sums(toy_ps2)

samplesover1000_all <- subset_samples(toy_ps2, sample_sums(toy_ps2) > 1000)

any(taxa_sums(samplesover1000_all) == 0)

sum(taxa_sums(samplesover1000_all) == 0)

prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)

any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
set.seed(81)

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

toy_ps3<- rarefy_samplesover1000_all
toy_ps3


