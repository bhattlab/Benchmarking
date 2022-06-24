library(tidyverse)
library(ggplot2)
library(here)

# goal - make a table that is phyla or genus to ncbi id to 16S copy number

rrndb <- read.csv(here("qPCR/16Scopynum/rrnDB-5.7_pantaxa_stats_NCBI.tsv"), header=TRUE, sep="\t")
genbank_taxonomy <- read.csv(here("qPCR/16Scopynum/genbank_taxonomy_array.tsv"), header=FALSE, sep="\t")
names(genbank_taxonomy) <- c("name", "taxid", "root", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")

conversion <- merge(genbank_taxonomy, rrndb, by="taxid", all.x=TRUE) %>% select(taxid, name.x, mean)
names(conversion) <- c("taxid", "name", "mean16S")

write.table(conversion, "taxonomic_copy_number.tsv", sep="\t", row.names=FALSE, quote=FALSE)

