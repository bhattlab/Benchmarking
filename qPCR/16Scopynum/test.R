library(here)
library(tidyverse)
library(reshape2)




##### MERGING TAXONOMIC ARRAY WITH KRAKEN OUTPUT TO FIND MATCHES #####
kraken_genus <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
tax_array <- read.csv(here("qPCR/16Scopynum/genbank_taxonomy_array.tsv"), header=FALSE, sep="\t")
names(tax_array) <- c("Name", "taxid", "root", "kingdom", "phylum", "order", "class", "family", "Genus", "species", "strain")
tax_array <- tax_array %>% select("Name", "taxid")
names(tax_array) <- c("Genus", "taxid")

kraken_genus$Genus <- row.names(kraken_genus)
row.names(kraken_genus) <- NULL
kraken_genus <- kraken_genus %>% select(Genus) %>% unique()

comparison <- merge(kraken_genus, tax_array, by="Genus", all.x=TRUE)
comparison %>% filter(is.na(taxid))

comparison <- comparison %>% mutate(taxid=ifelse(is.na(taxid), gsub(".*\\(", "", gsub("\\)", "", Genus)), taxid))
comparison <- comparison %>% mutate(Genus=gsub("\\(.*\\)", "", Genus))
comparison <- comparison %>% group_by(Genus) %>% filter(row_number() == 1)

copyndb <- read.csv(here("qPCR/taxonomic_copy_number.tsv"), sep="\t", header=TRUE)
colnames(copyndb)<- c("taxid", "Taxon", "mean16S")

test <- merge(comparison, copyndb, by="taxid", all.x=TRUE) %>% select("taxid", "Genus", "mean16S")

# re-read in taxonomy array to get domain info
tax_array <- read.csv(here("qPCR/16Scopynum/genbank_taxonomy_array.tsv"), header=FALSE, sep="\t")
names(tax_array) <- c("Name", "taxid", "root", "kingdom", "phylum", "order", "class", "family", "Genus", "species", "strain")
tax_array <- tax_array %>% select("taxid", "kingdom")
test <- merge(test, tax_array, by="taxid")

test <- mutate(test, mean16S=ifelse(is.na(mean16S), ifelse(kingdom == "Bacteria", 1.94, 0), mean16S))

write.table(test, here("qPCR/genus_copynumber_table.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

# FIGURING OUT HOW MUCH THE UNREPRESENTED MICROBES ACCOUNT FOR
kraken_genus <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
kraken_genus$Genus <- row.names(kraken_genus)
row.names(kraken_genus) <- NULL
kraken_genus <- kraken_genus %>% filter(Genus %in% uniq$Genus)
a <- colSums(Filter(is.numeric, kraken_genus))
head(kraken_genus)
b <- data.frame(a)

# FOR ANYTHING NOT IN RRNDB, ASSIGN 0 OR CONSTANT BASED ON NON-BACTERIAL OR BACTERIAL



###### phylum #####
kraken_phylum <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_phylum_percentage.txt"), sep="\t", header=TRUE)
tax_array <- read.csv(here("qPCR/16Scopynum/genbank_taxonomy_array.tsv"), header=FALSE, sep="\t")
names(tax_array) <- c("Name", "taxid", "root", "kingdom", "phylum", "order", "class", "family", "Genus", "species", "strain")
tax_array <- tax_array %>% select("Name", "taxid")
names(tax_array) <- c("Phylum", "taxid")

kraken_phylum$Phylum <- row.names(kraken_phylum)
row.names(kraken_phylum) <- NULL
kraken_phylum <- kraken_phylum %>% select(Phylum) %>% unique()


comparison <- merge(kraken_phylum, tax_array, by="Phylum", all.x=TRUE)
comparison %>% filter(is.na(taxid))

comparison <- comparison %>% mutate(taxid=ifelse(is.na(taxid), gsub(".*\\(", "", gsub("\\)", "", Phylum)), taxid))
comparison <- comparison %>% mutate(Phylum=gsub("\\(.*\\)", "", Phylum))
comparison <- comparison %>% group_by(Phylum) %>% filter(row_number() == 1)


# now try grouping with rrnDB
copyndb <- read.csv(here("qPCR/taxonomic_copy_number.tsv"), sep="\t", header=TRUE)
colnames(copyndb)<- c("taxid", "Taxon", "mean16S")

test <- merge(comparison, copyndb, by="taxid", all.x=TRUE) %>% select("taxid", "Phylum", "mean16S")



# re-read in taxonomy array to get domain info
tax_array <- read.csv(here("qPCR/16Scopynum/genbank_taxonomy_array.tsv"), header=FALSE, sep="\t")
names(tax_array) <- c("Name", "taxid", "root", "kingdom", "phylum", "order", "class", "family", "Genus", "species", "strain")
tax_array <- tax_array %>% select("taxid", "kingdom")
test <- merge(test, tax_array, by="taxid")

test <- mutate(test, mean16S=ifelse(is.na(mean16S), ifelse(kingdom == "Bacteria", 1.94, 0), mean16S))

write.table(test, here("qPCR/phylum_copynumber_table.tsv"), row.names=FALSE, quote=FALSE, sep="\t")


uniq <- test %>% filter(is.na(mean16S)) %>% select(Phylum) %>% unique()



