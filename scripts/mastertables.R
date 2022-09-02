library(ggplot2)
library(ggpubr)
library(here)
library(dplyr)
library(forcats)
library(gtable)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(ggnewscale)

##### SET UP DATA TABLES AND WRITE MASTER DATA TABLES #####
# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#read in the QSU relative data as three separate dataframes
raw <- read.csv(here("QSU_Data/raw_data_relative.csv"), header=TRUE) # raw per sample information
model <- read.csv(here("QSU_Data/raw_model_relative.csv"), header=TRUE) # means and confidence intervals for all conditions
sig <- read.csv(here("QSU_Data/raw_significance_relative.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests - # updated this with the 9/1 genus dna abundance

#Filter QSU relative data to remove outdated absolute data
model <- filter(model, feature!="MicrobesPerGram" & feature!="mean_16s_count" & feature!="DNAConcentration")
raw <- raw %>% select(-mean_16s_count)
raw <-raw %>% select(-DNAConcentration)
raw <-raw %>% select(-Patient, -Sample_Type, -Replication)
sig <- filter(sig, feature!="mean_16s_count")
colnames(model)[colnames(model) == "estimate"] <- "prediction"

# read in the QSU absolute data as three separate dataframes 
rawabs <- read.csv(here("QSU_Data/raw_data_all.csv"), header=TRUE) # raw per sample information
modelabs <- read.csv(here("QSU_Data/raw_data_model_means.csv"), header=TRUE) # means and confidence intervals for all conditions
sigabs <- read.csv(here("QSU_Data/raw_data_significance.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

#merge the absolute and relative QSU data
raw <- merge(raw,rawabs, by="Sample")
model <- rbind(model, modelabs)
sig <- rbind(sig,sigabs)

# read in the metagenomic preprocessing readcounts information
readcounts_DNA <- read.table(here("DNA/1.preprocess/preprocessing/01_processing/readcounts.tsv"), header=TRUE, sep="\t")
readcounts_DNA <- readcounts_DNA %>% select(Sample, raw_reads, dedup_reads, trimmed_reads,host_removed_reads,orphan_reads)
names(readcounts_DNA) <- c("Sample", "RawReads_DNA", "DeduplicatedReads_DNA", "TrimmedReads_DNA", "HostRemovedReads_DNA", "OrphanReads_DNA")

# read in the metatranscriptomic preprocessing readcounts information
readcounts_RNA <- read.table(here("RNA/01_rrna/readcounts.tsv"), header=TRUE, sep="\t")
readcounts_RNA <- readcounts_RNA %>% select(SampleID, raw_reads, dedup_reads, trimmed_reads,host_removed_reads,orphan_reads, After_RRNA)
names(readcounts_RNA) <- c("Sample", "RawReads_RNA", "DeduplicatedReads_RNA", "TrimmedReads_RNA", "HostRemovedReads_RNA", "OrphanReads_RNA", "rRNARemovedReads_RNA")
readcounts_RNA <- readcounts_RNA %>% mutate(Sample=gsub(" ", "-", Sample))

# edit existing dataframes such that column names indicate that all values are genomic
colnames(raw) <- paste(colnames(raw),"_DNA", sep="")
colnames(raw)[1] <- "Sample"
model <- model %>% mutate(measurement = "DNA")
sig <- sig %>% mutate(measurement = "DNA")

# merge preprocessing dataframes with raw dataframe
raw <- merge(raw, readcounts_DNA, by="Sample")
raw <- merge(raw, readcounts_RNA, by="Sample")

# NEXT, add in transcriptomic data

# read in the QSU absolute data as three separate dataframes 
rawrna <- read.csv(here("QSU_Data/raw_data_rna_all.csv"), header=TRUE) # raw per sample information
modelrna <- read.csv(here("QSU_Data/raw_data_rna_model_means.csv"), header=TRUE) # means and confidence intervals for all conditions
sigrna <- read.csv(here("QSU_Data/raw_data_rna_significance.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

# edit RNA dataframes such that column names indicate that all values are transcriptomic
colnames(rawrna) <- paste(colnames(rawrna),"_RNA", sep="")
colnames(rawrna)[1] <- "Sample"
colnames(rawrna)[2] <- "Patient"
colnames(rawrna)[3] <- "SampleType"
colnames(rawrna)[4] <- "Replication"
modelrna <- modelrna %>% mutate(measurement = "RNA")
sigrna <- sigrna %>% mutate(measurement = "RNA")

# merge all dataframes
raw <- merge(raw, rawrna, by="Sample")
model <- rbind(model, modelrna)
sig <- rbind(sig, sigrna)

# NEXT, read in raw sample data to get DNA and RNA concentrations
dna_meta <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
dna_meta <- dna_meta %>% select(SampleID, DNAConcentration)
colnames(dna_meta) <- c("Sample", "DNAConcentration")
dna_meta <- dna_meta %>% mutate(Sample=gsub("_", "-", Sample))
rna_meta <- read.csv(here("data/RNAExtraction.tsv"), sep="\t", header=TRUE)
rna_meta <- rna_meta %>% select(SampleID, RNAConcentration)
colnames(rna_meta) <- c("Sample", "RNAConcentration")
rna_meta <- rna_meta %>% mutate(Sample=gsub("_", "-", Sample))


raw <- merge(raw, dna_meta, by="Sample")
raw <- merge(raw, rna_meta, by="Sample")

#save to Github repo
write.table(raw, here("QSU_Data/raw_data_full.csv"), row.names = FALSE, sep=",", quote=FALSE)
write.table(model, here("QSU_Data/model_data_full.csv"), row.names = FALSE, sep=",", quote=FALSE)
write.table(sig, here("QSU_Data/sig_data_full.csv"), row.names = FALSE, sep=",", quote=FALSE)
