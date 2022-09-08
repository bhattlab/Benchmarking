library(ggpubr)
library(here)
library(dplyr)
library(forcats)
library(gtable)
library(tidyr)
library(reshape2)

# Read in main tables

# DNA tables
dna_raw <- read.csv(here("QSU_Data/raw_data_dna_full.csv"), sep=",", header=TRUE)
dna_model <- read.csv(here("QSU_Data/model_data_dna_full.csv"), sep=",", header=TRUE)
dna_sig <- read.csv(here("QSU_Data/sig_data_dna_full.tsv"), sep="\t", header=TRUE)


# RNA tables
rna_raw <- read.csv(here("QSU_Data/raw_data_rna_full.csv"), sep=",", header=TRUE)
rna_model <- read.csv(here("QSU_Data/model_data_rna_full.csv"), sep=",", header=TRUE)
rna_sig <- read.csv(here("QSU_Data/sig_data_rna_full.tsv"), sep="\t", header=TRUE)


# DNA Bray Curtis tables
dna_bc_raw <- read.csv(here("QSU_Data/raw_data_all_braycurtis.csv"), header=TRUE)
dna_bc_model <- read.csv(here("QSU_Data/raw_data_model_means_braycurtis.csv"), header=TRUE)
dna_bc_sig <- read.csv(here("QSU_Data/raw_data_significance_braycurtis.csv"), header=TRUE)


# RNA Bray Curtis tables
rna_bc_raw <- read.csv(here("QSU_Data/raw_data_rna_all_braycurtis.csv"), header=TRUE)
rna_bc_model <- read.csv(here("QSU_Data/raw_data_rna_model_means_braycurtis.csv"), header=TRUE)
rna_bc_sig <- read.csv(here("QSU_Data/raw_data_rna_significance_braycurtis.csv"), header=TRUE)


# build DNA raw table
dna_raw <- dna_raw %>% select("Sample", "Sample_Type","Replication",
                              "DNAConcentration", "mean_16s_count", "MicrobesPerGram", "ratio.of.B.F",
                              "RawReads_DNA", "DeduplicatedReads_DNA", "TrimmedReads_DNA", "HostRemovedReads_DNA", "OrphanReads_DNA",
                              "Absolute.Abundance..Firmicutes", "Absolute.Abundance..Bacteroidetes", "Absolute.Abundance..Actinobacteria", "Absolute.Abundance..Other.Bacteria","Absolute.Abundance..Viruses","Absolute.Abundance..Fungi", "Absolute.Abundance..Unclassified...Root", "Absolute.Abundance..Unclassified",
                              "Shannon.Entropy", "Inv.Simpson","Richness.0.01.",
                              "Relative.Abundance..Firmicutes","Relative.Abundance..Bacteroidetes", "Relative.Abundance..Actinobacteria", "Relative.Abundance..Other.Bacteria", "Relative.Abundance..Viruses", "Relative.Abundance..Fungi", "Relative.Abundance..Unclassified...Root", "Relative.Abundance..Unclassified")
                                
                                            
                                            
names(dna_raw) <- c("Sample", "Condition","Replication",
                    "DNAConcentration", "Mean16sCount", "MicrobesPerGram", "RatioBacteroidetesFirmicutes",
                    "RawReads", "DeduplicatedReads", "TrimmedReads", "HostRemovedReads", "OrphanReads",
                    "AbsoluteAbundanceFirmicutes", "AbsoluteAbundanceBacteroidetes", "AbsoluteAbundanceActinobacteria", "AbsoluteAbundanceOtherBacteria","AbsoluteAbundanceViruses","AbsoluteAbundanceFungi", "AbsoluteAbundanceUnclassifiedRoot", "AbsoluteAbundanceUnclassified",
                    "ShannonEntropy", "InverseSimpson","RichnessAbove0.01",
                    "RelativeAbundanceFirmicutes","RelativeAbundanceBacteroidetes", "RelativeAbundanceActinobacteria", "RelativeAbundanceOtherBacteria", "RelativeAbundanceViruses", "RelativeAbundanceFungi", "RelativeAbundanceUnclassifiedRoot", "RelativeAbundanceUnclassified")                                 

# build DNA model table
dna_model <- dna_model %>% filter(!feature %in% c("Richness 0.001%", "Richness 0.0001%", "Richness 0.00001%", "Richness Any%", "Richness 100M", "Richness 10M", "Richness 1M"))
names(dna_model) <- c("feature", "Condition", "prediction", "CI_low", "CI_high")

# build DNA sig table
dna_sig <- dna_sig %>% filter(!feature %in% c("Richness 0.001%", "Richness 0.0001%", "Richness 0.00001%", "Richness Any%", "Richness 100M", "Richness 10M", "Richness 1M"))
names(dna_sig) <- c("feature", "Condition1", "Condition2", "PercentChange", "p", "p.adj","p.format","p.signif","method")
dna_sig <- dna_sig %>% select("feature", "Condition1", "Condition2", "PercentChange", "p", "p.format")


# build RNA raw table
rna_raw <- rna_raw %>% select("Sample", "Sample_Type","Replication",
                              "RNAConcentration", 
                              "RawReads", "DeduplicatedReads", "TrimmedReads", "HostRemovedReads", "OrphanReads", "rRNARemovedReads",
                              "Shannon.Entropy", "Inv.Simpson","Richness.0.01.",
                              "Relative.Abundance..Firmicutes","Relative.Abundance..Bacteroidetes", "Relative.Abundance..Actinobacteria", "Relative.Abundance..Other.Bacteria", "Relative.Abundance..Viruses", "Relative.Abundance..Fungi", "Relative.Abundance..Unclassified...Root", "Relative.Abundance..Unclassified")



names(rna_raw) <- c("Sample", "Condition","Replication",
                    "RNAConcentration", 
                    "RawReads", "DeduplicatedReads", "TrimmedReads", "HostRemovedReads", "OrphanReads", "rRNARemovedReads",
                    "ShannonEntropy", "InverseSimpson","RichnessAbove0.01",
                    "RelativeAbundanceFirmicutes","RelativeAbundanceBacteroidetes", "RelativeAbundanceActinobacteria", "RelativeAbundanceOtherBacteria", "RelativeAbundanceViruses", "RelativeAbundanceFungi", "RelativeAbundanceUnclassifiedRoot", "RelativeAbundanceUnclassified")                                 

# build RNA model table
rna_model <- rna_model %>% filter(!feature %in% c("Richness 0.001%", "Richness 0.0001%", "Richness 0.00001%", "Richness Any%", "Ratio of Bacteroidetes to Firmicutes"))
names(rna_model) <- c("feature", "Condition", "prediction", "CI_low", "CI_high")

# build RNA sig table
rna_sig <- rna_sig %>% filter(!feature %in% c("Richness 0.001%", "Richness 0.0001%", "Richness 0.00001%", "Richness Any%", "Ratio of Bacteroidetes to Firmicutes"))
names(rna_sig) <- c("feature", "Condition1", "Condition2", "PercentChange", "p", "p.adj","p.format","p.signif","method")
rna_sig <- rna_sig %>% select("feature", "Condition1", "Condition2", "PercentChange", "p", "p.format")

# build Bray Curtis tables
dna_bc_sig <- dna_bc_sig %>% select("Feature", "group1", "group2", "p", "p.format")
rna_bc_sig <- rna_bc_sig %>% select("feature", "group1", "group2", "p", "p.format")
names(rna_bc_sig)[1] <- "Feature"


# NEXT, write out all the tables
write.table(dna_raw, here("outputs/tables/Supplementary/Metagenomic_fulldata.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(dna_model, here("outputs/tables/Supplementary/Metagenomic_model.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(dna_sig, here("outputs/tables/Supplementary/Metagenomic_significance.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

write.table(rna_raw, here("outputs/tables/Supplementary/Metatranscriptomic_fulldata.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(rna_model, here("outputs/tables/Supplementary/Metatranscriptomic_model.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(rna_sig, here("outputs/tables/Supplementary/Metatranscriptomic_significance.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

write.table(dna_bc_raw, here("outputs/tables/Supplementary/Metagenomic_bray_curtis_fulldata.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(dna_bc_model, here("outputs/tables/Supplementary/Metagenomic_bray_curtis_model.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(dna_bc_sig, here("outputs/tables/Supplementary/Metagenomic_bray_curtis_significance.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

write.table(rna_bc_raw, here("outputs/tables/Supplementary/Metatranscriptomic_bray_curtis_fulldata.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(rna_bc_model, here("outputs/tables/Supplementary/Metatranscriptomic_bray_curtis_model.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
write.table(rna_bc_sig, here("outputs/tables/Supplementary/Metatranscriptomic_bray_curtis_significance.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
