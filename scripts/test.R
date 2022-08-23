library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(reshape2)

metadata <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)

#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
#Separate the Condition column to create preservation method and temperature columns
metadata <- mutate(metadata, Preservation=substr(Condition,1,1))
metadata <- mutate(metadata, Temperature=substr(Condition,2,2))
metadata <- metadata %>% mutate(TemperatureLong = ifelse(Temperature == "F", "-80°C", ifelse(Temperature == "R", "23°C", "40°C"))) # write out temperature
metadata <- mutate(metadata, TemperatureLong=ifelse(Replicate == "R2", TemperatureLong, ""))
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))

phylum <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_phylum_percentage.txt"), sep="\t", header=TRUE)
phylum$Phylum <- row.names(phylum)
kraken_long <- melt(phylum, id.vars = "Phylum", variable.name = "Sample", value.name = "rel_abundance")
kraken_long <- mutate(kraken_long, SampleID=gsub("\\.", "_", Sample))
kraken_long <- merge(kraken_long, metadata, by="SampleID")

kraken_long <- kraken_long %>% filter(Phylum == "Firmicutes" | Phylum == "Bacteroidetes")
kraken_long <- kraken_long %>% filter(Donor == "D10")
kraken_long <- kraken_long %>% mutate(Grouping = ifelse(Condition == "NF", "Baseline", ifelse(Condition == "OF" | Condition == "ZF", "Kit", ifelse(Condition == "OR" | Condition == "ZR", "Time", "Temperature"))))

ggplot(kraken_long, aes(x=Grouping, y=rel_abundance, color = Phylum)) + 
  geom_point(size=2) + 
  theme_classic() +
  ylab("Relative Abundance") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=12))


summary <- kraken_long %>% 
  group_by(Phylum, Grouping) %>%
  summarise_at(vars(rel_abundance), list(Mean = mean, Max=max, Min=min))
