library(here)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)

mytheme<-theme_bw() + 
  theme(text = element_text(size=12), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())



# read in input table
counts <- read.table(here("DNA/4.SNV/1.seqtk_all/counts.tsv"), header=FALSE, sep="\t")
names(counts) <- c("Sample", "A_Count", "C_Count", "G_Count", "T_Count")

counts <- mutate(counts, SampleCode = gsub("_.*", "", Sample))

# group forward and reverse reads into a single row
counts <- counts %>% group_by(SampleCode) %>% summarise(A_Count = sum(A_Count), 
                                                        T_Count = sum(T_Count), 
                                                        C_Count = sum(C_Count), 
                                                        G_Count = sum(G_Count))

# calculate proportion from counts
counts <- mutate(counts, total=A_Count + C_Count + G_Count + T_Count)
counts <- mutate(counts, A_Proportion = A_Count / total)
counts <- mutate(counts, C_Proportion = C_Count / total)
counts <- mutate(counts, G_Proportion = G_Count / total)
counts <- mutate(counts, T_Proportion = T_Count / total)

# make a preservative column
counts <- mutate(counts, Preservative = gsub("D01-", "", gsub(".-R.", "", SampleCode)))
counts <- mutate(counts, Preservative = gsub("O", "OMNIgene", gsub("N", "None", gsub("Z", "Zymo", Preservative))))

# calculate base pair skew
counts <- mutate(counts, A_to_T = A_Count / T_Count)
counts <- mutate(counts, C_to_G = C_Count / G_Count)
base_ratios <- counts %>% select("SampleCode", "Preservative", "A_to_T", "C_to_G")
base_ratios <- melt(base_ratios, id.vars=c("SampleCode", "Preservative"))
names(base_ratios) <- c("SampleCode", "Preservative", "BasePair", "Skew")
ggplot(base_ratios, aes(x=Preservative, y=Skew)) + 
  geom_boxplot(aes(fill=BasePair), alpha=0.5, width=0.5) + 
  geom_point(aes(fill=BasePair), shape=21, color="black", position = position_dodge(width=0.5)) + 
  mytheme

# wide to long
counts <- melt(counts %>% select("SampleCode", "Preservative", "A_Proportion", "T_Proportion", "C_Proportion", "G_Proportion"),
               id.vars=c("SampleCode", "Preservative"))
names(counts) <- c("SampleCode", "Preservative", "Base", "Proportion")

ggplot(counts, aes(x=Preservative, y=Proportion)) + 
  geom_boxplot(aes(fill=Base), alpha=0.5, width=0.5) + 
  geom_point(aes(fill=Base), shape = 21, color="black", position=position_dodge(width=0.5)) +
  mytheme

ggsave(here("DNA/4.SNV/1.seqtk_all/complete_microbiome_composition.pdf"), w=4.5, h=3, dpi=300)
ggsave(here("DNA/4.SNV/1.seqtk_all/complete_microbiome_composition.jpg"), w=4.5, h=3, dpi=300)

