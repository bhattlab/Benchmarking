library(here)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggtext)

mytheme<-theme_bw() + 
  theme(text = element_text(size=12), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        plot.title = element_text(size = 12))


##### FULL SNV ANALYSIS #####

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

p1 <- ggplot(counts, aes(x=Preservative, y=Proportion)) + 
  geom_boxplot(aes(fill=Base), alpha=0.5, width=0.5) + 
  geom_point(aes(fill=Base), shape = 21, color="black", position=position_dodge(width=0.5)) +
  ggtitle("All Reads") +
  mytheme
p1
ggsave(here("DNA/4.SNV/1.seqtk_all/complete_microbiome_composition.pdf"), w=4.5, h=3, dpi=300)
ggsave(here("DNA/4.SNV/1.seqtk_all/complete_microbiome_composition.jpg"), w=4.5, h=3, dpi=300)


##### F PRAUSNITZII SNV ANALYSIS ####
# read in input table
counts <- read.table(here("DNA/4.SNV/2.seqtk_fprau/counts.tsv"), header=FALSE, sep="\t")
names(counts) <- c("Sample", "A_Count", "C_Count", "G_Count", "T_Count")

counts <- mutate(counts, SampleCode = gsub("_.*", "", Sample))

# group forward and reverse reads into a single row
counts <- counts %>% group_by(SampleCode) %>% summarise(A_Count = sum(A_Count), 
                                                        T_Count = sum(T_Count), 
                                                        C_Count = sum(C_Count), 
                                                        G_Count = sum(G_Count))

# calculate proportion from counts
counts <- mutate(counts, total=as.numeric(A_Count) + as.numeric(C_Count) + as.numeric(G_Count) + as.numeric(T_Count))
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

leg_labels <- c("A_Proportion", "T_Proportion", "C_Proportion", "G_Proportion")
names(leg_labels) <- c("A Proportion", "T Proportion", "C Proportion", "G Proportion")
p2 <- ggplot(counts, aes(x=Preservative, y=Proportion)) + 
  geom_boxplot(aes(fill=Base), alpha=0.5, width=0.5) + 
  geom_point(aes(fill=Base), shape = 21, color="black", position=position_dodge(width=0.5)) +
  mytheme + 
  ggtitle("*F. prausnitzii* Reads") + 
  theme(plot.title = element_markdown()) + 
  scale_fill_discrete(breaks = leg_labels, labels=names(leg_labels))
p2

ggsave(here("DNA/4.SNV/2.seqtk_fprau/fprau_microbiome_composition.pdf"), w=4.5, h=3, dpi=300)
ggsave(here("DNA/4.SNV/2.seqtk_fprau/fprau_microbiome_composition.jpg"), w=4.5, h=3, dpi=300)

z <- get_legend(p2)
plot_grid(p1 + theme(legend.position = "none"), 
          p2 + theme(legend.position = "none"), 
          z,
          nrow=1, ncol=3, rel_widths = c(1,1,0.4), scale = 0.9, labels = c("a", "b", ""))
ggsave(here("DNA/4.SNV/nucleotide_composition.pdf"), w=8, h=2.75, dpi=300)
ggsave(here("DNA/4.SNV/nucleotide_composition.jpg"), w=8, h=2.75, dpi=300)
