library(cowplot)
library(ggpubr)
library(gtools)
library(here)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(ggsignif)
library(vegan)
library(ggnewscale)
library(scales)


# Read in dataframe
precounts <- read.table(here("DNA/1.preprocess/preprocessing/01_processing/readcounts.tsv"), sep = "\t", header = T)

#Read pairs - multiply by 2
#Filter out controls
counts <- precounts %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)
counts <- counts %>%
  filter((Donor != "NCO" & Donor != "PCO")) %>%
  mutate(raw_reads = raw_reads * 2,
         trimmed_reads = trimmed_reads * 2,
         dedup_reads = dedup_reads * 2,
         host_removed_reads = host_removed_reads * 2,
         orphan_reads = orphan_reads * 2
  )

#Calculate # of human reads
counts <- mutate(counts, host_removed_reads=as.numeric(host_removed_reads))
counts <- mutate(counts, trimmed_reads=as.numeric(trimmed_reads))
dylcount <- counts %>%
  mutate(
    fracmicrobe = host_removed_reads/trimmed_reads,
    frachuman = 1-fracmicrobe
  )

hucount <- counts %>%
  mutate(
    human_removed = trimmed_reads - host_removed_reads,
    human_rm_frac = human_removed / trimmed_reads
  )

#Create plot order of conditions
counts <- mutate(hucount, PlotOrder=ifelse(Condition == "NF", 1, 
                                       ifelse(Condition == "OF", 2, 
                                              ifelse(Condition == "OR", 3, 
                                                     ifelse(Condition == "OH", 4, 
                                                            ifelse(Condition == "ZF", 5,
                                                                   ifelse(Condition == "ZR", 6, 7)))))))
#Color condition by
condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# Plot Fraction Human Reads Removed
ggplot(counts %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=human_rm_frac * 100)) +
  geom_jitter(alpha = 0.75, width = 0.3, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, aes(fill = Condition)) +
  scale_fill_manual(values = condition_palette) +
  theme_cowplot() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "% Human reads removed\nafter de-duplication"
  ) +
  background_grid(major = "y") 

ggsave(here("outputs/figures/hureadsrm_bycondition.jpg"), dpi=300, w=5, h=6)

