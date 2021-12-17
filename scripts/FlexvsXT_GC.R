library(ggplot2)
library(tidyverse)
library(paletteer) 
# library(ggsignif)
library(here)
# library(reshape2)
library(cowplot)
library(scales)
# library(ggpubr)

readcounts <- read.csv(here("XTvFlex/03_GC/01_subsampling/readcounts_total.tsv"), sep="\t", header = TRUE) # read in count data
readcounts <- mutate(readcounts, Donor = substr(Sample, 1, 3))
readcounts <- mutate(readcounts, Method = ifelse(grepl("merged", Sample, fixed=TRUE), "XT", "Flex"))
readcounts <- mutate(readcounts, Replicate = ifelse(grepl("R1-B", Sample, fixed=TRUE), "B", ifelse(grepl("R1-C", Sample, fixed=TRUE), "C", "A")))
readcounts <- mutate(readcounts, Label = paste(Donor, Method, Replicate))

ggplot(readcounts, aes(x=Label, y=host_removed_reads)) + 
  geom_col(aes(fill=Method)) + 
  theme_bw() + 
  scale_y_continuous(labels = comma) + 
  theme(axis.text.x = element_text(angle=90)) + 
  ylab("Reads") +
  theme(axis.title.x = element_blank())
  
ggsave(here("outputs/figures/XTvFlex_readcounts.pdf"), dpi=300, w=8, h=4)
ggsave(here("outputs/figures/XTvFlex_readcounts.jpg"), dpi=300, w=8, h=4)
