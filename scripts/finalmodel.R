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

# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

raw <- read.csv(here("QSU_Data/summary_model_raw.txt"), sep="\t", header=TRUE)
means <- read.csv(here("QSU_Data/summary_model_means.txt"), sep="\t", header=TRUE)

raw_long <- melt(raw, id.vars = c("Sample", "Patient", "Sample_Type"), variable.name="Microbe")
raw_long <- raw_long %>% mutate(Microbe = ifelse(grepl("Bacteroidetes", Microbe), "Bacteroidetes", "Firmicutes"))

raw_long <- raw_long %>%
  mutate(Sample_Type = fct_relevel(Sample_Type, "NF", "NF_OF_ZF", "NF_OF_OR_ZF_ZR", "NF_OF_OR_OH_ZF_ZR_ZH"))
bf_temps <- c("#933450FF", "#D99755FF")
names(bf_temps) <- c("Bacteroidetes", "Firmicutes")
ggplot(raw_long, aes(x=Sample_Type, y=value*100, color=Microbe)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width=0.3), shape=16, size=1.5) + 
  theme_bw() + 
  geom_errorbar(data=means, inherit.aes = FALSE, aes(x=Sample_Type, ymin=mean_ci_lower*100, ymax=mean_ci_upper*100, fill=Microbe), width=0.1, size=1, position=position_dodge(width=0.8)) +
  geom_point(data=means, inherit.aes = FALSE, aes(x=Sample_Type, y=mean*100, fill=Microbe, guides=NULL), width=0.1, size=2, position=position_dodge(width=0.8), show.legend = FALSE) +
  ylab("Relative Abundance (%)") +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) + 
  scale_x_discrete(labels=c("NF" = "No Preservative\nImmediately Frozen", 
                            "NF_OF_ZF" = "Any Preservative\nImmediately Frozen", 
                            "NF_OF_OR_ZF_ZR" = "Any Preservative\nAny Time before Freezing", 
                            "NF_OF_OR_OH_ZF_ZR_ZH" = "Any Preservative\nAny Time and Temperature\nbefore Freezing")) + 
  scale_color_manual(values=bf_temps)

ggsave(here("QSU_Data/fig5model.pdf"), dpi=300, w=10, h=4.5)
ggsave(here("QSU_Data/fig5model.jpg"), dpi=300, w=10, h=4.5)
