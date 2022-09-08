library(ggplot2)
library(tidyverse)
library(paletteer) 
library(ggsignif)
library(here)
library(ggpubr)


# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

dna <- read.csv(here("QSU_Data/raw_data_dna_full.csv"), head=TRUE) #FIXME WHy only 170 - ten samples have 0 concentration
dna <- dna %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove = FALSE)
dna <- dna %>% mutate(hiddenLabel=ifelse(Condition == "OR", "OMNIgene", ifelse(Condition== "ZR", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples

model <- read.csv(here("QSU_Data/model_data_dna_full.csv"), header=TRUE) # means and confidence intervals for all conditions
model <- model %>% mutate(Condition = Sample_Type)
sig <- read.csv(here("QSU_Data/sig_data_dna_full.tsv"), header=TRUE, sep="\t") # p-values and percent enrichment/depletion for all tests
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values

dna <- mutate(dna, PlotOrder=ifelse(Condition == "NF", 1, 
                                    ifelse(Condition == "OF", 2, 
                                           ifelse(Condition == "OR", 3, 
                                                  ifelse(Condition == "OH", 4,
                                                  ifelse(Condition == "ZF", 5,
                                                         ifelse(Condition == "ZR", 6, 7)))))))


p <- ggplot(dna, aes(x = reorder(Condition, PlotOrder), y=DNAConcentration, fill=Condition)) + 
  geom_jitter(width=0.2, aes(color=Condition), shape=16, size=1.5) +
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) +
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  theme_bw() +
  ylab("DNA Concentration (ng/ul)") +
  geom_errorbar(data=model %>% filter(feature == "DNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "DNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  stat_pvalue_manual(sig %>% filter(feature == "DNAConcentration") %>% filter(group1 == "NF" | group1 == "OF" | group1 == "ZF") %>% filter(p.adj <= 0.05),
                     tip.length=0, label = "p.signif", inherit.aes = FALSE, y.position = c(160, 140, 150)) +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p
b <- ggplot(dna %>% filter(Donor == "D01" & Replicate == "R1"), aes(x=reorder(Condition, PlotOrder), y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

dna_plot <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='lr')
dna_plot

ggsave(here("QSU_Data/Supplement_DNAConcentration.jpeg"), dpi=300, w=5, h=5)
ggsave(here("QSU_Data/Supplement_DNAConcentration.pdf"), dpi=300, w=5, h=5)

# RNA

rna <- read.csv(here("QSU_Data/raw_data_rna_full.csv"), head=TRUE) #FIXME WHy only 170 - ten samples have 0 concentration
rna <- rna %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove = FALSE)
rna <- rna %>% mutate(hiddenLabel=ifelse(Condition == "OR", "OMNIgene", ifelse(Condition== "ZR", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples

model <- read.csv(here("QSU_Data/model_data_rna_full.csv"), header=TRUE) # means and confidence intervals for all conditions
model <- model %>% mutate(Condition = Sample_Type)
sig <- read.csv(here("QSU_Data/sig_data_rna_full.tsv"), header=TRUE, sep="\t") # p-values and percent enrichment/depletion for all tests
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values

rna <- mutate(rna, PlotOrder=ifelse(Condition == "NF", 1, 
                                                        ifelse(Condition == "OF", 2, 
                                                               ifelse(Condition == "OR", 3, 
                                                                             ifelse(Condition == "ZF", 4,
                                                                                    ifelse(Condition == "ZR", 5, 6))))))

p <- ggplot(rna, aes(x = reorder(Condition, PlotOrder), y=RNAConcentration)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=3.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Condition), shape=16, size=1.5) + 
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  theme_bw() + 
  ylab("RNA Concentration (ng/ul)") + 
  ylim(0, 350) +
  geom_errorbar(data=model %>% filter(feature == "RNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "RNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  stat_pvalue_manual(sig %>% filter(feature == "RNAConcentration") %>% filter(p.adj <= 0.05), y.position = c(325, 350, 275, 275),
                     tip.length=0, label = "p.signif") +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p

foodf <- data.frame(xvals = c(0.3, 2, 4.6), labels=c("None", "OMNIgene", "Zymo"))

test <- ggplot(foodf, aes(x=xvals, y=0)) + 
  geom_text(aes(y=0, label=labels), fontface = "bold") + 
  ylim(-0.5, 0.5) +
  xlim(0,6) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme_void()
test

# b <- ggplot(rna %>% filter(Donor == "D01" & Replicate == "R1"), aes(x=Condition, y=0)) + 
#   geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
#   ylim(-0.5, 0.5) +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# b

rna_plot <- plot_grid(p,test, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='lr')
rna_plot

ggsave(here("QSU_Data/Supplement_RNAConcentration.jpeg"), dpi=300, w=5, h=5)
ggsave(here("QSU_Data/Supplement_RNAConcentration.pdf"), dpi=300, w=5, h=5)

plot_grid(dna_plot, rna_plot, ncol=2, nrow=1, scale=0.9, labels=c("a", "b"))

ggsave(here("QSU_Data/Supplement_Concentration.jpeg"), dpi=300, w=10, h=5)
ggsave(here("QSU_Data/Supplement_Concentration.pdf"), dpi=300, w=10, h=5)
