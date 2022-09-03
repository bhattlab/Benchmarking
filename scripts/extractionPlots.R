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

# dna <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
# 
# dna <- dna %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove = FALSE)
# dna <- mutate(dna, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
# dna$Condition <- factor(dna$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# dna$Donor <- factor(dna$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
# 
# dna <- mutate(dna, Preservative=ifelse(substr(Condition, 1,1) == "O", "Omnigene", ifelse(substr(Condition, 1, 1) == "Z", "Zymo", ifelse(substr(Condition, 1, 1) == "N", "No Preserative", "Control"))))
# dna <- mutate(dna, Temperature=ifelse(substr(Condition, 2,2) == "F", "Frozen", ifelse(substr(Condition, 2,2) == "H", "40C", ifelse(substr(Condition, 2, 2) == "R", "Room Temp", "Control"))))
# dna$Preservative <- factor(dna$Preservative, levels = c("No Preserative", "Omnigene", "Zymo"))
# dna$Temperature <- factor(dna$Temperature, levels = c("Frozen", "Room Temp", "40C"))
# dna <- dna %>% mutate(hiddenLabel=ifelse(Condition == "OR", "OMNIgene", ifelse(Condition== "ZR", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples
# 


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
                                                  ifelse(Condition == "ZF", 4,
                                                         ifelse(Condition == "ZR", 5, 6))))))


p <- ggplot(dna %>% filter(Condition != "Controls"), aes(x = Condition, y=DNAConcentration, fill=Condition)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Condition), shape=16, size=1.5) + 
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  theme_bw() +
  ylab("DNA Concentration (ng/ul)") + 
  geom_errorbar(data=model %>% filter(feature == "DNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "DNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  stat_pvalue_manual(sig %>% filter(feature == "DNAConcentration") %>% filter(group1 == "NF" | group1 == "OF" | group1 == "ZF") %>% filter(p.adj <= 0.05),
                     tip.length=0, label = "p.signif") 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p
b <- ggplot(dna %>% filter(Donor == "D01" & Replicate == "R1"), aes(x=Condition, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

conc <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='lr')
conc

ggsave(here("QSU_Data/Supplement_DNAConcentration.jpeg"), dpi=300, w=5, h=5)
ggsave(here("QSU_Data/Supplement_DNAConcentration.pdf"), dpi=300, w=5, h=5)

# RNA

# rna <- read.csv(here("data/RNAExtraction.tsv"), sep="\t", header=TRUE)
# 
# rna <- rna %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove = FALSE)
# rna <- mutate(rna, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
# rna$Condition <- factor(rna$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# rna$Donor <- factor(rna$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
# 
# rna <- mutate(rna, Preservative=ifelse(substr(Condition, 1,1) == "O", "Omnigene", ifelse(substr(Condition, 1, 1) == "Z", "Zymo", ifelse(substr(Condition, 1, 1) == "N", "No Preserative", "Control"))))
# rna <- mutate(rna, Temperature=ifelse(substr(Condition, 2,2) == "F", "Frozen", ifelse(substr(Condition, 2,2) == "H", "40C", ifelse(substr(Condition, 2, 2) == "R", "Room Temp", "Control"))))
# rna$Preservative <- factor(rna$Preservative, levels = c("No Preserative", "Omnigene", "Zymo"))
# rna$Temperature <- factor(rna$Temperature, levels = c("Frozen", "Room Temp", "40C"))
# rna <- rna %>% mutate(hiddenLabel=ifelse(Condition == "OR", "OMNIgene", ifelse(Condition== "ZR", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples
# 
# rna <- mutate(rna, RNAConcentration=ifelse(RNAConcentration == "low", 0, RNAConcentration)) %>% mutate(RNAConcentration=ifelse(RNAConcentration == "TOO LOW", 0, RNAConcentration))
# rna$RNAConcentration <- as.numeric(rna$RNAConcentration)
# rna <- rna %>% filter(Condition != "Controls") %>% filter(Condition != "OH")
# rna <- rna %>% mutate(SampleID = gsub("_", "-", SampleID))
# test1 <- rna$SampleID

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
  stat_pvalue_manual(sig %>% filter(feature == "RNAConcentration") %>% filter(p.adj <= 0.05), y.position = c(325, 350, 150, 150),
                     tip.length=0, label = "p.signif") +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p

b <- ggplot(rna %>% filter(Donor == "D01" & Replicate == "R1"), aes(x=Condition, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

conc <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='lr')
conc

ggsave(here("QSU_Data/Supplement_RNAConcentration.jpeg"), dpi=300, w=5, h=5)
ggsave(here("QSU_Data/Supplement_RNAConcentration.pdf"), dpi=300, w=5, h=5)

##### END
