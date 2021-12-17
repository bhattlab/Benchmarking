library(ggplot2)
library(tidyverse)
library(paletteer) 
library(ggsignif)
library(here)

conditionpalette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#97b557","#6b950f")
names(conditionpalette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")
condition_colors <- c("#bababa", paletteer_d("RColorBrewer::PuOr")[c(9,10,11,3,2,1)])
conditionlabels <- c("Immediate Frozen", "OMNI Frozen", "OMNI Room Temp", "OMNI 40C", "Zymo Frozen", "Zymo Room Temp", "Zymo 40C")

dna <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)

dna <- dna %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove = FALSE)
dna <- mutate(dna, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
dna$Condition <- factor(dna$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
dna$Donor <- factor(dna$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))

dna <- mutate(dna, Preservative=ifelse(substr(Condition, 1,1) == "O", "Omnigene", ifelse(substr(Condition, 1, 1) == "Z", "Zymo", ifelse(substr(Condition, 1, 1) == "N", "No Preserative", "Control"))))
dna <- mutate(dna, Temperature=ifelse(substr(Condition, 2,2) == "F", "Frozen", ifelse(substr(Condition, 2,2) == "H", "40C", ifelse(substr(Condition, 2, 2) == "R", "Room Temp", "Control"))))
dna$Preservative <- factor(dna$Preservative, levels = c("No Preserative", "Omnigene", "Zymo"))
dna$Temperature <- factor(dna$Temperature, levels = c("Frozen", "Room Temp", "40C"))


ggplot(dna %>% filter(Condition != "Controls"), aes(x=Donor, y=DNAConcentration, color=Condition, fill=Condition)) + 
  geom_point(alpha=0.8, pch=21, size=2) + 
  theme_classic() +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors)


ggplot(dna %>% filter(Condition != "Controls"), aes(x=Condition, y=DNAConcentration)) + 
  geom_jitter(pch=21, aes(fill=Condition)) + 
  geom_boxplot(alpha =0.80, outlier.shape=NA, aes(fill=Condition)) + 
  theme_classic() + 
  scale_color_manual(values = conditionpalette) +
  scale_fill_manual(values = conditionpalette) + 
  geom_signif(comparisons = list(c("OF", "OH")), map_signif_level = TRUE, y_position=180, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("OF", "OR")), map_signif_level = TRUE, y_position=170, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("OR", "OH")), map_signif_level = TRUE, y_position=160, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("ZF", "ZH")), map_signif_level = TRUE, y_position=60, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("ZF", "ZR")), map_signif_level = TRUE, y_position=50, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("ZR", "ZH")), map_signif_level = TRUE, y_position=40, color="black", tip_length = 0) + 
  guides(color=FALSE, fill=FALSE) + 
  scale_x_discrete(labels = conditionlabels) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust =1), axis.title.x = element_blank())

ggsave(here("outputs/figures/DNAConcentrations_summarybycondition.pdf"), dpi=300, w=4, h=5)




ggplot(dna %>% filter(Condition != "Controls"), aes(x=Donor, y=DNAConcentration, fill=Condition, color=Condition)) + 
  #geom_jitter(alpha=0.8, width=0.2) + 
  geom_point(alpha=0.9, pch=21, size=2) +
  theme_bw() + 
  #scale_fill_manual(values = condition_colors) + 
  #scale_color_manual(values = condition_colors) +
  scale_fill_manual(values=conditionpalette) + 
  scale_color_manual(values=conditionpalette)+
  #facet_wrap(Preservative~Temperature)
  facet_grid(rows=vars(Preservative), cols=vars(Temperature), scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1), axis.title.x = element_blank()) + 
  guides(color=FALSE, fill=FALSE) + 
  ylab("DNA Concentration (ng/ul)")

ggsave(here("outputs/figures/DNAConcentrations_facet.pdf"), dpi=300, w=6, h=4)


#### RNA
rna <- read.csv(here("data/RNAExtraction.tsv"), sep="\t", header=TRUE)
rna <- rna %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove = FALSE)
rna <- mutate(rna, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
rna$Condition <- factor(rna$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
rna$Donor <- factor(rna$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))

rna <- mutate(rna, Preservative=ifelse(substr(Condition, 1,1) == "O", "Omnigene", ifelse(substr(Condition, 1, 1) == "Z", "Zymo", ifelse(substr(Condition, 1, 1) == "N", "No Preserative", "Control"))))
rna <- mutate(rna, Temperature=ifelse(substr(Condition, 2,2) == "F", "Frozen", ifelse(substr(Condition, 2,2) == "H", "40C", ifelse(substr(Condition, 2, 2) == "R", "Room Temp", "Control"))))
rna$Preservative <- factor(rna$Preservative, levels = c("No Preserative", "Omnigene", "Zymo"))
rna$Temperature <- factor(rna$Temperature, levels = c("Frozen", "Room Temp", "40C"))
rna <- mutate(rna, RNAConcentration=ifelse(RNAConcentration == "low", 0, RNAConcentration)) %>% mutate(RNAConcentration=ifelse(RNAConcentration == "TOO LOW", 0, RNAConcentration))
rna$RNAConcentration <- as.numeric(rna$RNAConcentration)

ggplot(rna %>% filter(Condition != "Controls"), aes(x=Donor, y=RNAConcentration, fill=Condition, color=Condition)) + 
  #geom_jitter(alpha=0.8, width=0.2) + 
  geom_point(alpha=0.85, pch=21, size=2) +
  theme_bw() + 
  scale_fill_manual(values = conditionpalette) + 
  scale_color_manual(values = conditionpalette) +
  #facet_wrap(Preservative~Temperature)
  facet_grid(rows=vars(Preservative), cols=vars(Temperature), scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1), axis.title.x = element_blank()) + 
  guides(color=FALSE, fill=FALSE) + 
  ylab("RNA Concentration (ng/ul)")

ggsave(here("outputs/figures/RNAConcentrations_facet.pdf"), dpi=300, w=6, h=4)

ggplot(rna %>% filter(Condition != "Controls"), aes(x=Condition, y=RNAConcentration)) + 
  geom_jitter(pch=21, aes(fill=Condition)) + 
  geom_boxplot(alpha =0.85, outlier.shape=NA, aes(fill=Condition)) + 
  theme_classic() + 
  scale_color_manual(values = conditionpalette) +
  scale_fill_manual(values = conditionpalette) + 
  geom_signif(comparisons = list(c("OF", "OH")), map_signif_level = TRUE, y_position=180, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("OF", "OR")), map_signif_level = TRUE, y_position=170, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("OR", "OH")), map_signif_level = TRUE, y_position=160, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("ZF", "ZH")), map_signif_level = TRUE, y_position=60, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("ZF", "ZR")), map_signif_level = TRUE, y_position=50, color="black", tip_length = 0) +
  geom_signif(comparisons = list(c("ZR", "ZH")), map_signif_level = TRUE, y_position=40, color="black", tip_length = 0) + 
  guides(color=FALSE, fill=FALSE) + 
  scale_x_discrete(labels = conditionlabels) + 
  ylab("RNA Concentration (ng/ul)") +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust =1), axis.title.x = element_blank())

ggsave(here("outputs/figures/RNAConcentrations_summarybycondition.pdf"), dpi=300, w=4, h=5)
