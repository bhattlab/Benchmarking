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

##### SET UP #####
# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# read in the QSU data as three separate dataframes 
raw <- read.csv(here("QSU_Data/raw_data_rna_all.csv"), header=TRUE) # raw per sample information
model <- read.csv(here("QSU_Data/raw_data_rna_model_means.csv"), header=TRUE) # means and confidence intervals for all conditions
sig <- read.csv(here("QSU_Data/raw_data_rna_significance.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

# format and edit dataframes 
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values
raw <- raw %>% mutate(Condition = Sample_Type)
raw <- raw %>% mutate(Preservative = substr(Sample_Type, 1, 1)) # make preservative column (first character of sample code)
raw <- raw %>% mutate(Temperature = substr(Sample_Type, 2, 2)) # make temperature column (second character of sample code)
raw <- raw %>% mutate(Sample_Type = fct_relevel(Sample_Type, "NF", "OF", "OR", "ZF", "ZR", "ZH")) # force the ordering of the conditions
raw <- raw %>% mutate(TemperatureLong = ifelse(Temperature == "F", "-80C", ifelse(Temperature == "R", "23C", "40C"))) # write out temperature
raw <- raw %>% mutate(PreservativeLong = ifelse(Preservative == "N", "None", ifelse(Preservative == "O", "Omnigene", "Zymo"))) # write out preservative
raw <- raw %>% mutate(Label = paste(TemperatureLong, PreservativeLong, sep="\n")) # make a label of temperature and preserative
raw <- raw %>% mutate(hiddenLabel=ifelse(Sample_Type == "OR", "OMNIgene", ifelse(Sample_Type== "ZR", "Zymo", ifelse(Sample_Type == "NF", "None", "")))) # make a label just for the "R" samples
model <- model %>% mutate(Condition=Sample_Type)
colnames(model)[colnames(model) == "prediction"] <- "Mean"


##### TAXONOMIC RELATIVE ABUNDANCE STACKED BAR #####

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

genus <- read.csv(here("RNA/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
phylum_plottingorder <- read.csv(here("QSU_Data/genus_to_phylum_RNA.csv"), header=TRUE)
#Color palette
n_taxa <- 15
myCols <- colorRampPalette(brewer.pal(9, "Set1")) # WAS Set1
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

abundance_threshold <- sort(rowSums(genus), decreasing = T)[n_taxa]
bracken_plot <- genus[rowSums(genus) >= abundance_threshold,]
bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  100 - colSums(bracken_plot))))

bracken_plot$Genus <- row.names(bracken_plot)
bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", bracken_plot$Genus)
bracken_long <- melt(bracken_plot, id.vars = "Genus", variable.name = "Sample", value.name = "rel_abundance")
bracken_long <- mutate(bracken_long, Sample=gsub("\\.", "_", Sample))

# Merge in the metadata
colnames(metadata)[3]<-"Sample"
bracken_pheno <- merge(bracken_long, metadata, by = "Sample")

# Correct the plotting order
bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)

samplabels <- bracken_pheno$TemperatureLong
names(samplabels) <- bracken_pheno$Sample

bracken_pheno <- mutate(bracken_pheno, PlotOrder=ifelse(Condition == "NF", 1, 
                                                        ifelse(Condition == "OF", 2, 
                                                               ifelse(Condition == "OR", 3, 
                                                                      ifelse(Condition == "OH", 4, 
                                                                             ifelse(Condition == "ZF", 5,
                                                                                    ifelse(Condition == "ZR", 6, 7)))))))

#Filter out controls and other donors
bracken_pheno <- bracken_pheno %>% filter(Donor != "NCO" & Donor != "PCO")
bracken_pheno <- bracken_pheno %>% filter(Donor=="D01")
bracken_pheno <- mutate(bracken_pheno, Temperature=substr(Condition, 2,2))
bracken_pheno <- mutate(bracken_pheno, Genus=gsub("unclassified ", "", Genus))
bracken_pheno <- merge(bracken_pheno, phylum_plottingorder, by="Genus", all.x=TRUE)
bracken_pheno$Genus <- reorder(bracken_pheno$Genus, bracken_pheno$PlotOrder.y)
barplot_pal <- c("#7BBD5D", "#91C74E", "#A7D03E", "#D9C634", "#E8DB41", "#D12E2E", "#D95034", "#E0723B", "#E89441", "#2E3CD1","#3C3ED6","#4B41DA","#5943DF","#6846E3","#7648E8","#844AED","#934DF1","#A14FF6","#B052FA","#BE54FF", "#939393" )

bluegreens <- c("#dff2f1","#b2dfdb","#80cbc4","#4db6ac","#26a59a","#009788","#01887b","#00796b","#00695c","#004d40")
blues <- c("#e3f2fe","#bbdefb","#90caf9","#64b5f7","#42a5f5","#2096f3","#1f88e5","#1a76d2","#1665c0","#0d47a1", "#082b71")
purples <- c("#EDE7F6","#D1C4E9","#B39EDB","#9575CD","#7E58C2","#673AB7","#5E34B1","#512DA8","#4527A0","#311B92")
pinks <- c("#FDE4EC","#F8BBD0","#F48FB1","#F06293","#EB3F7A","#E91E63","#D81A60","#C2185B","#AD1456","#880E4F")
oranges <- c("#FFF3E0","#FFE0B2","#FFCC80","#FFB74D","#FFA726","#FF9801","#FB8B00","#F57C01","#EF6C00","#E65100")
# 1 2 4 8
barplot_pal <- c(bluegreens[3], rev(purples[3:4]), rev(oranges[2:5]), rev(blues[1:8]), "#939393")

r <-ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder.x), y=rel_abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  labs(
    x = "",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual("Genus", values = barplot_pal) +
  guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
  theme_bw() +
  scale_x_discrete(labels = samplabels) +
  theme(
    plot.title = element_text(face = "plain", size = 14),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
    strip.text = element_text(color = "black", size = 12)) + 
  scale_y_continuous(limits = c(-5, 100.1), expand = c(0, 0)) +
  new_scale_fill() +
  geom_tile(aes(x=Sample, y = -2, fill = Condition), show.legend = F) + 
  geom_tile(aes(x=Sample, y = -3, fill = Condition), show.legend = F) + 
  geom_tile(aes(x=Sample, y = -4, fill = Condition), show.legend = F) +
  scale_fill_manual(values = condition_palette) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
r

b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.05, 0.05) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ 
  theme(text = element_text(size=10))
b

stackedbargenus <- plot_grid(r,b, nrow=2, ncol=1, rel_heights=c(1, 0.1), align="v", axis='lr')
stackedbargenus

ggsave(here("QSU_Data/Figure4A.pdf"), dpi=300, h=4, w=6)
ggsave(here("QSU_Data/Figure4A.jpeg"), dpi=300, h=4, w=6)


##### SHANNON DIVERSITY #####
#Change column name
colnames(raw)[colnames(raw) == "Shannon.Entropy"] <- "Shannon"

#Shannon entropy
se <- ggplot(raw , aes(Sample_Type, Shannon)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=3.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Shannon Entropy"), inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Shannon Entropy"), inherit.aes=FALSE, aes(x=Condition, y=Mean), size=2.5) +
  ylim(0,4.5) +
  stat_pvalue_manual(sig %>% filter(feature == "Shannon Entropy") %>% filter(p.adj <= 0.05), y.position=c(4),
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Shannon Entropy") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
se

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") +
  ylim(-0.1, 0.1) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"))
b
shannon <- plot_grid(se,b, nrow=2, ncol=1, rel_heights=c(1,0.05), align="v", axis='l')
shannon

ggsave(here("QSU_Data/Fig4B_Shannon.jpeg"), dpi=300, width=5, h=5)
ggsave(here("QSU_Data/Fig4B_Shannon.pdf"), dpi=300, width=5, h=5)

##### RICHNESS #####
p <- ggplot(raw, aes(x = Sample_Type, y=Richness.0.01.)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=3.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type), shape=16, size=1.5) + 
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, y=Mean), size=2.5) +
  stat_pvalue_manual(sig %>% filter(feature == "Richness 0.01%") %>% filter(p.adj <= 0.05), y.position=c(350),
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Number of Genera above \n0.01% Abundance") + 
  ylim(0,400) + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

richness <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='l')
richness

ggsave(here("QSU_Data/Fig4C_Richness.jpeg"), dpi=300, w=5, h=5)
ggsave(here("QSU_Data/Fig4C_Richness.pdf"), dpi=300, w=5, h=5)

##### PHYLUM FOREST #####
raw <- mutate(raw, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                    ifelse(Sample_Type == "OF", 2, 
                                           ifelse(Sample_Type == "OR", 3, 
                                                         ifelse(Sample_Type == "ZF", 4,
                                                                ifelse(Sample_Type == "ZR", 5, 6))))))
model <- mutate(model, PlotOrder=ifelse(Condition == "NF", 1, 
                                        ifelse(Condition == "OF", 2, 
                                               ifelse(Condition == "OR", 3, 
                                                             ifelse(Condition == "ZF", 4,
                                                                    ifelse(Condition == "ZR", 5, 6))))))

bact <- ggplot(raw , aes(x=reorder(Sample_Type, PlotOrder), Relative.Abundance..Bacteroidetes*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Bacteroidetes"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Bacteroidetes"), inherit.aes=FALSE, aes(x=Condition, y=Mean*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Bacteroidetes") %>% filter(p.adj <= 0.05), y.position=c(90,95, 85),
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylim(0,100) + 
  ylab("Relative Abundance (%)") + 
  xlab("Bacteroidetes") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), axis.text.x = element_blank(), axis.ticks.x = element_blank())
bact

firm <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Firmicutes*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Firmicutes"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Firmicutes"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), y=Mean*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Firmicutes") %>% filter(p.adj <= 0.05), y.position=c(90, 95, 70,75), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylim(0,100) + 
  ylab("Relative Abundance") + 
  xlab("Firmicutes") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y=element_blank())
firm

act <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Actinobacteria*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Actinobacteria"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Actinobacteria"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), y=Mean*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Actinobacteria") %>% filter(p.adj <= 0.05), y.position=c(27, 22, 17), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylim(0,100) + 
  ylab("Relative Abundance") + 
  xlab("Actinobacteria") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y=element_blank())
act

vir <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Viruses*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Viruses"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Viruses"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), y=Mean*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Viruses") %>% filter(p.adj <= 0.05), y.position=c(14, 15, 13, 13), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylim(0,15) + 
  ylab("Relative Abundance") + 
  xlab("Viruses") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())
vir

fungi <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Fungi*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Fungi"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Fungi"), inherit.aes=FALSE, aes(x=reorder(Condition, PlotOrder), y=Mean*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Fungi") %>% filter(p.adj <= 0.05), y.position=c(13, 12, 11, 12), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylim(0,15) + 
  ylab("Relative Abundance") + 
  xlab("Fungi") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y=element_blank())
fungi

pdiff <- plot_grid(bact, firm, act, vir, fungi, nrow=1, ncol=5, rel_widths = c(1.3, 1, 1, 1.1, 1))
pdiff

ggsave(here("QSU_Data/Fig4_PhylumEnrichmentDraft.pdf"), dpi=300, h=3.2, w=4)
ggsave(here("QSU_Data/Fig4_PhylumEnrichmentDraft.jpeg"), dpi=300, h=3.2, w=4)


##### BRAY CURTIS #####
bc_raw <- read.csv(here("QSU_Data/raw_data_rna_all_braycurtis.csv"), header=TRUE)
bc_raw <- bc_raw %>% mutate(UniqueComparison = ifelse(Protocol_1 > Protocol_2, paste(Patient, Protocol_1, Protocol_2, sep="_"), paste(Patient, Protocol_2, Protocol_1, sep="_")))
bc_raw <- bc_raw %>% select(Patient, bcdist, UniqueComparison)
bc_unique <- distinct(bc_raw, UniqueComparison, bcdist, .keep_all = TRUE)

bc_model <- read.csv(here("QSU_Data/raw_data_rna_model_means_braycurtis.csv"), header=TRUE)
bc_model <- bc_model %>% separate(Protocols, c("group1", "group2"), remove=FALSE) %>% filter(group1 == "NF") %>% filter(group1 != group2)
bc_model <- bc_model %>% mutate(Condition=group2)

bc_sig <- read.csv(here("QSU_Data/raw_data_rna_significance_braycurtis.csv"), header=TRUE)
bc_sig_toNF <- bc_sig %>% filter(feature == "BrayCurtisBetweenConditions")
bc_sig_toNF <- bc_sig_toNF %>% mutate(group1=gsub("_NF", "", group1)) %>% mutate(group2=gsub("_NF", "", group2))

# across condition Bray Curtis
bc_across <- bc_unique %>% separate(UniqueComparison, c("Donor", "group1", "group2"), remove=FALSE)
bc_across <- bc_across %>% filter(group1 != group2) %>% filter(group1 == "NF" | group2 == "NF")
bc_across <- bc_across %>% mutate(modelID = paste(group2, group1, sep="_"))
bc_across <- bc_across %>% mutate(group1 = fct_relevel(group1,  "OF", "OR", "ZF", "ZR", "ZH"))

bc_plot <- ggplot(bc_across, aes(x=group1, y=bcdist)) + 
  geom_vline(aes(xintercept=2.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=group1),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=estimate), size=2.5) +
  ylim(0,1) +
  stat_pvalue_manual(bc_sig_toNF %>% filter(p.adj <= 0.05), y.position=c(0.95, .85, .9), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Bray-Curtis Dissimilarity") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
bc_plot

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1") %>% filter(Sample_Type != "NF"), aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") +
  ylim(-0.05, 0.05) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"))
b
bc_full <- plot_grid(bc_plot,b, nrow=2, ncol=1, rel_heights=c(1,0.05), align="v", axis='lr')
bc_full

ggsave(here("QSU_Data/Fig4_BrayCurtis.jpeg"), dpi=300, w=5, h=5)
ggsave(here("QSU_Data/Fig4_BrayCurtis.pdf"), dpi=300, w=5, h=5)

##### HEATMAP #####

sig <- sig %>% mutate(PercentFormatted = as.numeric(gsub("%.*", "", percentchange)))
sig <- sig %>% mutate(Condition = group2)
sig <- sig %>% mutate(y.position = ifelse(Condition == "OF", 2, 
                                                      ifelse(Condition=="ZF", 1, 
                                                             ifelse(Condition=="OR", 3, 
                                                                           ifelse(Condition == "ZR", 2, 1)))))

genus_sig <- sig %>% filter(grepl(":", feature)) %>% filter(!grepl("Relative", feature))
raw_abundance <- model %>% filter(grepl(":", feature)) %>% filter(!grepl("Relative", feature))
raw_abundance <- raw_abundance %>% filter(Condition == "NF") %>% select(feature, Mean)
names(raw_abundance) <- c("feature", "NFMean")
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: unclassified Clostridiales.*", "Firmicutes: unclassified Clostridiales (miscel...", feature))
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: unclassified Firmicutes sensu.*", "Firmicutes: unclassified Firmicutes sensu stri...", feature))
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Bacteroidetes: Prevotellaceae: environmental.*", "Bacteroidetes: Prevotellaceae: environmental s...", feature))
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: Lachnospiraceae: environmental.*", "Firmicutes: Lachnospiraceae: environmental sam...", feature))


raw_abundance <- raw_abundance %>% mutate(Phylum = gsub(":.*", "", feature))
raw_abundance <- raw_abundance %>% 
  group_by(Phylum) %>% 
  arrange(Phylum, desc(NFMean))
raw_abundance$x.position <- as.numeric(row.names(raw_abundance))


genus_sig <- merge(genus_sig, raw_abundance, by="feature", all.x = TRUE)
genus_sig <- genus_sig %>% mutate(Phylum = gsub(":.*", "", feature))
genus_sig <- genus_sig %>% mutate(Phylum = ifelse(Phylum == "null", "Virus", Phylum))
genus_sig <- genus_sig %>% mutate(Plabel = ifelse(p.format <= 0.05, "*", ""))
genus_sig <- genus_sig %>% mutate(Genus = gsub(" sensu.*", "", gsub("unclassified ", "other ", gsub(".*: ", "", gsub(" \\(misc.*", "", gsub(": environmental.*", "", feature))))))
genus_sig <- genus_sig %>% mutate(PercentFormatted2 = ifelse(PercentFormatted > 100, 100, PercentFormatted))
preservative <- ggplot(genus_sig %>% filter(Condition == "OF" | Condition == "ZF"), aes(x=x.position, y=y.position, fill=PercentFormatted2)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1, 2), labels=c("Zymo -80°C", "OMNI -80°C")) + 
  labs(title = "Preservative Effect (Relative to no preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-100, 100)) +
  labs(fill="Percent\nEnrichment")+
  geom_point(data = genus_sig %>% filter(Condition == "OF" | Condition == "ZF") %>% filter(Plabel == "*"), size=1, color="white") + 
  theme(panel.border = element_rect(size=1)) + 
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=12.5) +
  geom_vline(xintercept=13.5) +
  geom_vline(xintercept=44.5) +
  geom_vline(xintercept=45.5) +
  geom_vline(xintercept=47.5) +
  geom_vline(xintercept=50.5) 
preservative

z <- get_legend(preservative)

temperature <- ggplot(genus_sig %>% filter(Condition == "OR" | Condition == "ZH" | Condition == "ZR"), aes(x=x.position, y=y.position, fill=PercentFormatted2)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0), breaks=genus_sig$x.position, labels=genus_sig$Genus) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1,2,3), labels=c("Zymo 40°C", "Zymo 23°C", "OMNI 23°C")) + 
  labs(title = "Temperature Effect (Relative to -80°C for each preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(),
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-100, 100)) +
  geom_point(data = genus_sig %>% filter(Condition == "OR" | Condition == "ZH" | Condition == "ZR") %>% filter(Plabel == "*"), size=1, color="white") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ 
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=12.5) +
  geom_vline(xintercept=13.5) +
  geom_vline(xintercept=44.5) +
  geom_vline(xintercept=45.5) +
  geom_vline(xintercept=47.5) +
  geom_vline(xintercept=50.5) +
  theme(panel.border = element_rect(size=1))
temperature

raw_abund <- ggplot(genus_sig %>% filter(Condition == "ZF"), aes(x=x.position, y=y.position, fill=NFMean)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1), labels=c("Abundance")) + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("#F3F3F3", "black")) +
  labs(fill="Baseline\nAbundance")+ 
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=12.5) +
  geom_vline(xintercept=13.5) +
  geom_vline(xintercept=44.5) +
  geom_vline(xintercept=45.5) +
  geom_vline(xintercept=47.5) +
  geom_vline(xintercept=50.5) +
  theme(panel.border = element_rect(size=1))
raw_abund

z2 <- get_legend(raw_abund)

phylum_scale <- ggplot(genus_sig %>% filter(Condition == "ZF"), aes(x=x.position, y=y.position, fill=Phylum)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1), labels=c("Phylum")) + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(fill="Phylum") + 
  geom_vline(xintercept=4.5) +
  geom_vline(xintercept=12.5) +
  geom_vline(xintercept=13.5) +
  geom_vline(xintercept=44.5) +
  geom_vline(xintercept=45.5) +
  geom_vline(xintercept=47.5) +
  geom_vline(xintercept=50.5) +
  theme(panel.border = element_rect(size=1))
phylum_scale 

z3 <- get_legend(phylum_scale)

main <- plot_grid(phylum_scale + theme(legend.position = "none"),
          raw_abund + theme(legend.position ="none"),  
          preservative + theme(legend.position = "none"), 
          temperature + theme(legend.position = "none"), 
          ncol=1, nrow=4, rel_widths = c(1,1,1,1), rel_heights = c(1, 1, 1, 3),
          align="v", axis="l")

plot_grid(main, z3, z2, z, nrow=1, ncol=4, rel_widths = c(1,.1,.1,.1), scale=0.9)

ggsave(here("QSU_Data/Fig4_Heatmap.pdf"), dpi=300, w=11, h=4)
ggsave(here("QSU_Data/Fig4_Heatmap.jpeg"), dpi=300, w=11, h=4)


#### SUPPLEMENTS ####

#3A: Stacked bar plot for Donor 1
# read in metadata
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

genus <- read.csv(here("RNA/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
phylum_plottingorder <- read.csv(here("QSU_Data/genus_to_phylum_RNA.csv"), header=TRUE)
#Color palette
n_taxa <- 15
myCols <- colorRampPalette(brewer.pal(9, "Set1")) # WAS Set1
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

abundance_threshold <- sort(rowSums(genus), decreasing = T)[n_taxa]
bracken_plot <- genus[rowSums(genus) >= abundance_threshold,]
bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  100 - colSums(bracken_plot))))

bracken_plot$Genus <- row.names(bracken_plot)
bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", bracken_plot$Genus)
bracken_long <- melt(bracken_plot, id.vars = "Genus", variable.name = "Sample", value.name = "rel_abundance")
bracken_long <- mutate(bracken_long, Sample=gsub("\\.", "_", Sample))

# Merge in the metadata
colnames(metadata)[3]<-"Sample"
bracken_pheno <- merge(bracken_long, metadata, by = "Sample")
#bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))

# Correct the plotting order
bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)

samplabels <- bracken_pheno$TemperatureLong
names(samplabels) <- bracken_pheno$Sample

bracken_pheno <- mutate(bracken_pheno, PlotOrder=ifelse(Condition == "NF", 1, 
                                                        ifelse(Condition == "OF", 2, 
                                                               ifelse(Condition == "OR", 3, 
                                                                      ifelse(Condition == "OH", 4, 
                                                                             ifelse(Condition == "ZF", 5,
                                                                                    ifelse(Condition == "ZR", 6, 7)))))))

#Filter out controls and other donors
bracken_pheno <- bracken_pheno %>% filter(Donor != "NCO" & Donor != "PCO")
bracken_pheno <- mutate(bracken_pheno, Temperature=substr(Condition, 2,2))
bracken_pheno <- mutate(bracken_pheno, Genus=gsub("unclassified ", "", Genus))
bracken_pheno <- merge(bracken_pheno, phylum_plottingorder, by="Genus", all.x=TRUE)
bracken_pheno$Genus <- reorder(bracken_pheno$Genus, bracken_pheno$PlotOrder.y)
barplot_pal <- c("#7BBD5D", "#91C74E", "#A7D03E", "#D9C634", "#E8DB41", "#D12E2E", "#D95034", "#E0723B", "#E89441", "#2E3CD1","#3C3ED6","#4B41DA","#5943DF","#6846E3","#7648E8","#844AED","#934DF1","#A14FF6","#B052FA","#BE54FF", "#939393" )

bracken_pheno <- mutate(bracken_pheno, DonorLabel=ifelse(Donor == "D10", "Donor 10", paste("Donor", substr(Donor, 3,3))))
bracken_pheno$DonorLabel <- factor(bracken_pheno$DonorLabel, levels = c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6", "Donor 7", "Donor 8", "Donor 9", "Donor 10"))
bracken_pheno <- bracken_pheno %>% mutate(PreservationLong = ifelse(Preservation == "N", "No Preservative", ifelse(Preservation == "O", "OMNI", "Zymo")))
bracken_pheno <- bracken_pheno %>% mutate(TemperatureLongAll = ifelse(Temperature == "F", "-80°C", ifelse(Temperature == "R", "23°C", "40°C"))) # write out temperature
bracken_pheno <- bracken_pheno %>% mutate(LegendLabel = paste(PreservationLong, TemperatureLongAll, sep=" "))

bluegreens <- c("#dff2f1","#b2dfdb","#80cbc4","#4db6ac","#26a59a","#009788","#01887b","#00796b","#00695c","#004d40")
blues <- c("#e3f2fe","#bbdefb","#90caf9","#64b5f7","#42a5f5","#2096f3","#1f88e5","#1a76d2","#1665c0","#0d47a1", "#082b71")
purples <- c("#EDE7F6","#D1C4E9","#B39EDB","#9575CD","#7E58C2","#673AB7","#5E34B1","#512DA8","#4527A0","#311B92")
pinks <- c("#FDE4EC","#F8BBD0","#F48FB1","#F06293","#EB3F7A","#E91E63","#D81A60","#C2185B","#AD1456","#880E4F")
oranges <- c("#FFF3E0","#FFE0B2","#FFCC80","#FFB74D","#FFA726","#FF9801","#FB8B00","#F57C01","#EF6C00","#E65100")
barplot_pal <- c(bluegreens[3], rev(purples[3:4]), rev(oranges[2:5]), rev(blues[1:8]), "#939393")

r <-ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder.x), y=rel_abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  labs(
    x = "",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual("Genus", values = barplot_pal) +
  guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
  theme_bw() +
  scale_x_discrete(labels = NULL) +
  theme(
    plot.title = element_text(face = "plain", size = 14),
    legend.text = element_text(size = 10),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
    strip.text = element_text(color = "black", size = 12)) + 
  scale_y_continuous(limits = c(-5, 100.1), expand = c(0, 0)) +
  facet_wrap(~DonorLabel, ncol = 5, scales = "free") + 
  new_scale_fill() +
  geom_tile(aes(x=Sample, y = -2, fill = Condition), show.legend = F) + 
  geom_tile(aes(x=Sample, y = -3, fill = Condition), show.legend = F) + 
  geom_tile(aes(x=Sample, y = -4, fill = Condition), show.legend = F) +
  scale_fill_manual(values = condition_palette) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
r

a <- get_legend(r)
condition_labels2 <- c("No Preservative -80°C","OMNI -80°C","OMNI 23°C","Zymo -80°C","Zymo 23°C","Zymo 40°C")
names(condition_labels2) <- c( "NF","OF", "OR", "ZF", "ZR", "ZH")
condition_palette <- c("#acaaaf","#c3a5d0","#9a6faa","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("NF", "OF", "OR", "ZF", "ZR", "ZH")

dummy <- ggplot(bracken_pheno %>% filter(Condition != "OH"), aes(x=reorder(Sample, PlotOrder.x), y=rel_abundance)) +
  geom_bar(stat="identity", aes(fill=Condition)) + 
  scale_fill_manual(values=condition_palette, labels=condition_labels2) + 
  theme(legend.direction="horizontal", legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  theme(text = element_text(size=10))
dummy
b <- get_legend(dummy)
stacked_bar <- plot_grid(r, b, nrow=2, ncol=1, rel_heights = c(1, 0.08), rel_widths = c(1, 1))
stacked_bar
ggsave(here("QSU_Data/Supplement_RNA_stackedbar.pdf"), dpi=300, h=6, w=15)
ggsave(here("QSU_Data/Supplement_RNA_stackedbar.jpeg"), dpi=300, h=6, w=15)

#Plot Figure 4
# a <- plot_grid(stacked_bar,shannon,nrow=1,ncol=2,scale=0.9,rel_widths=c(1, 0.3),labels=c("a","b"), align = "v")
# a
# b <- plot_grid(richness,bc_full,pdiff,nrow=1,ncol=3, scale=0.9, rel_widths=c(0.7, 0.6,1.1),labels=c("c","d","e"), align="v", axis="t")
# b
# three<-plot_grid(a, b, nrow=2, ncol=1, rel_widths=c(1, 1), align="h", axis="lr")
# three
a <- plot_grid(stacked_bar,pdiff,nrow=1,ncol=2,scale=0.95,rel_widths=c(1.5,1.1),labels=c("a","b"), align = "v")
a
b <- plot_grid(shannon,richness,bc_full,nrow=1,ncol=3, scale=0.9, rel_widths=c(1,1,1),labels=c("c","d","e"), align="v", axis="t")
b
three<-plot_grid(a, b, nrow=2, ncol=1, rel_widths=c(1, 1), align="h", axis="lr")
three

ggsave(here("QSU_Data/Figure4.pdf"), dpi=300, h=9, w=17)
ggsave(here("QSU_Data/Figure4.jpeg"), dpi=300, h=9, w=17)
