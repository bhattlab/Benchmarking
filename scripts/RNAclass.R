library(ggplot2)
library(tidyverse)
library(paletteer)
library(ggsignif)
library(here)
library(RColorBrewer)
library(reshape2)
library(metagenomeSeq)
library(vegan)
library(cowplot)
library(ggpubr)
library(ggnewscale)
library(Maaslin2)

##### READ COUNT PLOTS #####
preprocess_reads <- read.table(here("RNA/01_processing/readcounts.tsv"), sep="\t", header = TRUE)
preprocess_reads <- mutate(preprocess_reads, Sample=gsub("_.*", "", Sample))
preprocess_reads <- preprocess_reads %>% separate(Sample, c("Donor", "Condition", "Replicate", "Lib"), remove=FALSE)
preprocess_reads <- preprocess_reads %>% mutate(SampleID=paste(Donor, Condition, Replicate))

preprocess_grouped <- preprocess_reads %>% 
  group_by(SampleID) %>% 
  summarise_at(vars(raw_reads, dedup_reads, trimmed_reads, host_removed_reads, orphan_reads, dedup_frac, trimmed_frac, host_removed_frac, orphan_frac), sum)

rrna_reads <- read.table(here("RNA/01_rrna/aggregate/readcounts.tsv"), sep="\t", header=FALSE)
names(rrna_reads) <- c("Sample", "After_RRNA")
rrna_reads <- rrna_reads %>% separate(Sample, c("Donor", "Condition", "Replicate", "Dir"), remove=FALSE) %>% 
  filter(Dir == "1") %>% 
  mutate(SampleID=paste(Donor, Condition, Replicate)) %>% select(SampleID, After_RRNA)

readcounts <- merge(preprocess_grouped, rrna_reads, by="SampleID")

readcounts.melt.count <- melt(readcounts[,c('SampleID', 'raw_reads', 
                                            'dedup_reads', 'trimmed_reads', 
                                            'host_removed_reads', 'After_RRNA')], id.vars = 'SampleID') 

readcounts <- readcounts %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
g.count <- ggplot(readcounts.melt.count, aes(x=SampleID, y=value, fill=variable)) + 
  geom_bar(stat='identity' , position='dodge') + 
  theme_bw() + 
  scale_fill_brewer(palette = 'Set2') +
  labs(title='Readcounts at each processing step', 
       y = 'Read count') + 
  guides(fill=guide_legend(title="Processing level")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous(labels=scales::comma)
g.count


nsamp <- nrow(readcounts)
plot.width <- 2 + (nsamp/2)
plot.height <- 6

pdf(here("outputs/figures/RNA_readcounts.pdf"), height = plot.height, width=plot.width)
print(g.count)
dev.off()

write.table(readcounts, here("RNA/01_rrna/readcounts.tsv"), sep="\t", row.names = FALSE)

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

ggplot(readcounts, aes(x=raw_reads, y=After_RRNA, color=Condition)) + 
  geom_point() + 
  theme_bw() + 
  scale_x_continuous(labels=scales::comma) + 
  scale_y_continuous(labels=scales::comma) +
  scale_color_manual(values=condition_palette)

ggplot(readcounts, aes(x=Condition, y=raw_reads, fill=Condition)) + 
  geom_jitter(color="grey", width=0.2) + 
  geom_boxplot(alpha = 0.8, outlier.shape=None) + 
  theme_bw() + 
  scale_y_continuous(labels=scales::comma) +
  scale_fill_manual(values=condition_palette) + 
  stat_compare_means(method="wilcox.test", comparisons=list(c("OF", "NF")))

##### CLASSIFIED READS #####
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
#Mai did and no work metadata <- metadata %>% separate(Condition, c("Preservation", "Temperature"), remove=FALSE)
metadata <- mutate(metadata, Preservation=substr(Condition,1,1))
metadata <- mutate(metadata, Temperature=substr(Condition,2,2))

#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
metadata <- mutate(metadata, SampleCode=SampleID)

genus_kraken <- read.csv(here("RNA/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
genus_kraken <- t(genus_kraken %>% filter(row.names(genus_kraken) %in% c("unclassified")))
genus_kraken <- data.frame(genus_kraken)
names <- rownames(genus_kraken)
rownames(genus_kraken) <- NULL
genus_kraken <- cbind(names, genus_kraken)
names(genus_kraken) <- c("Sample", "Unclassified")
genus_kraken <- mutate(genus_kraken, Sample=gsub("\\.", "_", Sample))

# Merge in the metadata
colnames(metadata)[3]<-"Sample"
genus_kraken <- merge(genus_kraken, metadata, by = "Sample")

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

ggplot(genus_kraken, aes(x=Condition, y=Unclassified, fill = Condition)) + 
  geom_jitter(width = 0.2) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + 
  scale_fill_manual(values=condition_palette) + 
  ylim(0,100) +
  ylab("Unclassified Reads at Genus Level (%)") +
  theme_bw()

ggsave(here("outputs/figures/RNA_percentunclassified.jpg"), dpi=300, w = 10, h = 5.3)
ggsave(here("outputs/figures/RNA_percentunclassified.pdf"), dpi=300, w = 10, h = 5.3)

# read in DNA for comparison
genus_kraken_dna <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
genus_kraken_dna <- t(genus_kraken_dna %>% filter(row.names(genus_kraken_dna) %in% c("unclassified")))
genus_kraken_dna <- data.frame(genus_kraken_dna)
names <- rownames(genus_kraken_dna)
rownames(genus_kraken_dna) <- NULL
genus_kraken_dna <- cbind(names, genus_kraken_dna)
names(genus_kraken_dna) <- c("Sample", "Unclassified_DNA")
genus_kraken_dna <- mutate(genus_kraken_dna, Sample=gsub("\\.", "_", Sample))

genus_kraken_all <- merge(genus_kraken, genus_kraken_dna, by= "Sample")

ggplot(genus_kraken_all, aes(x=Unclassified_DNA, y=Unclassified)) + 
  geom_abline(intercept = 0, slope = 1, color = "grey") + 
  geom_point(aes(color=Condition)) + 
  theme_bw() + 
  scale_color_manual(values = condition_palette) + 
  xlab("Unclassified Reads (DNA %)") + 
  ylab("Unclassified Reads (RNA %)") + 
  ylim(0,60) + 
  xlim(0,60)

ggsave(here("outputs/figures/DNA_RNA_percentunclassified.jpg"), dpi=300, w = 8, h = 6)
ggsave(here("outputs/figures/DNA_RNA_percentunclassified.pdf"), dpi=300, w = 8, h = 6)


# compare to readcounts
preprocess_reads <- read.table(here("RNA/01_processing/readcounts.tsv"), sep="\t", header = TRUE)
preprocess_reads <- mutate(preprocess_reads, Sample=gsub("_.*", "", Sample))
preprocess_reads <- preprocess_reads %>% separate(Sample, c("Donor", "Condition", "Replicate", "Lib"), remove=FALSE)
preprocess_reads <- preprocess_reads %>% mutate(SampleID=paste(Donor, Condition, Replicate))

preprocess_grouped <- preprocess_reads %>% 
  group_by(SampleID) %>% 
  summarise_at(vars(raw_reads, dedup_reads, trimmed_reads, host_removed_reads, orphan_reads, dedup_frac, trimmed_frac, host_removed_frac, orphan_frac), sum)

rrna_reads <- read.table(here("RNA/01_rrna/aggregate/readcounts.tsv"), sep="\t", header=FALSE)
names(rrna_reads) <- c("Sample", "After_RRNA")
rrna_reads <- rrna_reads %>% separate(Sample, c("Donor", "Condition", "Replicate", "Dir"), remove=FALSE) %>% 
  filter(Dir == "1") %>% 
  mutate(SampleID=paste(Donor, Condition, Replicate)) %>% select(SampleID, After_RRNA)

readcounts <- merge(preprocess_grouped, rrna_reads, by="SampleID")

readcounts.melt.count <- melt(readcounts[,c('SampleID', 'raw_reads', 
                                            'dedup_reads', 'trimmed_reads', 
                                            'host_removed_reads', 'After_RRNA')], id.vars = 'SampleID') 

readcounts <- readcounts %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
readcounts <- mutate(readcounts, SampleID=gsub(" ", "_", SampleID)) %>% select(SampleID, After_RRNA)
names(readcounts) <- c("Sample", "After_RRNA")


count_by_classification <- merge(readcounts, genus_kraken_all, by="Sample")


ggplot(count_by_classification, aes(x=After_RRNA, y=Unclassified)) + 
  geom_point(aes(color=Condition)) + 
  scale_color_manual(values=condition_palette) + 
  theme_bw() + 
  xlab("Number of Reads") + 
  ylab("Unclassified RNA reads (%)") + 
  scale_x_continuous(labels=scales::comma)

ggsave(here("outputs/figures/DNA_RNA_reads_vs_classification.jpg"), dpi=300, w = 8, h = 6)
ggsave(here("outputs/figures/DNA_RNA_reads_vs_classification.pdf"), dpi=300, w = 8, h = 6)
##### GENUS LEVEL TAXONOMIC RELATIVE ABUNDANCE STACKED BAR PLOT, FACET BY DONOR #####

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
#Mai did and no work metadata <- metadata %>% separate(Condition, c("Preservation", "Temperature"), remove=FALSE)
metadata <- mutate(metadata, Preservation=substr(Condition,1,1))
metadata <- mutate(metadata, Temperature=substr(Condition,2,2))

#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
metadata <- mutate(metadata, SampleCode=SampleID)

genus <- read.csv(here("RNA/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)

#Color pal
n_taxa <- 20
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
bracken_long <- mutate(bracken_long, Sample=gsub("A_.*", "A", Sample))
bracken_long <- mutate(bracken_long, Sample=gsub("B_.*", "B", Sample))



# Merge in the metadata
colnames(metadata)[3]<-"Sample"
bracken_pheno <- merge(bracken_long, metadata, by = "Sample")
#bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))

# Correct the plotting order
bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)

samplabels <- bracken_pheno$Label
names(samplabels) <- bracken_pheno$Sample


condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

bracken_pheno <- bracken_pheno %>% filter(Donor != "NCO" & Donor != "PCO")
#bracken_pheno$Condition <- factor(bracken_pheno$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
bracken_pheno <- mutate(bracken_pheno, Temperature=substr(Condition, 2,2))
bracken_pheno <- mutate(bracken_pheno, PlotOrder=ifelse(Condition == "NF", 1, 
                                                        ifelse(Condition == "OF", 2, 
                                                               ifelse(Condition == "OR", 3, 
                                                                      ifelse(Condition == "OH", 4, 
                                                                             ifelse(Condition == "ZF", 5,
                                                                                    ifelse(Condition == "ZR", 6, 7)))))))

bracken_pheno <- mutate(bracken_pheno, DonorLabel=ifelse(Donor == "D10", "Donor 10", paste("Donor", substr(Donor, 3,3))))
bracken_pheno$DonorLabel <- factor(bracken_pheno$DonorLabel, levels = c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6", "Donor 7", "Donor 8", "Donor 9", "Donor 10"))

ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder), y=rel_abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  labs(
    x = "",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual("TaxaPal", values = barplot_pal) +
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
  scale_y_continuous(limits = c(-5, 100.1), expand = c(0, 0))  + 
  facet_wrap(~DonorLabel, ncol = 5, scales = "free") + 
  #theme(legend.position = "none") + 
  new_scale_fill() + 
  geom_tile(aes(x=Sample, y = -2, fill = Condition)) + 
  geom_tile(aes(x=Sample, y = -3, fill = Condition)) + 
  geom_tile(aes(x=Sample, y = -4, fill = Condition)) + 
  scale_fill_manual(values = condition_palette) +
  theme(axis.text.x = element_blank())

ggsave(here("outputs/figures/RNA_facetedgenusabundance_nicelabels.jpg"), dpi=300, w = 15, h = 5.3)
ggsave(here("outputs/figures/RNA_facetedgenusabundance_nicelabels.pdf"), dpi=300, w = 15, h = 5.3)

##### PHYLUM LEVEL TAXONOMIC RELATIVE ABUNDANCE STACKED BAR PLOT FACET BY DONOR #####

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
#Mai did and no work metadata <- metadata %>% separate(Condition, c("Preservation", "Temperature"), remove=FALSE)
metadata <- mutate(metadata, Preservation=substr(Condition,1,1))
metadata <- mutate(metadata, Temperature=substr(Condition,2,2))

#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
metadata <- mutate(metadata, SampleCode=SampleID)

genus <- read.csv(here("RNA/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE)

#Color pal
n_taxa <- 10
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
bracken_long <- mutate(bracken_long, Sample=gsub("A_.*", "A", Sample))
bracken_long <- mutate(bracken_long, Sample=gsub("B_.*", "B", Sample))



# Merge in the metadata
colnames(metadata)[3]<-"Sample"
bracken_pheno <- merge(bracken_long, metadata, by = "Sample")
#bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))

# Correct the plotting order
bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)

samplabels <- bracken_pheno$Label
names(samplabels) <- bracken_pheno$Sample


condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

bracken_pheno <- bracken_pheno %>% filter(Donor != "NCO" & Donor != "PCO")
#bracken_pheno$Condition <- factor(bracken_pheno$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
bracken_pheno <- mutate(bracken_pheno, Temperature=substr(Condition, 2,2))
bracken_pheno <- mutate(bracken_pheno, PlotOrder=ifelse(Condition == "NF", 1, 
                                                        ifelse(Condition == "OF", 2, 
                                                               ifelse(Condition == "OR", 3, 
                                                                      ifelse(Condition == "OH", 4, 
                                                                             ifelse(Condition == "ZF", 5,
                                                                                    ifelse(Condition == "ZR", 6, 7)))))))

bracken_pheno <- mutate(bracken_pheno, DonorLabel=ifelse(Donor == "D10", "Donor 10", paste("Donor", substr(Donor, 3,3))))
bracken_pheno$DonorLabel <- factor(bracken_pheno$DonorLabel, levels = c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6", "Donor 7", "Donor 8", "Donor 9", "Donor 10"))

ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder), y=rel_abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  labs(
    x = "",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual("TaxaPal", values = barplot_pal) +
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
  scale_y_continuous(limits = c(-5, 100.1), expand = c(0, 0))  + 
  facet_wrap(~DonorLabel, ncol = 5, scales = "free") + 
  #theme(legend.position = "none") + 
  new_scale_fill() + 
  geom_tile(aes(x=Sample, y = -2, fill = Condition)) + 
  geom_tile(aes(x=Sample, y = -3, fill = Condition)) + 
  geom_tile(aes(x=Sample, y = -4, fill = Condition)) + 
  scale_fill_manual(values = condition_palette) +
  theme(axis.text.x = element_blank())

ggsave(here("outputs/figures/RNA_facetedphylumabundance_nicelabels.jpg"), dpi=300, w = 15, h = 5.3)
ggsave(here("outputs/figures/RNA_facetedphylumabundance_nicelabels.pdf"), dpi=300, w = 15, h = 5.3)
##### ALPHA DIVERSITY: SHANNON  #####
benchmark_s <- read.table(here("RNA/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t")
benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")
benchmark_groups <- mutate(benchmark_groups, Preservation=substr(Condition,1,1))
benchmark_groups <- mutate(benchmark_groups, Temperature=substr(Condition,2,2))
benchmark_groups <- mutate(benchmark_groups, sample=SampleID)

shannon_div_s <- diversity(t(benchmark_s), index = "shannon")
div <- data.frame("shannon_div" = shannon_div_s, "sample" = names(shannon_div_s))

div <- mutate(div, sample=gsub("\\.", "_", sample))
div <- mutate(div, sample=gsub("A_.*", "A", sample))
div <- mutate(div, sample=gsub("B_.*", "B", sample))
div_meta <- merge(div, benchmark_groups, by="sample")

# Average pairs of libraries from same sample and replicate
div_grouped_libraries <- div_meta %>% dplyr::select(sample, shannon_div, SampleID, Donor, Condition, Replicate, Preservation, Temperature) %>% 
  group_by(SampleID) %>% 
  summarise_at(vars(shannon_div), list(shannon_div = mean))

div_grouped_libraries <- div_grouped_libraries %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
div_grouped_libraries$Condition <- factor(div_grouped_libraries$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
div_grouped_libraries$Donor <- factor(div_grouped_libraries$Donor, levels = c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
div_grouped_libraries <- mutate(div_grouped_libraries, Preservation=substr(Condition,1,1))
div_grouped_libraries <- mutate(div_grouped_libraries, Temperature=substr(Condition,2,2))


div_grouped_replicates <- div_grouped_libraries %>% dplyr::select(SampleID, shannon_div, Donor, Condition, Replicate, Preservation, Temperature) %>%
  group_by(Donor, Condition) %>%
  summarise_at(vars(shannon_div), list(shannon_div = mean))

div_grouped_replicates <- mutate(div_grouped_replicates, Preservation=substr(Condition,1,1))
div_grouped_replicates <- mutate(div_grouped_replicates, Temperature=substr(Condition,2,2))

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

div_meta <- div_meta %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "ZF", "ZR", "ZH"))
div_grouped_replicates <- div_grouped_replicates %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))

preserve_effect <- ggplot(div_grouped_replicates %>% filter(Condition %in% c("NF", "OF", "ZF")), aes(Condition, shannon_div)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85, aes(fill=Condition), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
  labs(x = "",
       y = "Shannon Diversity", 
       fill = "", 
       title = "Preservative Effect") + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF", "OF"), c("OF", "ZF"), c("NF", "ZF")), 
                     tip.length = 0, label = "p.signif")
preserve_effect

# has issues with paired test because some donors are not present in both OF and OR conditions - D10 and D05 are missing. Removed them for this calculation
temperature_effect <- ggplot(div_grouped_replicates %>% filter(Donor != "D05" && Condition!= "OF", Donor != "D10" && Condition != "OF"), aes(Condition, shannon_div)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85, aes(fill=Condition), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
  labs(x = "",
       y = "Shannon Diversity", 
       fill = "", 
       title = "Temperature Effect") + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("OF", "OR"), c("ZF", "ZR"), c("ZF", "ZH")),
                     tip.length = 0, label="p.signif")
temperature_effect

plot_grid(preserve_effect, temperature_effect, rel_widths = c(1,1.9))

ggsave(here("outputs/figures/RNA_alphadiv_paired_by_condition.jpg"), dpi=300, w=8, h=4)
ggsave(here("outputs/figures/RNA_alphadiv_paired_by_condition.pdf"), dpi=300, w=8, h=4)

##### ALPHA DIVERSITY: RICHNESS #####
benchmark_s <- read.table(here("RNA/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t")
benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")
benchmark_groups <- mutate(benchmark_groups, Preservation=substr(Condition,1,1))
benchmark_groups <- mutate(benchmark_groups, Temperature=substr(Condition,2,2))
benchmark_groups <- mutate(benchmark_groups, sample=SampleID)

# filter out all genera <0.01 abundance
benchmark_s <- mutate_all(benchmark_s, funs(ifelse(. < 0.01, 0, .)))

specnumber <- specnumber(t(benchmark_s))
div <- data.frame("SpeciesNumber" = specnumber, "sample" = names(specnumber))
div <- mutate(div, sample=gsub("\\.", "_", sample))
div <- mutate(div, sample=gsub("A_.*", "A", sample))
div <- mutate(div, sample=gsub("B_.*", "B", sample))
div_meta <- merge(div, benchmark_groups, by = "sample")

# Average pairs of libraries from same sample and replicate
div_grouped_libraries <- div_meta %>% dplyr::select(sample, SpeciesNumber, SampleID, Donor, Condition, Replicate, Preservation, Temperature) %>% 
  group_by(SampleID) %>% 
  summarise_at(vars(SpeciesNumber), list(SpeciesNumber = mean))

div_grouped_libraries <- div_grouped_libraries %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
div_grouped_libraries$Condition <- factor(div_grouped_libraries$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
div_grouped_libraries$Donor <- factor(div_grouped_libraries$Donor, levels = c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
div_grouped_libraries <- mutate(div_grouped_libraries, Preservation=substr(Condition,1,1))
div_grouped_libraries <- mutate(div_grouped_libraries, Temperature=substr(Condition,2,2))


div_grouped_replicates <- div_grouped_libraries %>% dplyr::select(SampleID, SpeciesNumber, Donor, Condition, Replicate, Preservation, Temperature) %>%
  group_by(Donor, Condition) %>%
  summarise_at(vars(SpeciesNumber), list(SpeciesNumber = mean))

div_grouped_replicates <- mutate(div_grouped_replicates, Preservation=substr(Condition,1,1))
div_grouped_replicates <- mutate(div_grouped_replicates, Temperature=substr(Condition,2,2))

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

div_meta <- div_meta %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "ZF", "ZR", "ZH"))
div_grouped_replicates <- div_grouped_replicates %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))

preserve_effect <- ggplot(div_grouped_replicates %>% filter(Condition %in% c("NF", "OF", "ZF")), aes(Condition, SpeciesNumber)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85, aes(fill=Condition), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
  labs(x = "",
       y = "Genus Number", 
       fill = "", 
       title = "Preservative Effect") + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF", "OF"), c("OF", "ZF"), c("NF", "ZF")), 
                     tip.length = 0, label = "p.signif")
preserve_effect

temperature_effect <- ggplot(div_grouped_replicates %>% filter(Donor != "D05" && Condition!= "OF", Donor != "D10" && Condition != "OF"), aes(Condition, SpeciesNumber)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85, aes(fill=Condition), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
  labs(x = "",
       y = "Genus Number", 
       fill = "", 
       title = "Temperature Effect") + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("OF", "OR"),  c("ZF", "ZR"), c("ZF", "ZH")), # list(c("OF", "OR"), c("ZF", "ZR"), c("ZF", "ZH"))
                     tip.length = 0, label="p.signif")
temperature_effect

plot_grid(preserve_effect, temperature_effect, rel_widths = c(1,1.9))

ggsave(here("outputs/figures/RNA_richness_0.01abund.jpg"), dpi=300, w=8, h=4)
ggsave(here("outputs/figures/RNA_richness_0.01abund.pdf"), dpi=300, w=8, h=4)
##### BETA DIVERSITY: BRAY CURTIS PCoA #####
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
#Mai did and no work metadata <- metadata %>% separate(Condition, c("Preservation", "Temperature"), remove=FALSE)
metadata <- mutate(metadata, Preservation=substr(Condition,1,1))
metadata <- mutate(metadata, Temperature=substr(Condition,2,2))

#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
metadata <- mutate(metadata, SampleCode=SampleID)

# Read in count data
cts <- read.csv(here("RNA/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_genus_reads.txt"), sep="\t", header = TRUE) 

cts <- cts %>% select(!starts_with("NCO") & !starts_with("PCO"))

mr <- newMRexperiment(cts)
p <- cumNormStatFast(mr)
mr_css <- cumNorm(mr, p = p)
counts_css <- MRcounts(mr_css, norm = T, log = T)
vare_dis <- vegdist(t(counts_css), method = "bray")

# calculate MDS
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)
mds_values <- mds$points
wa_scores <- wascores(mds_values, t(cts)) # check that cts is the right thing to pass in here
wa_scores <- data.frame(sample = rownames(wa_scores),
                        x = wa_scores[,1],
                        y = wa_scores[,2])
n_taxa <- 10
wa_scores_1<- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
wa_scores_2<- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
wa_scores_final <- rbind(wa_scores_1, wa_scores_2)

# calculate percentage of variation that each MDS axis accounts for
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
mds_var_per

mds_var_per_df <- data_frame(mds_var_per)
names(mds_var_per_df) <- c("PercentVariation")
mds_var_per_df$Axis <- seq.int(nrow(mds_var_per_df))
ggplot(mds_var_per_df, aes(x = Axis, y = PercentVariation)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  labs(
    x = "Dimension",
    y = "Percent Variation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "plain", size = 14),
    axis.text.x = element_text(size = 10, vjust = 0.5),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))

ggsave(here("outputs/figures/RNA_MDS_Scree_genus.pdf"), dpi=300, w=7, h=5)

#Plot
mds_data <- data.frame(Sample = gsub("\\.", "_", rownames(mds_values)),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
metadata <- filter(metadata, Donor!="NCO" & Donor!="PCO")
mds_data <-rename (mds_data, SampleID=Sample)
mds_data <- mutate(mds_data, SampleID=gsub("\\.", "_", SampleID))
mds_data <- mutate(mds_data, SampleID=gsub("A_.*", "A", SampleID))
mds_data <- mutate(mds_data, SampleID=gsub("B_.*", "B", SampleID))

mds_meta <- merge(mds_data, metadata, by = "SampleID")

maicolors <- paletteer_d("ggthemes::Tableau_10")
mds_plot <- ggplot(mds_meta, aes(x, y)) +
  stat_ellipse(aes(color = Donor), type = 't', size = 0.5, show.legend = F) +
  geom_point(size = 1, alpha = 0.7, aes(color=Donor)) +
  scale_color_manual(values = maicolors) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       title = "Genus-Level Bray Curtis",
       color = "") +
  theme_classic() +
  coord_fixed() +
  background_grid()+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

mds_plot

ggsave(here("outputs/figures/RNAGenus_Bray_Curtis.pdf"), dpi=300, w=5, h=4)
ggsave(here("outputs/figures/RNAGenus_Bray_Curtis.jpeg"), dpi=300, w=5, h=4)

#MDS colored by preservation method -NOT DONE
condition_colors <- paletteer_d("nbapalettes::lakers_alt")
mds_plot <- ggplot(mds_meta, aes(x, y)) +
  geom_point(size = 1, alpha = 0.7, aes(color=Preservation)) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       title = "Species-Level Bray Curtis",
       color = "") +
  theme_classic() +
  coord_fixed() +
  background_grid()+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

mds_plot

ggsave(here("outputs/figures/RNASpecies_Bray_Curtis_Preservation.pdf"), dpi=300, w=5, h=6)

##### BETA DIVERISTY: PHYLUM-LEVEL BRAY CURTIS BOXPLOT #####
# Read in count data
dist_matrix <- read.csv(here("RNA/02_kraken2_classification/processed_results/braycurtis_matrices/braycurtis_distance_phylum.txt"), sep="\t", header = TRUE)

# convert row names to a main column
dist_matrix <- data.frame(names = row.names(dist_matrix), dist_matrix)
rownames(dist_matrix) <- NULL

# filter rows and columns to exclude controls and contaminated sample
dist_matrix <- dist_matrix %>% select(!starts_with("NCO") & !starts_with("PCO")) # filter columns
dist_matrix <- mutate(dist_matrix, names=gsub("\\.", "_", names))
dist_matrix <- mutate(dist_matrix, names=gsub("A_.*", "A", names))
dist_matrix <- mutate(dist_matrix, names=gsub("B_.*", "B", names))

dist_matrix <- dist_matrix %>% filter(!stringr::str_detect(names, 'NCO|PCO'))

# convert data to long format
data_long <- melt(data = dist_matrix, id.vars = "names", variable.name="comparison", value.name = "bc")
names(data_long) <- c("Sample1", "Sample2", "BrayCurtis")

data_long <- mutate(data_long, Sample2=gsub("\\.", "-", Sample2))
data_long <- mutate(data_long, Sample2=gsub("A_.*", "A", Sample2))
data_long <- mutate(data_long, Sample2=gsub("B_.*", "B", Sample2))
data_long <- data_long %>% separate(Sample1, c("Donor1", "Condition1", "Replicate1", "Library1"), remove=FALSE)
data_long <- data_long %>% separate(Sample2, c("Donor2", "Condition2", "Replicate2", "Library2"), remove=FALSE)

# filter out comparisons from within replicate, two different libraries
data_long <- mutate(data_long, Sample1_nonrep_code=paste(Donor1, Condition1, Replicate1))
data_long <- mutate(data_long, Sample2_nonrep_code=paste(Donor2, Condition2, Replicate2))
data_long <- data_long %>% filter(Sample1_nonrep_code != Sample2_nonrep_code)

# create a unique identifier for the sample comparisons and filter out redundant comparisons (since this is an all by all, the table is symmetric and has duplicates)
data_long <- mutate(data_long, comparisonID = ifelse(as.character(Sample1) > as.character(Sample2), paste(Sample1, Sample2), paste(Sample2, Sample1)))
data_long_filtered <- distinct(data_long, comparisonID, .keep_all = TRUE)

# filter out self-comparisons
data_long_filtered <- data_long_filtered %>% filter(BrayCurtis > 0)

# taking means across all comparisons done with each tech rep
data_long_filtered <- mutate(data_long_filtered, nonrepID=gsub("\\.R[123]", "", comparisonID)) # create uniq id not including rep
data_long_meta <- data_long_filtered %>% select(Donor1, Condition1, Donor2, Condition2, nonrepID) %>% distinct(nonrepID, .keep_all = TRUE)
distance_averages <- data_long_filtered %>%
  group_by(nonrepID) %>%
  summarise_at(vars(BrayCurtis), list(BrayCurtis = mean))

distance_averages <- merge(distance_averages, data_long_meta)

# make dataframe that only includes within-donor comparisons, but not same condition comparisons
distance_averages_within_donor <- distance_averages %>% filter(Donor1 == Donor2) %>% filter(Condition1 != Condition2) 
distance_averages_within_donor <- mutate(distance_averages_within_donor, condition_codes = ifelse(as.character(Condition1) > as.character(Condition2), paste(Condition2, "vs", Condition1), paste(Condition1, "vs", Condition2)))

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("NF vs OH", "NF vs OR", "NF vs OF", "NF", "NF vs ZF", "NF vs ZR", "NF vs ZH")

distance_averages_within_donor <- distance_averages_within_donor %>% mutate(condition_codes = fct_relevel(condition_codes, "NF vs OF", "NF vs OR", "NF vs ZF", "NF vs ZR", "NF vs ZH"))

# FIXME: Stats don't work here due to some uneven sample sizes
temperature_effect <- ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF vs OF", "NF vs OR", "NF vs OH", "NF vs ZF", "NF vs ZR", "NF vs ZH")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85,  color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  labs(x = "",
       y = "Phylum-Level Bray Curtis Dissimilarity", 
       fill = "", 
       title = "Temperature Effect") + 
  scale_fill_manual(values = condition_palette) + 
  theme_cowplot(12) +
  theme(legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm"), 
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1)) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF vs OF", "NF vs OR"), c("NF vs ZF", "NF vs ZR"), c("NF vs ZR", "NF vs ZH"), c("NF vs ZF", "NF vs ZH")), 
                     tip.length = 0, label = "p.signif", group.by = "Donor1")

temperature_effect

# FIXME: stats don't work here due to some uneven sample sizes
kit_effect <- ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF vs OF", "NF vs ZF")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85,  color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  labs(x = "",
       y = "Phylum-Level Bray Curtis Dissimilarity", 
       fill = "", 
       title = "Preservative Effect") + 
  scale_fill_manual(values = condition_palette) + 
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF vs OF", "NF vs ZF")), 
                     tip.length = 0, label = "p.signif", group.by = "Donor1")

kit_effect

plot_grid(kit_effect, temperature_effect, rel_widths = c(1,2.35))
ggsave(here("outputs/figures/RNA_phylabraycurtis_boxplot.jpg"), dpi=300, w = 8, h = 5)
ggsave(here("outputs/figures/RNA_phylabraycurtis_boxplot.pdf"), dpi=300, w = 8, h = 5)

##### BETA DIVERSITY: GENUS-LEVEL BRAY CURTIS BOXPLOT #####
# Read in count data
dist_matrix <- read.csv(here("RNA/02_kraken2_classification/processed_results/braycurtis_matrices/braycurtis_distance_genus.txt"), sep="\t", header = TRUE)

# convert row names to a main column
dist_matrix <- data.frame(names = row.names(dist_matrix), dist_matrix)
rownames(dist_matrix) <- NULL

# filter rows and columns to exclude controls and contaminated sample
dist_matrix <- dist_matrix %>% select(!starts_with("NCO") & !starts_with("PCO")) # filter columns
dist_matrix <- mutate(dist_matrix, names=gsub("\\.", "_", names))
dist_matrix <- mutate(dist_matrix, names=gsub("A_.*", "A", names))
dist_matrix <- mutate(dist_matrix, names=gsub("B_.*", "B", names))

dist_matrix <- dist_matrix %>% filter(!stringr::str_detect(names, 'NCO|PCO'))

# convert data to long format
data_long <- melt(data = dist_matrix, id.vars = "names", variable.name="comparison", value.name = "bc")
names(data_long) <- c("Sample1", "Sample2", "BrayCurtis")

data_long <- mutate(data_long, Sample2=gsub("\\.", "-", Sample2))
data_long <- mutate(data_long, Sample2=gsub("A_.*", "A", Sample2))
data_long <- mutate(data_long, Sample2=gsub("B_.*", "B", Sample2))
data_long <- data_long %>% separate(Sample1, c("Donor1", "Condition1", "Replicate1", "Library1"), remove=FALSE)
data_long <- data_long %>% separate(Sample2, c("Donor2", "Condition2", "Replicate2", "Library2"), remove=FALSE)

# filter out comparisons from within replicate, two different libraries
data_long <- mutate(data_long, Sample1_nonrep_code=paste(Donor1, Condition1, Replicate1))
data_long <- mutate(data_long, Sample2_nonrep_code=paste(Donor2, Condition2, Replicate2))
data_long <- data_long %>% filter(Sample1_nonrep_code != Sample2_nonrep_code)

# create a unique identifier for the sample comparisons and filter out redundant comparisons (since this is an all by all, the table is symmetric and has duplicates)
data_long <- mutate(data_long, comparisonID = ifelse(as.character(Sample1) > as.character(Sample2), paste(Sample1, Sample2), paste(Sample2, Sample1)))
data_long_filtered <- distinct(data_long, comparisonID, .keep_all = TRUE)

# filter out self-comparisons
data_long_filtered <- data_long_filtered %>% filter(BrayCurtis > 0)

# taking means across all comparisons done with each tech rep
data_long_filtered <- mutate(data_long_filtered, nonrepID=gsub("\\.R[123]", "", comparisonID)) # create uniq id not including rep
data_long_meta <- data_long_filtered %>% select(Donor1, Condition1, Donor2, Condition2, nonrepID) %>% distinct(nonrepID, .keep_all = TRUE)
distance_averages <- data_long_filtered %>%
  group_by(nonrepID) %>%
  summarise_at(vars(BrayCurtis), list(BrayCurtis = mean))

distance_averages <- merge(distance_averages, data_long_meta)

# make dataframe that only includes within-donor comparisons, but not same condition comparisons
distance_averages_within_donor <- distance_averages %>% filter(Donor1 == Donor2) %>% filter(Condition1 != Condition2) 
distance_averages_within_donor <- mutate(distance_averages_within_donor, condition_codes = ifelse(as.character(Condition1) > as.character(Condition2), paste(Condition2, "vs", Condition1), paste(Condition1, "vs", Condition2)))

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("NF vs OH", "NF vs OR", "NF vs OF", "NF", "NF vs ZF", "NF vs ZR", "NF vs ZH")

distance_averages_within_donor <- distance_averages_within_donor %>% mutate(condition_codes = fct_relevel(condition_codes, "NF vs OF", "NF vs OR", "NF vs ZF", "NF vs ZR", "NF vs ZH"))

# FIXME stats don't work because of uneven sample size
temperature_effect <- ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF vs OF", "NF vs OR", "NF vs OH", "NF vs ZF", "NF vs ZR", "NF vs ZH")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85,  color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  labs(x = "",
       y = "Genus-Level Bray Curtis Dissimilarity", 
       fill = "", 
       title = "Temperature Effect") + 
  scale_fill_manual(values = condition_palette) + 
  theme_cowplot(12) +
  theme(legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm"), 
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1)) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF vs OF", "NF vs OR"), c("NF vs OR", "NF vs OH"), c("NF vs OF", "NF vs OH"), c("NF vs ZF", "NF vs ZR"), c("NF vs ZR", "NF vs ZH"), c("NF vs ZF", "NF vs ZH")), 
                     tip.length = 0, label = "p.signif", group.by = "Donor1")

temperature_effect

# FIXME stats don't work because of uneven sample size
kit_effect <- ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF vs OF", "NF vs ZF")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85,  color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  labs(x = "",
       y = "Genus-Level Bray Curtis Dissimilarity", 
       fill = "", 
       title = "Preservative Effect") + 
  scale_fill_manual(values = condition_palette) + 
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF vs OF", "NF vs ZF")), 
                     tip.length = 0, label = "p.signif", group.by = "Donor1")

kit_effect

plot_grid(kit_effect, temperature_effect, rel_widths = c(1,2.35))
ggsave(here("outputs/figures/RNA_braycurtis_genusboxplot.jpg"), dpi=300, w = 8, h = 5)
ggsave(here("outputs/figures/RNA_braycurtis_genusboxplot.pdf"), dpi=300, w = 8, h = 5)


##### GENUS LEVEL DIFFERENTIAL ABUNDANCE: MAASLIN2 MODEL #####
#Installation of Maaslin2
#if(!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install(version = "3.12") 
#BiocManager::install("Maaslin2", version="3.12")

#Other installation method
#install.packages("devtools")
# library(devtools)
#devtools::install_github("biobakery/maaslin2")
# to get help on usage of Maaslin2 ?

library(Maaslin2)

#Reading in metadata dataframe and creating metadata file for Maaslin2
#benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
#benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#benchmark_groups <- mutate(benchmark_groups, sample=gsub("_", ".", SampleID))
#benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")
#benchmark_groups <- mutate(benchmark_groups, Kit=substr(Condition,1,1))
#benchmark_groups <- mutate(benchmark_groups, Temperature=substr(Condition,2,2))
#benchmark_groups <- mutate(benchmark_groups, Kit=ifelse(Kit=="N", "No Preservative", ifelse(Kit=="O", "Omnigene", "Zymo")))
#benchmark_groups <- mutate(benchmark_groups, Temperature=ifelse(Temperature=="F", "Frozen", ifelse(Temperature=="R", "Room Temp", "Hot")))
#rownames(benchmark_groups) <- benchmark_groups$sample

#Exporting metdata
#write.table(benchmark_groups, here("outputs/tables/DNA_metadata.tsv"), row.names = FALSE, sep="\t", quote=FALSE)  

#Transpose the classification file
kraken_genus_og <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)

kraken_genus <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE))
#kraken_genus <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE))

#Running Maaslin2
#fit_data <- Maaslin2(
  #kraken_genus, benchmark_groups, here('DNA/3.maaslin'), transform = "NONE",
  #fixed_effects = c('Kit', 'Temperature'),
  #random_effects = c('Donor'),
  #normalization = 'NONE',
  #standardize = FALSE,
  #min_abundance = 0.01,
  #min_prevalence = 0.1, 
  #reference = 'Kit,No Preservative;Temperature,Frozen')

#Exporting Maaslin2 results
#write.csv(fit_data, here("DNA/3.maaslin/DNA_maaslin.csv"), row.names = FALSE)


##Maaslin analysis on only Zymo conditions
#Filter metdata to only include Zymo to do Zymo only Maaslin analysis
benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample
benchmark_zymo <- filter(benchmark_groups, Kit=="Zymo")


##Maaslin2 analysis of ZF v ZH
#Filter metadata to only include Zymo kit
#benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
#rownames(benchmark_groups) <- benchmark_groups$sample
#zymoresults <- filter(benchmark_groups, Kit=="Zymo")

#Running Maaslin2
fit_data <- Maaslin2(
  #kraken_genus, zymoresults, here('DNA/3.maaslin'), transform = "NONE",
  kraken_genus, benchmark_zymo, here('DNA/3.maaslin/zymo'), transform = "NONE",
  fixed_effects = c('Temperature'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Temperature,Frozen')

#Filter for significant and large effect sizes
#results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)

#Exporting Maaslin2 results
write.csv(fit_data, here("DNA/3.maaslin/zymo/DNA_maaslin_zymo.csv"), row.names = FALSE)

#Filter Maaslin results to include significant effect sizes
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#Comparison of ZR and ZF
zymoresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Frozen", "Hot"))

ggplot(zymoresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

ggsave(here("outputs/figures/DNA_ZFZH_maaslin2.pdf"), dpi=300, w=5, h=6)
ggsave(here("outputs/figures/DNA_ZFZH_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)

##Maaslin2 analysis of OF v OH
#Filter metadata to only include Omni kit
#benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
#rownames(benchmark_groups) <- benchmark_groups$sample
#omniresults <- filter(benchmark_groups, Kit=="Omnigene")

#Running Maaslin2
#fit_data <- Maaslin2(
 # kraken_genus, omniresults, here('DNA/3.maaslin'), transform = "NONE",

##Maaslin analysis on only Omni conditions
#Filter metdata to only include Omni to do Omni only Maaslin analysis
benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample
benchmark_omni <- filter(benchmark_groups, Kit=="Omnigene")

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, benchmark_omni, here('DNA/3.maaslin/omni'), transform = "NONE",
  fixed_effects = c('Temperature'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Temperature,Frozen')

#Filter for significant and large effect sizes
#results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#Exporting Maaslin2 results
write.csv(fit_data, here("DNA/3.maaslin/omni/DNA_maaslin_omni.csv"), row.names = FALSE)

#Filter Maaslin results to include significant effect sizes
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#Comparison of ZR and ZF
omniresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Frozen", "Hot"))

ggplot(omniresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf"), dpi=300, w=5, h=6)
ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)


##Maaslin2 analysis of NF v ZF v OF
#Filter metadata to only include frozen temperature
#benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
#rownames(benchmark_groups) <- benchmark_groups$sample
#frozenresults <- filter(benchmark_groups, Temperature=="Frozen")

#Running Maaslin2
#fit_data <- Maaslin2(
  #kraken_genus, frozenresults, here('DNA/3.maaslin'), transform = "NONE",
##Maaslin analysis on only frozen samples
#Filter metdata to only include Frozen
benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample
benchmark_frozen <- filter(benchmark_groups, Temperature=="Frozen")

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, benchmark_frozen, here('DNA/3.maaslin/frozen'), transform = "NONE",
  fixed_effects = c('Kit'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Kit,No Preservative')

# #Filter for significant and large effect sizes only for NF v ZF
# results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# nfzfresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Zymo"))
# 
# ggplot(nfzfresults, aes(x=reorder(feature, -coef), y=coef)) +

#Exporting Maaslin2 results
write.csv(fit_data, here("DNA/3.maaslin/frozen/DNA_maaslin_frozen.csv"), row.names = FALSE)

#Filter Maaslin results to include significant effect sizes
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#Comparison of ZR and ZF
frozenresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Zymo"))

ggplot(frozenresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

# ggsave(here("outputs/figures/DNA_NFZF_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_NFZF_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# 
# #Filter for significant and large effect sizes only for NF v OF
# results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# nfofresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Omnigene"))
# 
# ggplot(nfofresults, aes(x=reorder(feature, -coef), y=coef)) +
#   geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
#   theme_bw() + 
#   labs(y = "Effect Size", x = "Genus", fill = "") +
#   coord_flip() + 
#   #scale_fill_manual(values = flexxt_palette) + 
#   theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")
# 
# ggsave(here("outputs/figures/DNA_NFOF_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_NFOF_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# 
# 
# ##Generating heatmap for Kit and Temp comparison
# #Read in dataframes
# zymoresults <- read.table(here("DNA/3.maaslin/zymo/all_results.tsv"), sep="\t", header=TRUE)
# omniresults <- read.table(here("DNA/3.maaslin/omni/all_results.tsv"), sep="\t", header=TRUE)
# frozenresults <- read.table(here("DNA/3.maaslin/frozen/all_results.tsv"), sep="\t", header=TRUE)
# 
# zymoresults <- filter(zymoresults, value=="Hot") %>% mutate(index=paste(metadata, "Zymo"))
# omniresults <- filter(omniresults, value=="Hot") %>% mutate(index=paste(metadata, "Omni"))
# frozenresults <- mutate(frozenresults, index=paste(metadata,value))
# 
# #Concatenate the dataframes
# results <- rbind(zymoresults, frozenresults, omniresults)
# write.table(results, here("DNA/3.maaslin/concat_results.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
# 
# #Generating the heatmap
# results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
# results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl(".environmental", feature))
# sigresults <- filter(results, abs(coef)>0.25)
# uniquetaxa <- unique(sigresults$feature)
# results <- filter(results, feature %in% uniquetaxa)
# results <- mutate(results, coef2=ifelse(coef>3, 3, ifelse(coef < -3, -3, coef)))
# limit <- max(abs(results$coef2)) * c(-1, 1)
# 
# 
# 
# ggplot(results, aes(x=index, y=feature, fill=coef2)) + 
ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf"), dpi=300, w=5, h=6)
ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)

##Generating heatmap for Maaslin comparison of Kit and Temperature
#Reading in dataframes and adding in column of comparison to each Maaslin analysis
frozenresults <- read.table(here("DNA/3.maaslin/frozen/all_results.tsv"), sep="\t", header=TRUE)
omniresults <- read.table(here("DNA/3.maaslin/omni/all_results.tsv"), sep="\t", header=TRUE)
zymoresults <- read.table(here("DNA/3.maaslin/zymo/all_results.tsv"), sep="\t", header=TRUE)

frozenresults <- mutate(frozenresults, index=paste(metadata,value))
omniresults <- filter(omniresults, value=="Hot") %>% mutate(index=paste(metadata,"Omni"))
zymoresults <- filter(zymoresults, value=="Hot") %>% mutate(index=paste(metadata,"Zymo"))

#Concatenate the dataframes together
results <- rbind(frozenresults,omniresults,zymoresults)

#Making the heatmap
results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
limit <- max(abs(results$coef)) * c(-1, 1)
results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl("environmental.samples", feature))
results <- mutate(results, coef2=ifelse(coef<3, coef, 4))


ggplot(results, aes(x=index, y=feature, fill=coef)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(~metadata, scales="free_x") + 
  # scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  # geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
  # scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  # theme(panel.grid = element_blank(),
  #       panel.border = element_blank()) + 
  # theme(axis.text.y = element_text(size = 15)) + 
  #scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "white", size = 3, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 12)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 

ggsave(here("outputs/figures/DNA_heatmap.pdf"), dpi=300, w=12, h=12)
ggsave(here("outputs/figures/DNA_heatmap.jpeg"), dpi=300, w=12, h=12)

ggplot(results, aes(x=coef)) + 
  geom_histogram()+
  xlim(c(-1,1))

#Generating the heatmap with gram stain and phyla
results <- read.table(here("DNA/3.maaslin/concat_results.tsv"), sep="\t", header=TRUE)
results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl(".environmental", feature))
sigresults <- filter(results, abs(coef)>0.25)
uniquetaxa <- unique(sigresults$feature)
results <- filter(results, feature %in% uniquetaxa)
results <- mutate(results, coef2=ifelse(coef>3, 3, ifelse(coef < -3, -3, coef)))
limit <- max(abs(results$coef2)) * c(-1, 1)
#Merge with gram stain and phyla dataframe
phylagram <- read.csv(here("outputs/tables/Maaslin_heatmap_126.csv"), sep=",", header=TRUE)
phylagram <- merge(results, phylagram, by='feature')

#Plotting phyla heatmap
ggplot(phylagram, aes(x=index, y=feature, fill=coef2)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(Phylum~metadata, scales="free", space="free_y") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 13)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 

ggsave(here("outputs/figures/DNA_heatmap_phyla.pdf"), dpi=300, w=12, h=14)
ggsave(here("outputs/figures/DNA_heatmap_phyla.jpeg"), dpi=300, w=12, h=14)

#Plotting gram stain heatmap
ggplot(phylagram, aes(x=index, y=feature, fill=coef2)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(GramStatus~metadata, scales="free", space="free_y") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 13)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 

ggsave(here("outputs/figures/DNA_heatmap_gram.pdf"), dpi=300, w=12, h=14)
ggsave(here("outputs/figures/DNA_heatmap_gram.jpeg"), dpi=300, w=12, h=12)


#Grouped Maaslin analysis genus-level
#Read in genus-level Kraken data and metadata
kraken_genus <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE))
benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, benchmark_groups, here('DNA/3.maaslin/genuscomplete'), transform = "NONE",
  fixed_effects = c('Condition'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = c(('Condition,NF')))

#Reading in Maaslin results
genusresults <- read.table(here("DNA/3.maaslin/genuscomplete/all_results.tsv"), sep="\t", header=TRUE)

#Filter Maaslin results to include significant effect sizes
genusresults <- mutate(genusresults, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
#limit <- max(abs(genusresults$coef)) * c(-1, 1)
genusresults <- filter(genusresults, !grepl("unclassified", feature)) %>% filter(!grepl("environmental.samples", feature))
genusresults <- mutate(genusresults, coef2=ifelse(coef<3, coef, 4))
genusresults <- mutate(genusresults, index=paste(metadata,value))
#Filter out low effect size differences
genusresults <- genusresults %>% filter(abs(coef2)>0.5)
limit <- max(abs(genusresults$coef2))* c(-1,1)

ggplot(genusresults, aes(x=index, y=feature, fill=coef2)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(~metadata, scales="free", space="free_y") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "dark grey", size = 5, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 13)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 

ggsave(here("outputs/figures/DNA_groupedheatmap_genus_coef.pdf"), dpi=300, w=10, h=12)
ggsave(here("outputs/figures/DNA_groupedheatmap_genus_coef.jpeg"), dpi=300, w=10, h=12)

###Genus-level heatmap faceted by phyla or gram status
#Read in data on phyla and gram status
phylagram <- read.csv(here("outputs/tables/Maaslin_heatmap_32.csv"), sep=",", header=TRUE)
phylagram <- merge(genusresults, phylagram, by='feature')

#Plotting gram stain heatmap
ggplot(phylagram, aes(x=index, y=feature, fill=coef2)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(GramStatus~metadata, scales="free", space="free_y") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "dark grey", size = 5, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 13)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 
ggsave(here("outputs/figures/DNA_groupedheatmap_genus_coefgram.pdf"), dpi=300, w=10, h=12)
ggsave(here("outputs/figures/DNA_groupedheatmap_genus_coefgram.jpeg"), dpi=300, w=10, h=12)

#Plotting by phylum heatmap
ggplot(phylagram, aes(x=index, y=feature, fill=coef2)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(Phylum~metadata, scales="free", space="free_y") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "dark grey", size = 5, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 13)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 
ggsave(here("outputs/figures/DNA_groupedheatmap_genus_coefphyla.pdf"), dpi=300, w=10, h=12)
ggsave(here("outputs/figures/DNA_groupedheatmap_genus_coefphyla.jpeg"), dpi=300, w=10, h=12)

##### PHYLUM LEVEL TAXONOMIC RELATIVE ABUNDANCE STACKED BAR PLOT, FACET BY DONOR#####
# read in metadata
# metadata <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
# 
# #Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
# metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
# #Modify Condition column, so that anything labeled with B# is changed to Controls
# metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
# #Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
# metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
# #Separate the Condition column to create preservation method and temperature columns
# #Mai did and no work metadata <- metadata %>% separate(Condition, c("Preservation", "Temperature"), remove=FALSE)
# metadata <- mutate(metadata, Preservation=substr(Condition,1,1))
# metadata <- mutate(metadata, Temperature=substr(Condition,2,2))
# 
# #Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
# metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
# #Modify Condition column, so that anything labeled with B# is changed to Controls
# metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
# #Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
# metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
# metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
# metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# 
# phyla <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE)
# 
# #Color pal
# n_taxa <- 20
# myCols <- colorRampPalette(brewer.pal(9, "Set1")) # WAS Set1
# barplot_pal <- myCols(n_taxa)
# barplot_pal <- sample(barplot_pal)
# barplot_pal[n_taxa + 1] <- "gray"
# 
# abundance_threshold <- sort(rowSums(phyla), decreasing = T)[n_taxa]
# bracken_plot <- phyla[rowSums(phyla) >= abundance_threshold,]
# bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  100 - colSums(bracken_plot))))
# 
# bracken_plot$Phyla <- row.names(bracken_plot)
# bracken_plot$Phyla <- gsub("\\(miscellaneous\\)", "", bracken_plot$Phyla)
# bracken_long <- melt(bracken_plot, id.vars = "Phyla", variable.name = "Sample", value.name = "rel_abundance")
# bracken_long <- mutate(bracken_long, Sample=gsub("\\.", "_", Sample))
# 
# # Merge in the metadata
# colnames(metadata)[3]<-"Sample"
# 
# bracken_pheno <- merge(bracken_long, metadata, by = "Sample")
# #bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))
# 
# # Correct the plotting order
# bracken_pheno$Phyla <- factor(bracken_pheno$Phyla, levels = bracken_plot$Phyla)
# 
# samplabels <- bracken_pheno$Label
# names(samplabels) <- bracken_pheno$Sample
# 
# condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
# names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")
# 
# bracken_pheno <- bracken_pheno %>% filter(Donor != "NCO" & Donor != "PCO")
# #bracken_pheno$Condition <- factor(bracken_pheno$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# bracken_pheno <- mutate(bracken_pheno, Temperature=substr(Condition, 2,2))
# bracken_pheno <- mutate(bracken_pheno, PlotOrder=ifelse(Condition == "NF", 1, 
#                                                         ifelse(Condition == "OF", 2, 
#                                                                ifelse(Condition == "OR", 3, 
#                                                                       ifelse(Condition == "OH", 4, 
#                                                                              ifelse(Condition == "ZF", 5,
#                                                                                     ifelse(Condition == "ZR", 6, 7)))))))
# 
# bracken_pheno <- mutate(bracken_pheno, DonorLabel=ifelse(Donor == "D10", "Donor 10", paste("Donor", substr(Donor, 3,3))))
# bracken_pheno$DonorLabel <- factor(bracken_pheno$DonorLabel, levels = c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6", "Donor 7", "Donor 8", "Donor 9", "Donor 10"))
# 
# ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder), y=rel_abundance, fill=Phyla)) +
#   geom_bar(stat="identity") +
#   labs(
#     x = "",
#     y = "Relative Abundance (%)"
#   ) +
#   scale_fill_manual("Phylum", values = barplot_pal) +
#   guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
#   theme_bw() +
#   scale_x_discrete(labels = samplabels) + 
#   theme(
#     plot.title = element_text(face = "plain", size = 14),
#     legend.text = element_text(size = 10),
#     #axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     panel.grid = element_blank(), 
#     panel.border = element_blank(),
#     strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
#     strip.text = element_text(color = "red", size = 12)) +
#   scale_y_continuous(limits = c(-5, 100.1), expand = c(0, 0))  + 
#   facet_wrap(~DonorLabel, ncol = 5, scales = "free") + 
#   #theme(legend.position = "none") + 
#   new_scale_fill() + 
#   geom_tile(aes(x=Sample, y = -2, fill = Condition)) + 
#   geom_tile(aes(x=Sample, y = -3, fill = Condition)) + 
#   geom_tile(aes(x=Sample, y = -4, fill = Condition)) + 
#   scale_fill_manual(values = condition_palette)

#ggsave(here("outputs/figures/DNA_facetedphylaabundance_nicelabels.jpg"), dpi=300, w = 15, h = 6)
#ggsave(here("outputs/figures/DNA_facetedphylaabundance_nicelabels.pdf"), dpi=300, w = 15, h = 6)

##### PHYLUM ALPHA DIVERSITY: SHANNON #####
# benchmark_s <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t")
# benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
# #Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
# benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
# benchmark_groups <- mutate(benchmark_groups, sample=gsub("_", ".", SampleID))
# benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")
# benchmark_groups <- mutate(benchmark_groups, Preservation=substr(Condition,1,1))
# benchmark_groups <- mutate(benchmark_groups, Temperature=substr(Condition,2,2))
# 
# shannon_div_s <- diversity(t(benchmark_s), index = "shannon")
# div <- data.frame("shannon_div" = shannon_div_s, "sample" = names(shannon_div_s))
# div_meta <- merge(div, benchmark_groups, by = "sample")
# div_meta <- filter(div_meta, sample != "D03.ZH.R2")
# 
# names(div_meta)
# div_grouped_replicates <- div_meta %>% dplyr::select(sample, shannon_div, SampleID, Donor, Condition, Replicate, Preservation, Temperature) %>%
#   group_by(Donor, Condition) %>%
#   summarise_at(vars(shannon_div), list(shannon_div = mean))
# 
# div_grouped_replicates <- mutate(div_grouped_replicates, Preservation=substr(Condition,1,1))
# div_grouped_replicates <- mutate(div_grouped_replicates, Temperature=substr(Condition,2,2))
# 
# condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
# names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")
# 
# div_meta <- div_meta %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# div_grouped_replicates <- div_grouped_replicates %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
# 
# preserve_effect <- ggplot(div_grouped_replicates %>% filter(Condition %in% c("NF", "OF", "ZF")), aes(Condition, shannon_div)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85, aes(fill=Condition), color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
#   labs(x = "",
#        y = "Shannon Diversity", 
#        fill = "", 
#        title = "Preservative Effect") + 
#   ylim(0.7,1.25) +
#   scale_color_manual(values = condition_palette) + 
#   scale_fill_manual(values = condition_palette) +
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(size = 12),
#         legend.position = "none",
#         title = element_text(size = 13),
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
#   stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF", "OF"), c("OF", "ZF"), c("NF", "ZF")), 
#                      tip.length = 0, label = "p.signif")
# preserve_effect
# 
# temperature_effect <- ggplot(div_grouped_replicates , aes(Condition, shannon_div)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85, aes(fill=Condition), color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
#   labs(x = "",
#        y = "Shannon Diversity", 
#        fill = "", 
#        title = "Temperature Effect") + 
#   ylim(0.7,1.25) +
#   scale_color_manual(values = condition_palette) + 
#   scale_fill_manual(values = condition_palette) +
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(size = 12),
#         legend.position = "none",
#         title = element_text(size = 13),
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
#   stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("OF", "OR"), c("OF", "OH"), c("ZF", "ZR"), c("ZF", "ZH")), 
#                      tip.length = 0, label="p.signif")
# temperature_effect
# 
# plot_grid(preserve_effect, temperature_effect, rel_widths = c(1,1.9))

#ggsave(here("outputs/figures/DNA_alphadiv_phyla.jpg"), dpi=300, w = 7, h = 4)
#ggsave(here("outputs/figures/DNA_alphadiv_phyla.pdf"), dpi=300, w = 7, h = 4)




##### BETA DIVERSITY:PHYLYM-LEVEL BRAY CURTIS BOXPLOT#####
# Read in count data
# dist_matrix <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/braycurtis_matrices/braycurtis_distance_species.txt"), sep="\t", header = TRUE)
# 
# # convert row names to a main column
# dist_matrix <- data.frame(names = row.names(dist_matrix), dist_matrix)
# rownames(dist_matrix) <- NULL
# 
# # filter rows and columns to exclude controls and contaminated sample
# dist_matrix <- dist_matrix %>% select(!starts_with("NCO") & !starts_with("PCO") & !matches("D03.ZH.R2")) # filter columns
# dist_matrix <- mutate(dist_matrix, names = gsub('-', '.', names))
# dist_matrix <- dist_matrix %>% filter(!stringr::str_detect(names, 'NCO|PCO') & !stringr::str_detect(names, 'D03.ZH.R2'))
# 
# # convert data to long format
# data_long <- melt(data = dist_matrix, id.vars = "names", variable.name="comparison", value.name = "bc")
# names(data_long) <- c("Sample1", "Sample2", "BrayCurtis")
# 
# data_long <- data_long %>% separate(Sample1, c("Donor1", "Condition1", "Replicate1"), remove=FALSE)
# data_long <- data_long %>% separate(Sample2, c("Donor2", "Condition2", "Replicate2"), remove=FALSE)
# 
# # create a unique identifier for the sample comparisons and filter out redundant comparisons (since this is an all by all, the table is symmetric and has duplicates)
# data_long <- mutate(data_long, comparisonID = ifelse(as.character(Sample1) > as.character(Sample2), paste(Sample1, Sample2), paste(Sample2, Sample1)))
# data_long_filtered <- distinct(data_long, comparisonID, .keep_all = TRUE)
# 
# # filter out self-comparisons
# data_long_filtered <- data_long_filtered %>% filter(BrayCurtis > 0)
# 
# # taking means across all comparisons done with each tech rep
# data_long_filtered <- mutate(data_long_filtered, nonrepID=gsub("\\.R[123]", "", comparisonID)) # create uniq id not including rep
# data_long_meta <- data_long_filtered %>% select(Donor1, Condition1, Donor2, Condition2, nonrepID) %>% distinct(nonrepID, .keep_all = TRUE)
# distance_averages <- data_long_filtered %>%
#   group_by(nonrepID) %>%
#   summarise_at(vars(BrayCurtis), list(BrayCurtis = mean))
# 
# distance_averages <- merge(distance_averages, data_long_meta)
# 
# # make dataframe that only includes within-donor comparisons, but not same condition comparisons
# distance_averages_within_donor <- distance_averages %>% filter(Donor1 == Donor2) %>% filter(Condition1 != Condition2) 
# distance_averages_within_donor <- mutate(distance_averages_within_donor, condition_codes = ifelse(as.character(Condition1) > as.character(Condition2), paste(Condition2, "vs", Condition1), paste(Condition1, "vs", Condition2)))
# 
# condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
# names(condition_palette) <- c("NF vs OH", "NF vs OR", "NF vs OF", "NF", "NF vs ZF", "NF vs ZR", "NF vs ZH")
# 
# distance_averages_within_donor <- distance_averages_within_donor %>% mutate(condition_codes = fct_relevel(condition_codes, "NF vs OF", "NF vs OR", "NF vs OH", "NF vs ZF", "NF vs ZR", "NF vs ZH"))
# 
# 
# temperature_effect <- ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF vs OF", "NF vs OR", "NF vs OH", "NF vs ZF", "NF vs ZR", "NF vs ZH")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85,  color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#   labs(x = "",
#        y = "Species-Level Bray Curtis Dissimilarity", 
#        fill = "", 
#        title = "Temperature Effect") + 
#   ylim(0, 0.5) + 
#   scale_fill_manual(values = condition_palette) + 
#   theme_cowplot(12) +
#   theme(legend.position = "none",
#         title = element_text(size = 13),
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm"), 
#         axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1)) +
#   stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF vs OF", "NF vs OR"), c("NF vs OR", "NF vs OH"), c("NF vs OF", "NF vs OH"), c("NF vs ZF", "NF vs ZR"), c("NF vs ZR", "NF vs ZH"), c("NF vs ZF", "NF vs ZH")), 
#                      tip.length = 0, label = "p.signif", group.by = "Donor1")
# 
# temperature_effect
# 
# 
# kit_effect <- ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF vs OF", "NF vs ZF")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85,  color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#   labs(x = "",
#        y = "Species-Level Bray Curtis Dissimilarity", 
#        fill = "", 
#        title = "Preservative Effect") + 
#   ylim(0, 0.5) + 
#   scale_fill_manual(values = condition_palette) + 
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
#         legend.position = "none",
#         title = element_text(size = 13),
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
#   stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF vs OF", "NF vs ZF")), 
#                      tip.length = 0, label = "p.signif", group.by = "Donor1")
# 
# kit_effect
# 
# plot_grid(kit_effect, temperature_effect, rel_widths = c(1,2.35))
# #ggsave(here("outputs/figures/DNA_braycurtis_boxplot.jpg"), dpi=300, w = 8, h = 5)
# #ggsave(here("outputs/figures/DNA_braycurtis_boxplot.pdf"), dpi=300, w = 8, h = 5)
# ##### PHYLUM DIFFERENTIAL ABUNDANCE: MAASLIN2
# 
# library(Maaslin2)
# 
# #Transpose the classification file
# kraken_phyla_og <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE)
# 
# kraken_phyla <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE))
# 
# 
# ##Maaslin analysis on only Zymo conditions
# #Filter metdata to only include Zymo to do Zymo only Maaslin analysis
# benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# rownames(benchmark_groups) <- benchmark_groups$sample
# benchmark_zymo <- filter(benchmark_groups, Kit=="Zymo")
# 
# 
# ##Maaslin2 analysis of ZF v ZH
# #Filter metadata to only include Zymo kit
# #benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# #rownames(benchmark_groups) <- benchmark_groups$sample
# #zymoresults <- filter(benchmark_groups, Kit=="Zymo")
# 
# #Running Maaslin2
# fit_data <- Maaslin2(
#   #kraken_genus, zymoresults, here('DNA/3.maaslin'), transform = "NONE",
#   kraken_genus, benchmark_zymo, here('DNA/3.maaslin/zymo'), transform = "NONE",
#   fixed_effects = c('Temperature'),
#   random_effects = c('Donor'),
#   normalization = 'NONE',
#   standardize = FALSE,
#   min_abundance = 0.01,
#   min_prevalence = 0.1, 
#   reference = 'Temperature,Frozen')
# 
# #Filter for significant and large effect sizes
# #results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# 
# #Exporting Maaslin2 results
# write.csv(fit_data, here("DNA/3.maaslin/zymo/DNA_maaslin_zymo.csv"), row.names = FALSE)
# 
# #Filter Maaslin results to include significant effect sizes
# results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# #Comparison of ZR and ZF
# zymoresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Frozen", "Hot"))
# 
# ggplot(zymoresults, aes(x=reorder(feature, -coef), y=coef)) +
#   geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
#   theme_bw() + 
#   labs(y = "Effect Size", x = "Genus", fill = "") +
#   coord_flip() + 
#   #scale_fill_manual(values = flexxt_palette) + 
#   theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")
# 
# ggsave(here("outputs/figures/DNA_ZFZH_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_ZFZH_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# 
# ##Maaslin2 analysis of OF v OH
# #Filter metadata to only include Omni kit
# #benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# #rownames(benchmark_groups) <- benchmark_groups$sample
# #omniresults <- filter(benchmark_groups, Kit=="Omnigene")
# 
# #Running Maaslin2
# #fit_data <- Maaslin2(
# # kraken_genus, omniresults, here('DNA/3.maaslin'), transform = "NONE",
# 
# ##Maaslin analysis on only Omni conditions
# #Filter metdata to only include Omni to do Omni only Maaslin analysis
# benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# rownames(benchmark_groups) <- benchmark_groups$sample
# benchmark_omni <- filter(benchmark_groups, Kit=="Omnigene")
# 
# #Running Maaslin2
# fit_data <- Maaslin2(
#   kraken_genus, benchmark_omni, here('DNA/3.maaslin/omni'), transform = "NONE",
#   fixed_effects = c('Temperature'),
#   random_effects = c('Donor'),
#   normalization = 'NONE',
#   standardize = FALSE,
#   min_abundance = 0.01,
#   min_prevalence = 0.1, 
#   reference = 'Temperature,Frozen')
# 
# #Filter for significant and large effect sizes
# #results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# #Exporting Maaslin2 results
# write.csv(fit_data, here("DNA/3.maaslin/omni/DNA_maaslin_omni.csv"), row.names = FALSE)
# 
# #Filter Maaslin results to include significant effect sizes
# results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# #Comparison of ZR and ZF
# omniresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Frozen", "Hot"))
# 
# ggplot(omniresults, aes(x=reorder(feature, -coef), y=coef)) +
#   geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
#   theme_bw() + 
#   labs(y = "Effect Size", x = "Genus", fill = "") +
#   coord_flip() + 
#   #scale_fill_manual(values = flexxt_palette) + 
#   theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")
# 
# ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# 
# 
# ##Maaslin2 analysis of NF v ZF v OF
# #Filter metadata to only include frozen temperature
# #benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# #rownames(benchmark_groups) <- benchmark_groups$sample
# #frozenresults <- filter(benchmark_groups, Temperature=="Frozen")
# 
# #Running Maaslin2
# #fit_data <- Maaslin2(
# #kraken_genus, frozenresults, here('DNA/3.maaslin'), transform = "NONE",
# ##Maaslin analysis on only frozen samples
# #Filter metdata to only include Frozen
# benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# rownames(benchmark_groups) <- benchmark_groups$sample
# benchmark_frozen <- filter(benchmark_groups, Temperature=="Frozen")
# 
# #Running Maaslin2
# fit_data <- Maaslin2(
#   kraken_genus, benchmark_frozen, here('DNA/3.maaslin/frozen'), transform = "NONE",
#   fixed_effects = c('Kit'),
#   random_effects = c('Donor'),
#   normalization = 'NONE',
#   standardize = FALSE,
#   min_abundance = 0.01,
#   min_prevalence = 0.1, 
#   reference = 'Kit,No Preservative')
# 
# # #Filter for significant and large effect sizes only for NF v ZF
# # results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# # nfzfresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Zymo"))
# # 
# # ggplot(nfzfresults, aes(x=reorder(feature, -coef), y=coef)) +
# 
# #Exporting Maaslin2 results
# write.csv(fit_data, here("DNA/3.maaslin/frozen/DNA_maaslin_frozen.csv"), row.names = FALSE)
# 
# #Filter Maaslin results to include significant effect sizes
# results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# #Comparison of ZR and ZF
# frozenresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Zymo"))
# 
# ggplot(frozenresults, aes(x=reorder(feature, -coef), y=coef)) +
#   geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
#   theme_bw() + 
#   labs(y = "Effect Size", x = "Genus", fill = "") +
#   coord_flip() + 
#   #scale_fill_manual(values = flexxt_palette) + 
#   theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")
# 
# # ggsave(here("outputs/figures/DNA_NFZF_maaslin2.pdf"), dpi=300, w=5, h=6)
# # ggsave(here("outputs/figures/DNA_NFZF_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# # 
# # #Filter for significant and large effect sizes only for NF v OF
# # results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
# # nfofresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Omnigene"))
# # 
# # ggplot(nfofresults, aes(x=reorder(feature, -coef), y=coef)) +
# #   geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
# #   theme_bw() + 
# #   labs(y = "Effect Size", x = "Genus", fill = "") +
# #   coord_flip() + 
# #   #scale_fill_manual(values = flexxt_palette) + 
# #   theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")
# # 
# # ggsave(here("outputs/figures/DNA_NFOF_maaslin2.pdf"), dpi=300, w=5, h=6)
# # ggsave(here("outputs/figures/DNA_NFOF_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# # 
# # 
# # ##Generating heatmap for Kit and Temp comparison
# # #Read in dataframes
# # zymoresults <- read.table(here("DNA/3.maaslin/zymo/all_results.tsv"), sep="\t", header=TRUE)
# # omniresults <- read.table(here("DNA/3.maaslin/omni/all_results.tsv"), sep="\t", header=TRUE)
# # frozenresults <- read.table(here("DNA/3.maaslin/frozen/all_results.tsv"), sep="\t", header=TRUE)
# # 
# # zymoresults <- filter(zymoresults, value=="Hot") %>% mutate(index=paste(metadata, "Zymo"))
# # omniresults <- filter(omniresults, value=="Hot") %>% mutate(index=paste(metadata, "Omni"))
# # frozenresults <- mutate(frozenresults, index=paste(metadata,value))
# # 
# # #Concatenate the dataframes
# # results <- rbind(zymoresults, frozenresults, omniresults)
# # write.table(results, here("DNA/3.maaslin/concat_results.tsv"), row.names = FALSE, sep="\t", quote=FALSE)
# # 
# # #Generating the heatmap
# # results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
# # results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl(".environmental", feature))
# # sigresults <- filter(results, abs(coef)>0.25)
# # uniquetaxa <- unique(sigresults$feature)
# # results <- filter(results, feature %in% uniquetaxa)
# # results <- mutate(results, coef2=ifelse(coef>3, 3, ifelse(coef < -3, -3, coef)))
# # limit <- max(abs(results$coef2)) * c(-1, 1)
# # 
# # 
# # 
# # ggplot(results, aes(x=index, y=feature, fill=coef2)) + 
# ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_OFOH_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
# 
# ##Generating heatmap for Maaslin comparison of Kit and Temperature
# #Reading in dataframes and adding in column of comparison to each Maaslin analysis
# frozenresults <- read.table(here("DNA/3.maaslin/frozen/all_results.tsv"), sep="\t", header=TRUE)
# omniresults <- read.table(here("DNA/3.maaslin/omni/all_results.tsv"), sep="\t", header=TRUE)
# zymoresults <- read.table(here("DNA/3.maaslin/zymo/all_results.tsv"), sep="\t", header=TRUE)
# 
# frozenresults <- mutate(frozenresults, index=paste(metadata,value))
# omniresults <- filter(omniresults, value=="Hot") %>% mutate(index=paste(metadata,"Omni"))
# zymoresults <- filter(zymoresults, value=="Hot") %>% mutate(index=paste(metadata,"Zymo"))
# 
# #Concatenate the dataframes together
# results <- rbind(frozenresults,omniresults,zymoresults)
# 
# #Making the heatmap
# results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
# limit <- max(abs(results$coef)) * c(-1, 1)
# results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl("environmental.samples", feature))
# results <- mutate(results, coef2=ifelse(coef<3, coef, 4))
# 
# 
# ggplot(results, aes(x=index, y=feature, fill=coef)) + 
#   geom_tile() + 
#   theme_bw() +
#   scale_x_discrete(position = "top") +
#   facet_grid(~metadata, scales="free_x") + 
#   # scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
#   # geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
#   # scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
#   # theme(panel.grid = element_blank(),
#   #       panel.border = element_blank()) + 
#   # theme(axis.text.y = element_text(size = 15)) + 
#   #scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
#   geom_point(aes(alpha=signif), color = "white", size = 3, show.legend = F) + 
#   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank()) + 
#   theme(axis.text.y = element_text(size = 12)) + 
#   theme(axis.text.x = element_text(size = 15)) + 
#   theme(axis.title.y = element_blank()) + 
#   theme(axis.title.x = element_blank()) + 
#   theme(axis.ticks.x = element_blank()) + 
#   theme(
#     strip.text.x = element_text(
#       size = 18, color = "black", face = "bold"
#     ),
#     strip.background = element_rect(
#       color="white", fill="white", size=1.5, linetype="solid"
#     ),
#     strip.placement = 'outside'
#   ) 
# 
# ggsave(here("outputs/figures/DNA_heatmap.pdf"), dpi=300, w=12, h=12)
# ggsave(here("outputs/figures/DNA_heatmap.jpeg"), dpi=300, w=12, h=12)
# 
# ggplot(results, aes(x=coef)) + 
#   geom_histogram()+
#   xlim(c(-1,1))
# 
# #Generating the heatmap with gram stain and phyla
# results <- read.table(here("DNA/3.maaslin/concat_results.tsv"), sep="\t", header=TRUE)
# results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
# results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl(".environmental", feature))
# sigresults <- filter(results, abs(coef)>0.25)
# uniquetaxa <- unique(sigresults$feature)
# results <- filter(results, feature %in% uniquetaxa)
# results <- mutate(results, coef2=ifelse(coef>3, 3, ifelse(coef < -3, -3, coef)))
# limit <- max(abs(results$coef2)) * c(-1, 1)
# #Merge with gram stain and phyla dataframe
# phylagram <- read.csv(here("outputs/tables/Maaslin_heatmap_126.csv"), sep=",", header=TRUE)
# phylagram <- merge(results, phylagram, by='feature')
# 
# #Plotting phyla heatmap
# ggplot(phylagram, aes(x=index, y=feature, fill=coef2)) + 
#   geom_tile() + 
#   theme_bw() +
#   scale_x_discrete(position = "top") +
#   facet_grid(Phylum~metadata, scales="free", space="free_y") + 
#   scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
#   geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
#   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank()) + 
#   theme(axis.text.y = element_text(size = 13)) + 
#   theme(axis.text.x = element_text(size = 13)) + 
#   theme(axis.title.y = element_blank()) + 
#   theme(axis.title.x = element_blank()) + 
#   theme(axis.ticks.x = element_blank()) + 
#   theme(
#     strip.text.x = element_text(
#       size = 18, color = "black", face = "bold"
#     ),
#     strip.background = element_rect(
#       color="white", fill="white", size=1.5, linetype="solid"
#     ),
#     strip.placement = 'outside'
#   ) 
# 
# ggsave(here("outputs/figures/DNA_heatmap_phyla.pdf"), dpi=300, w=12, h=14)
# ggsave(here("outputs/figures/DNA_heatmap_phyla.jpeg"), dpi=300, w=12, h=14)
# 
# #Plotting gram stain heatmap
# ggplot(phylagram, aes(x=index, y=feature, fill=coef2)) + 
#   geom_tile() + 
#   theme_bw() +
#   scale_x_discrete(position = "top") +
#   facet_grid(GramStatus~metadata, scales="free", space="free_y") + 
#   scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
#   geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
#   scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank()) + 
#   theme(axis.text.y = element_text(size = 13)) + 
#   theme(axis.text.x = element_text(size = 13)) + 
#   theme(axis.title.y = element_blank()) + 
#   theme(axis.title.x = element_blank()) + 
#   theme(axis.ticks.x = element_blank()) + 
#   theme(
#     strip.text.x = element_text(
#       size = 18, color = "black", face = "bold"
#     ),
#     strip.background = element_rect(
#       color="white", fill="white", size=1.5, linetype="solid"
#     ),
#     strip.placement = 'outside'
#   ) 
# 
# ggsave(here("outputs/figures/DNA_heatmap_gram.pdf"), dpi=300, w=12, h=14)
# ggsave(here("outputs/figures/DNA_heatmap_gram.jpeg"), dpi=300, w=12, h=12)
##### PHYLUM LEVEL DIFFERENTIAL ABUNDANCE: MAASLIN2 MODEL #####
# library(Maaslin2)
# 
# #Grouped analysis: Phylum-level
# #Read in classification and metadata files
# kraken_phylum <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE))
# 
# benchmark_groups <- read.csv(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
# rownames(benchmark_groups) <- benchmark_groups$sample
# 
# #Running Maaslin2
# fit_data <- Maaslin2(
#   kraken_phylum, benchmark_groups, here('DNA/3.maaslin/complete'), transform = "NONE",
#   fixed_effects = c('Condition'),
#   random_effects = c('Donor'),
#   normalization = 'NONE',
#   standardize = FALSE,
#   min_abundance = 0.01,
#   min_prevalence = 0.1, 
#   reference = c(('Condition,NF')))
# 
# #Read in Maaslin results
# results <- read.table(here("DNA/3.maaslin/complete/all_results.tsv"), sep="\t", header=TRUE)
# 
# #Filter Maaslin results to include significant effect sizes
# results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
# limit <- max(abs(results$coef)) * c(-1, 1)
# results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl("environmental.samples", feature))
# results <- mutate(results, coef2=ifelse(coef<3, coef, 4))
# allresults <- mutate(results, index=paste(metadata,value))
# 
#  ggplot(allresults, aes(x=index, y=feature, fill=coef2)) + 
#    geom_tile() + 
#    theme_bw() +
#    scale_x_discrete(position = "top") +
#    facet_grid(~metadata, scales="free", space="free_y") + 
#    scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
#    geom_point(aes(alpha=signif), color = "dark grey", size = 5, show.legend = F) + 
#    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
#    theme(panel.grid = element_blank(),
#          panel.border = element_blank()) + 
#    theme(axis.text.y = element_text(size = 13)) + 
#    theme(axis.text.x = element_text(size = 13)) + 
#    theme(axis.title.y = element_blank()) + 
#    theme(axis.title.x = element_blank()) + 
#    theme(axis.ticks.x = element_blank()) + 
#    theme(
#      strip.text.x = element_text(
#        size = 18, color = "black", face = "bold"
#      ),
#      strip.background = element_rect(
#        color="white", fill="white", size=1.5, linetype="solid"
#      ),
#      strip.placement = 'outside'
#    ) 
# 
# ggsave(here("outputs/figures/DNA_groupedheatmap_phylum.pdf"), dpi=300, w=10, h=8)
# ggsave(here("outputs/figures/DNA_groupedheatmap_phylum.jpeg"), dpi=300, w=10, h=8)

##### READ KMER DIVERSITY ####
k51 <- read.table(here("RNA/04_sourmash/04_sourmash_compare/compare_k51.csv"), sep=",", header=TRUE)
colnames(k51) <- gsub(".*D", "D", colnames(k51))
colnames(k51) <- gsub("_con.*", "", colnames(k51))
rownames(k51) <- names(k51)
k51 <- cbind(Sample = rownames(k51), k51)
rownames(k51) <- NULL

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

data_long <- melt(k51, id.vars=c("Sample")) # reshape to long format
names(data_long) <- c("Sample1", "Sample2", "Identity")

data_long <- data_long %>% separate(Sample1, c("Donor1", "Condition1", "Replicate1"), remove=FALSE)
data_long <- data_long %>% separate(Sample2, c("Donor2", "Condition2", "Replicate2"), remove=FALSE)

# make a table for within-sample comparisons
precision_table <- data_long %>% filter(Donor1 == Donor2) %>%
                                 filter(Condition1 == Condition2) %>% 
                                 filter(Replicate1 != Replicate2)

# filter table to remove double measures
precision_table <- mutate(precision_table, uniqueID=ifelse(as.character(Sample1) > as.character(Sample2), paste(Sample1, Sample2), paste(Sample2, Sample1)))
precision_table <- precision_table %>% group_by(uniqueID) %>% filter(row_number() == 1)

# new table calculating averages within donor/condition combo
precision_values <- precision_table %>% group_by(Donor1, Condition1) %>% 
  summarise_at("Identity", list(Identity= mean))

precision_values <- mutate(precision_values, PlotOrder=ifelse(Condition1 == "NF", 1, 
                                                        ifelse(Condition1 == "OF", 2, 
                                                               ifelse(Condition1 == "OR", 3, 
                                                                      ifelse(Condition1 == "OH", 4, 
                                                                             ifelse(Condition1 == "ZF", 5,
                                                                                    ifelse(Condition1 == "ZR", 6, 7)))))))

ggplot(precision_values, aes(x=reorder(Condition1, PlotOrder), y=Identity)) + 
  geom_jitter(width=0.1) +
  geom_boxplot(alpha=0.5, outlier.shape=NA, aes(fill=Condition1)) + 
  theme_bw() + 
  scale_fill_manual(values=condition_palette) + 
  xlab("Condition") + 
  ylab("Similarity") + 
  ggtitle("Within-Condition Similarity (k = 51)") + 
  stat_compare_means(method = "t.test", comparisons=list(c("OF", "OR"), c("ZF", "ZR"), c("ZF", "ZH"), c("NF", "OF"), c("NF", "ZF")), 
                     tip.length = 0, label="p.signif") + 
  theme(legend.position = "none")

ggsave(here("outputs/figures/RNA_kmerdistance.pdf"), dpi=300, w=10, h=6)
ggsave(here("outputs/figures/RNA_kmerdistance.jpeg"), dpi=300, w=10, h=6)
