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

##### EDITING METADATA ######
#Split SampleID into donor, condition, replicate 
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

##### SPECIES-LEVEL RELATIVE ABUNDANCE #####
# Read in species percentage bracken data
species <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_species_percentage.txt"), sep="\t", header=TRUE)

#Color palette
n_taxa <- 30
myCols <- colorRampPalette(brewer.pal(9, "Set1"))
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

abundance_threshold <- sort(rowSums(species), decreasing = T)[n_taxa]
bracken_plot <- species[rowSums(species) >= abundance_threshold,]
bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  100 - colSums(bracken_plot))))

bracken_plot$Species <- row.names(bracken_plot)
bracken_plot$Species <- gsub("\\(miscellaneous\\)", "", bracken_plot$Species)
bracken_long <- melt(bracken_plot, id.vars = "Species", variable.name = "Sample", value.name = "rel_abundance")
bracken_long <- mutate(bracken_long, Sample=gsub("\\.", "_", Sample))

# Merging the metadata with the bracken_long
colnames(metadata)[3]<-"Sample"
bracken_pheno <- merge(bracken_long, metadata, by = "Sample")
#bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))

# set factor in correct plotting order
bracken_pheno$Species <- factor(bracken_pheno$Species, levels = bracken_plot$Species)

#Stacked bar plot for all samples
plot_bracken <- function(counts, title){
  g <- ggplot(counts, aes(x=Sample, y=rel_abundance, fill=Species)) +
    geom_bar(stat="identity") +
    labs(
      title = title,
      x = "",
      y = "Relative Abundance (%)"
    ) +
    scale_fill_manual(values = barplot_pal) +
    guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
    #theme_cowplot(12) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "plain", size = 14),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))
  
  return(g)
}

plot_bracken(bracken_pheno, "Species-Level Relative Abundance") + theme(legend.position = "none")
#ggsave(here("outputs/figures/OverallDNA_species_abundance.pdf"), dpi=300, w=20, h=4)

##### GENUS-LEVEL RELATIVE ABUNDANCE #####
# Read in genus percentage bracken data
genus <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)

#Color palate
n_taxa <- 20
myCols <- colorRampPalette(brewer.pal(9, "Set1"))
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))

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

# plot in order of decreasing relative abundance of desired taxon
#bracken_pheno$label <- factor(bracken_pheno$label, levels = c("Donor1 F1", "Donor1 F2","Donor1 F3","Donor1 FO1","Donor1 FO2","Donor1 FO3", "Donor1 C1","Donor1 C2","Donor1 C3","Donor1 CL1","Donor1 CL2","Donor1 CL3", "Donor1 H1","Donor1 H2","Donor1 H3", "Donor2 F1","Donor2 F2","Donor2 F3", "Donor2 FO1", "Donor2 FO2","Donor2 FO3","Donor2 C1","Donor2 C2","Donor2 C3", "Donor2 H1", "Donor2 H2", "Donor2 H3"))

plot_bracken <- function(counts, title){
  g <- ggplot(counts, aes(x=Sample, y=rel_abundance, fill=Genus)) +
    geom_bar(stat="identity") +
    labs(
      title = title,
      x = "",
      y = "Relative Abundance (%)"
    ) +
    scale_fill_manual(values = barplot_pal) +
    guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
    #theme_cowplot(12) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "plain", size = 14),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))
  
  return(g)
}

plot_bracken(bracken_pheno, "Genus-Level Relative Abundance") #+ theme(legend.position = "none")

#ggsave(here("outputs/figures/OverallDNA_genus_abundance.pdf"), dpi=300, w=20, h=5)

#####GENUS-LEVEL FACET BY DONOR######
metadata <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
genus <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)



#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
metadata <- metadata %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
#Modify Condition column, so that anything labeled with B# is changed to Controls
metadata <- mutate(metadata, Condition=ifelse(Condition %in% c("B1", "B2", "B3", "B4"), "Controls", Condition))
#Within the DNA dataframe and Condition/Donor column, factor() alters the sorting of the variables in Condition/Donor - does not change the data frame
metadata$Donor <- factor(metadata$Donor, levels = c("NCO", "PCO", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
metadata <- mutate(metadata, Label=ifelse(Replicate == "R2", Condition, ""))
metadata$Condition <- factor(metadata$Condition, levels = c("Controls", "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))

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

#Color pal
n_taxa <- 20
myCols <- colorRampPalette(brewer.pal(9, "Set1")) # WAS Set1
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

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
  scale_y_continuous(limits = c(-3, 100.1), expand = c(0, 0))  + 
  facet_wrap(~DonorLabel, ncol = 5, scales = "free") + 
  theme(legend.position = "none") + 
  new_scale_fill() + 
  geom_tile(aes(x=Sample, y = -2, fill = Condition)) + 
  scale_fill_manual(values = condition_palette)

ggsave(here("outputs/figures/DNA_facetedgenusabundance_nicelabels.jpg"), dpi=300, w = 12, h = 5)
ggsave(here("outputs/figures/DNA_facetedgenusabundance_nicelabels.pdf"), dpi=300, w = 12, h = 5)


##### BRAY CURTIS PCoA #####
# Read in count data
cts <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_species_reads.txt"), sep="\t", header = TRUE) 
cts <- cts %>% select(!starts_with("NCO") & !starts_with("PCO") & !matches("D03.ZH.R2"))

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

#ggsave(here("outputs/figures/DNA_MDS_Scree.pdf"), dpi=300, w=7, h=5)

#Plot
mds_data <- data.frame(Sample = gsub("\\.", "_", rownames(mds_values)),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
metadata <- filter(metadata, Donor!="NCO" & Donor!="PCO")
mds_data <-rename (mds_data, Sample=SampleID)
mds_meta <- merge(mds_data, metadata, by = "Sample")

maicolors <- paletteer_d("ggthemes::Tableau_10")
mds_plot <- ggplot(mds_meta, aes(x, y)) +
  stat_ellipse(aes(color = Donor), type = 't', size = 0.5, show.legend = F) +
  geom_point(size = 1, alpha = 0.7, aes(color=Donor)) +
  scale_color_manual(values = maicolors) +
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

#ggsave(here("outputs/figures/DNASpecies_Bray_Curtis.pdf"), dpi=300, w=5, h=6)

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

#ggsave(here("outputs/figures/DNASpecies_Bray_Curtis_Preservation.pdf"), dpi=300, w=5, h=6)

#####Shannon Diversity Plots#####
benchmark_s <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_species_percentage.txt"), sep="\t")
benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
benchmark_groups <- mutate(benchmark_groups, sample=gsub("_", ".", SampleID))
benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")

shannon_div_s <- diversity(t(benchmark_s), index = "shannon")
div <- data.frame("shannon_div" = shannon_div_s, "sample" = names(shannon_div_s))
div_meta <- merge(div, benchmark_groups, by = "sample")

maicolors <- paletteer_d("ggthemes::Tableau_10")
div_plot <- ggplot(div_meta, aes(Donor, shannon_div)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.6),
              alpha = 0.85, aes(fill = Donor), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Donor), alpha = 0.5) +
  labs(x = "",
       y = "Shannon Diversity",
       fill = "") +
  scale_color_manual(values = maicolors) +
  scale_y_continuous(limits = c(2, 6.6)) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) 
div_plot

# ggsave(here("outputs/figures/DNA_alphadiv.pdf"), dpi=300, w=5, h=5)
# ggsave(here("outputs/figures/DNA_alphadiv.jpg"), dpi=300, w=5, h=5)

#####Attempt at Maaslin2 for diff abundance#####
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.12") 
BiocManager::install("Maaslin2", version="3.12")

library(Maaslin2)

#Other installation method to be able to use References
#install.packages("devtools")
# library(devtools)
devtools::install_github("biobakery/maaslin2")
# 
# library(Maaslin2)
# to get help on usage of Maaslin2 ?Maaslin2

#Reading in metadata dataframe and creating metadata file for Maaslin2
benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
benchmark_groups <- mutate(benchmark_groups, sample=gsub("_", ".", SampleID))
benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")
benchmark_groups <- mutate(benchmark_groups, Kit=substr(Condition,1,1))
benchmark_groups <- mutate(benchmark_groups, Temperature=substr(Condition,2,2))
benchmark_groups <- mutate(benchmark_groups, Kit=ifelse(Kit=="N", "No Preservative", ifelse(Kit=="O", "Omnigene", "Zymo")))
benchmark_groups <- mutate(benchmark_groups, Temperature=ifelse(Temperature=="F", "Frozen", ifelse(Temperature=="R", "Room Temp", "Hot")))
rownames(benchmark_groups) <- benchmark_groups$sample

#Exporting metdata
write.table(benchmark_groups, here("outputs/tables/DNA_metadata.tsv"), row.names = FALSE, sep="\t", quote=FALSE)  

#Transpose the classification file
kraken_genus <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE))

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, benchmark_groups, here('DNA/3.maaslin'), transform = "NONE",
  fixed_effects = c('Kit', 'Temperature'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Kit,No Preservative;Temperature,Frozen')

#from Dylan
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#results <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Illumina DNA Prep", "Nextera XT"))
zymoresults <- filter(results, value=='Zymo')
zymoresults <- mutate(zymoresults, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Zymo"))
#results <- mutate(results, feature=gsub(" miscellaneous", "", gsub("\\.", " ", feature)))

ggplot(zymoresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

# ggsave(here("outputs/figures/DNA_NFZymo_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_NFZymo_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)

#Omnigene vs NF
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#results <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Illumina DNA Prep", "Nextera XT"))
omniresults <- filter(results, value=='Omnigene')
omniresults <- mutate(omniresults, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Omni"))
#results <- mutate(results, feature=gsub(" miscellaneous", "", gsub("\\.", " ", feature)))

ggplot(omniresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

# ggsave(here("outputs/figures/DNA_NFOmni_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_NFOmni_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)

#Temperature 
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
#results <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Illumina DNA Prep", "Nextera XT"))
hotresults <- filter(results, value=='Hot')
hotresults <- mutate(hotresults, DirectionEnriched = ifelse(coef < 0, "Frozen", "Hot"))
#results <- mutate(results, feature=gsub(" miscellaneous", "", gsub("\\.", " ", feature)))

ggplot(hotresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

# ggsave(here("outputs/figures/DNA_NFOmni_maaslin2.pdf"), dpi=300, w=5, h=6)
# ggsave(here("outputs/figures/DNA_NFOmni_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)
