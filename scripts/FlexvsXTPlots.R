library(ggplot2)
library(tidyverse)
library(paletteer) 
library(ggsignif)
library(here)
library(RColorBrewer)
library(reshape2)

library(metagenomeSeq) # need to figure out how to install
library(devtools)
library(vegan)
library(DESeq2)
library(genefilter)

n_taxa <- 30
myCols <- colorRampPalette(brewer.pal(9, "Set1"))
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"


##### SPECIES-LEVEL RELATIVE ABUNDANCE #####
# Read in species percentage bracken data
species <- read.csv(here("XTvFlex/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_species_percentage.txt"), sep="\t", header=TRUE)

# read in sample metadata; change replicate letters to numbers
groupings <- read.table(here("XTvFlex/groups_and_conditions.tsv"), sep="\t", header = TRUE)
groupings <- mutate(groupings, Replicate=ifelse(Replicate == 'B', 2, ifelse(Replicate == 'C', 3, ifelse(Replicate == "A", 1, 4))))
groupings <- mutate(groupings, groupedID=ifelse(Replicate == 1, "XT 1", ifelse(Replicate == '2', "XT 2", ifelse(Replicate == '3', "XT 3", "Flex"))))


abundance_threshold <- sort(rowSums(species), decreasing = T)[n_taxa]
bracken_plot <- species[rowSums(species) >= abundance_threshold,]
bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  100 - colSums(bracken_plot))))

bracken_plot$Species <- row.names(bracken_plot)
bracken_plot$Species <- gsub("\\(miscellaneous\\)", "", bracken_plot$Species)
bracken_long <- melt(bracken_plot, id.vars = "Species", variable.name = "Sample", value.name = "rel_abundance")
bracken_long <- mutate(bracken_long, Sample=gsub("\\.", "-", Sample))

# merge in pheno data
bracken_pheno <- merge(bracken_long, groupings, by = "Sample")
bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))

# set factor in correct plotting order
bracken_pheno$Species <- factor(bracken_pheno$Species, levels = bracken_plot$Species)

plot_bracken <- function(counts, title){
  g <- ggplot(counts, aes(x=label, y=rel_abundance, fill=Species)) +
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

plot_bracken_facet <- function(counts, title){
  g <- ggplot(counts, aes(x=groupedID, y=rel_abundance, fill=Species)) +
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
    scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))  + 
    facet_wrap(~Donor, ncol = 5) 
  
  return(g)
}
plot_bracken_facet(bracken_pheno, "Species-Level Relative Abundance") 
ggsave(here("outputs/figures/XTvFlex_species_composition_facet.pdf"), dpi=300, w=9, h=5)
ggsave(here("outputs/figures/XTvFlex_species_composition_facet.jpg"), dpi=300, w=9, h=5)


##### GENUS-LEVEL RELATIVE ABUNDANCE #####

n_taxa <- 20
myCols <- colorRampPalette(brewer.pal(9, "Set1"))
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

# Read in species percentage bracken data
genus <- read.csv(here("XTvFlex/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)

# read in sample metadata; change replicate letters to numbers
groupings <- read.table(here("XTvFlex/groups_and_conditions.tsv"), sep="\t", header = TRUE)
groupings <- mutate(groupings, Replicate=ifelse(Replicate == 'B', 2, ifelse(Replicate == 'C', 3, ifelse(Replicate == "A", 1, 4))))
groupings <- mutate(groupings, groupedID=ifelse(Replicate == 1, "XT 1", ifelse(Replicate == '2', "XT 2", ifelse(Replicate == '3', "XT 3", "Flex"))))


abundance_threshold <- sort(rowSums(genus), decreasing = T)[n_taxa]
bracken_plot <- genus[rowSums(genus) >= abundance_threshold,]
bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  100 - colSums(bracken_plot))))

bracken_plot$Genus <- row.names(bracken_plot)
bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", bracken_plot$Genus)
bracken_long <- melt(bracken_plot, id.vars = "Genus", variable.name = "Sample", value.name = "rel_abundance")
bracken_long <- mutate(bracken_long, Sample=gsub("\\.", "-", Sample))

# merge in pheno date
bracken_pheno <- merge(bracken_long, groupings, by = "Sample")
bracken_pheno <- mutate(bracken_pheno, label=paste(Donor, groupedID))

# set factor in correct plotting order
bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)

# plot in order of decreasing relative abundance of desired taxon
#bracken_pheno$label <- factor(bracken_pheno$label, levels = c("Donor1 F1", "Donor1 F2","Donor1 F3","Donor1 FO1","Donor1 FO2","Donor1 FO3", "Donor1 C1","Donor1 C2","Donor1 C3","Donor1 CL1","Donor1 CL2","Donor1 CL3", "Donor1 H1","Donor1 H2","Donor1 H3", "Donor2 F1","Donor2 F2","Donor2 F3", "Donor2 FO1", "Donor2 FO2","Donor2 FO3","Donor2 C1","Donor2 C2","Donor2 C3", "Donor2 H1", "Donor2 H2", "Donor2 H3"))

plot_bracken <- function(counts, title){
  g <- ggplot(counts, aes(x=label, y=rel_abundance, fill=Genus)) +
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

plot_bracken(bracken_pheno, "Genus-Level Relative Abundance") + theme(legend.position = "none")

ggsave(here("outputs/figures/XTvFlex_genus_composition.pdf"), dpi=300, w=7, h=5)
ggsave(here("outputs/figures/XTvFlex_genus_composition.jpg"), dpi=300, w=7, h=5)


# to do
plot_bracken_facet <- function(counts, title){
  g <- ggplot(counts, aes(x=groupedID, y=rel_abundance, fill=Genus)) +
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
    scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))  + 
    facet_wrap(~Donor, ncol = 5) 
  
  return(g)
}
plot_bracken_facet(bracken_pheno, "Genus-Level Relative Abundance") 
ggsave(here("outputs/figures/XTvFlex_genus_composition_facet.pdf"), dpi=300, w=9, h=5)
ggsave(here("outputs/figures/XTvFlex_genus_composition_facet.jpg"), dpi=300, w=9, h=5)


##### BRAY CURTIS #####
## css normalize
cts <- read.csv(here("XTvFlex/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_species_reads.txt"), sep="\t", header = TRUE) # read in count data

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

ggsave(here("outputs/figures/mds_scree.pdf"), dpi=300, w=7, h=5)

# plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
# mds_meta <- merge(mds_data, groupings, by = "sample")
# mds_meta$group <- factor(mds_meta$group, levels = c("D1F", "D1FO","D1C","D1CL","D1H","D2F","D2FO", "D2C","D2H"))
# 
# reds <- brewer.pal(5, name = "Reds")
# blues <- brewer.pal(4, name = "Blues")
# mycolors <- c(reds, blues)
# 
# mds_meta <- mutate(mds_meta, label=group)
# 
# # lengthen labels
# mds_meta$label <- gsub("D1", "Donor1 ", mds_meta$label)
# mds_meta$label <- gsub("D2", "Donor2 ", mds_meta$label)
# mds_meta$label <- gsub("C", "Cold", mds_meta$label)
# mds_meta$label <- gsub("L", ", Long BB", mds_meta$label)
# mds_meta$label <- gsub("F", "Fresh", mds_meta$label)
# mds_meta$label <- gsub("O", ", OMNIgene", mds_meta$label)
# mds_meta$label <- gsub("H", "Hot", mds_meta$label)
# mds_meta$label <- factor(mds_meta$label, levels = c("Donor1 Fresh", "Donor1 Fresh, OMNIgene", "Donor1 Cold", "Donor1 Cold, Long BB", "Donor1 Hot", "Donor2 Fresh", "Donor2 Fresh, OMNIgene", "Donor2 Cold", "Donor2 Hot"))
# 
# 
# 
# mds_plot <- ggplot(mds_meta, aes(x, y, fill = label)) +
#   geom_point(size = 2, alpha = 0.95, pch = 21) +
#   scale_fill_manual(values = mycolors) +
#   labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
#        y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
#        title = "Species-Level Bray Curtis",
#        color = "") +
#   theme_classic() +
#   background_grid() +
#   coord_fixed() +
#   theme(legend.position = "bottom", legend.title = element_blank()) + 
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))
# 
# mds_plot
# ggsave("/Users/dylanmaghini/scg4/projects/africa/01_benchmarking/02_classification/plots/braycurtis_mds.jpg", width = 8, height = 5)