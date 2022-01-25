library(ggplot2)
library(tidyverse)
library(paletteer) 
library(ggsignif)
library(here)
library(RColorBrewer)
library(reshape2)
library(metagenomeSeq)
library(devtools)
library(vegan)
library(DESeq2)
library(genefilter)
library(cowplot)
library(ggpubr)
library(PairedData)
library(forcats)

# save a paletteer palette, and subset it to the first ten values
flexxt_palette <- c("#6c7dff", "#ff8c6c")

# add "names" to the colors in the palette, corresponding to the donor names
names(flexxt_palette) <- c("Illumina DNA Prep", "Nextera XT")

##### SPECIES-LEVEL RELATIVE ABUNDANCE #####
n_taxa <- 30
myCols <- colorRampPalette(brewer.pal(9, "Set1"))
barplot_pal <- myCols(n_taxa)
barplot_pal <- sample(barplot_pal)
barplot_pal[n_taxa + 1] <- "gray"

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
#ggsave(here("outputs/figures/XTvFlex_species_composition_facet.pdf"), dpi=300, w=9, h=5)
#ggsave(here("outputs/figures/XTvFlex_species_composition_facet.jpg"), dpi=300, w=9, h=5)


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

#ggsave(here("outputs/figures/XTvFlex_genus_composition.pdf"), dpi=300, w=7, h=5)
#ggsave(here("outputs/figures/XTvFlex_genus_composition.jpg"), dpi=300, w=7, h=5)


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
#ggsave(here("outputs/figures/XTvFlex_genus_composition_facet.pdf"), dpi=300, w=9, h=5)
#ggsave(here("outputs/figures/XTvFlex_genus_composition_facet.jpg"), dpi=300, w=9, h=5)


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

#ggsave(here("outputs/figures/mds_scree.pdf"), dpi=300, w=7, h=5)

# plot
mds_data <- data.frame(Sample = gsub("\\.", "-", rownames(mds_values)),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, groupings, by = "Sample")

mds_plot <- ggplot(mds_meta, aes(x, y)) +
  geom_point(size = 3, alpha = 0.7, aes(shape=Method, color=Donor)) +
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

#ggsave(here("outputs/figures/xt_braycurtis_mds.pdf"), dpi=300, w=5, h=6)
#ggsave(here("outputs/figures/xt_braycurtis_mds.jpg"), dpi=300, w=5, h=6)



##### DIFFERENTIAL ABUNDANCE #####


### NOTE : THIS ISN'T PAIRED, BUT WITH ALL XT SAMPLES ###
count_data <- read.csv(here("XTvFlex/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_genus_reads.txt"), sep="\t", header = TRUE) # read in count data
head(count_data)
row.names(count_data)
metaData <- groupings <- read.table(here("XTvFlex/groups_and_conditions.tsv"), sep="\t", header = TRUE)
metaData <- mutate(metaData, Sample=gsub("-", ".", Sample))
head(metaData)

# filter metadata to only include first XT replicates for each sample
metaData <- metaData %>% filter(Method == "Flex" | Replicate == "A")
metaData <- mutate(metaData, Donor=as.factor(Donor))
metaData <- mutate(metaData, Metehod=as.factor(Method))

count_data <- count_data %>% dplyr::select(!ends_with("B")) %>% dplyr::select(!ends_with("C"))


dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metaData, design = ~ Method)
dds <- estimateSizeFactors(dds, type = "poscounts")
idx <- genefilter(counts(dds, normalized = TRUE), pOverA(0.2, 500))  # 20% need to be over 500
dds <- dds[idx, ]
dds <- DESeq(dds, parallel = T)
res <- results(dds, name="Method_XT_vs_Flex", alpha = 0.05)
resultsNames(dds)

resOrdered <- res %>%
  as_tibble(rownames = "genus") %>%
  arrange(pvalue) %>%
  mutate(group = ifelse(log2FoldChange < 0, "Flex", "XT"))

res_df_filt <- resOrdered %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  arrange(log2FoldChange)
res_df_filt$genus <- factor(res_df_filt$genus,
                            levels = as.character(res_df_filt$genus))

deseq_genera <- ggplot(res_df_filt, aes(log2FoldChange, genus, fill = group)) +
  #geom_point(size = 2, color = "black", pch = 21) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_cowplot() +
  theme(axis.text.y = element_text(face = "italic")) +
  #scale_fill_manual(values = hotcold) +
  # scale_x_continuous(breaks = seq(-2, 2, 1)) +
  labs(
    x = "Log2 Fold Change",
    y = "Genus",
    fill = ""
  ) +
  background_grid() +
  theme(legend.position = "top",
        legend.justification = "center")

deseq_genera
ggsave(here("outputs/figures/XTvFlex_genus_deseq.pdf"), dpi=300, w=6, h=5)

### TRYIING DESEQ WITH PAIRED ANALYSIS
count_data <- read.csv(here("XTvFlex/02_kraken2_classification/processed_results_krakenonly/taxonomy_matrices/kraken_genus_reads.txt"), sep="\t", header = TRUE) # read in count data
head(count_data)
row.names(count_data)
metaData <- groupings <- read.table(here("XTvFlex/groups_and_conditions.tsv"), sep="\t", header = TRUE)
metaData <- mutate(metaData, Sample=gsub("-", ".", Sample))


# filter metadata and count data to only include first XT replicates for each sample
metaData <- metaData %>% filter(Method == "Flex" | Replicate == "A")
metaData <- mutate(metaData, Donor=as.factor(Donor))
metaData <- mutate(metaData, Method=as.factor(Method))
count_data <- count_data %>% dplyr::select(!ends_with("B")) %>% dplyr::select(!ends_with("C"))


dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metaData, design = ~ Donor + Method)
dds <- estimateSizeFactors(dds, type = "poscounts")
idx <- genefilter(counts(dds, normalized = TRUE), pOverA(0.2, 500))  # 20% need to be over 500
dds <- dds[idx, ]
dds <- DESeq(dds, parallel = T)
res <- results(dds, name="Method_XT_vs_Flex", alpha = 0.05)
resultsNames(dds)

resOrdered <- res %>%
  as_tibble(rownames = "genus") %>%
  arrange(pvalue) %>%
  mutate(group = ifelse(log2FoldChange < 0, "Flex", "XT"))

res_df_filt <- resOrdered %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  arrange(log2FoldChange)
res_df_filt$genus <- factor(res_df_filt$genus,
                            levels = as.character(res_df_filt$genus))

deseq_genera <- ggplot(res_df_filt, aes(log2FoldChange, genus, fill = group)) +
  #geom_point(size = 2, color = "black", pch = 21) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_cowplot() +
  theme(axis.text.y = element_text(face = "italic")) +
  #scale_fill_manual(values = hotcold) +
  # scale_x_continuous(breaks = seq(-2, 2, 1)) +
  labs(
    x = "Log2 Fold Change",
    y = "Genus",
    fill = ""
  ) +
  background_grid() +
  theme(legend.position = "top",
        legend.justification = "center")

deseq_genera
#ggsave(here("outputs/figures/XTvFlex_genus_deseq_paired.pdf"), dpi=300, w=6, h=5)


### tryig to shrink log fold change and plot
plotMA(res, ylim=c(-2,2), alpha = 0.001, xlab="Mean of Normalized Counts", ylab="Log2 Fold Change")


##### ALPHA DIVERSITY #####
# Species level

benchmark_s <- read.table(here("XTvFlex/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_species_percentage.txt"), sep="\t")
benchmark_groups <- read.table(here("XTvFlex/groups_and_conditions.tsv"), sep="\t", header = TRUE)
benchmark_groups <- mutate(benchmark_groups, sample=gsub("-", ".", Sample))

benchmark_groups <- benchmark_groups %>% filter(Method == "Flex" | Replicate == "A")
benchmark_groups <- mutate(benchmark_groups, Donor=as.factor(Donor))
benchmark_groups <- mutate(benchmark_groups, Metehod=as.factor(Method))

benchmark_s <- benchmark_s %>% select(!ends_with("B")) %>% select(!ends_with("C"))

shannon_div_s <- diversity(t(benchmark_s), index = "shannon")
div <- data.frame("shannon_div" = shannon_div_s, "sample" = names(shannon_div_s))
div_meta <- merge(div, benchmark_groups, by = "sample")

# pval with wilcox
pval <- compare_means(shannon_div ~ Method, data = div_meta, method = "wilcox.test", p.adjust.method = "fdr")
pval

# trying with fisher t test
pval <- compare_means(shannon_div ~ Method, data = div_meta, method = "t.test", p.adjust.method = "fdr")
pval

## trying paired plots
ggpaired(div_meta, x="Method", y="shannon_div", fill="Method", line.color="gray", id="Donor") + 
  stat_compare_means(paired=TRUE) + 
  theme_cowplot() + 
  ylim(c(4,5)) +
  labs(x = "", 
       y = "Shannon Diveristy", 
       title = "Shannon Diversity (Species)", 
       fill = "") +
  theme(axis.title.x = element_blank()) + theme(legend.position="none")

ggsave(here("outputs/figures/XTvFlex_alphadiv.jpg"), dpi=300, w=5, h=5)
ggsave(here("outputs/figures/XTvFlex_alphadiv.pdf"), dpi=300, w=5, h=5)
##### MAASLIN2 #####
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("Maaslin2", version = "3.12", force = TRUE)
library(Maaslin2)

inputdata <- t(read.table(here("XTvFlex/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE))
metadata <- read.table(here("XTvFlex/groups_and_conditions.tsv"), sep="\t", header=TRUE)
metadata <- mutate(metadata, Sample=gsub('-', '.', Sample))
rownames(metadata) <- metadata$Sample
metadata <- filter(metadata, Replicate!="B")
metadata <- filter(metadata, Replicate!="C")
metadata <- mutate(metadata, Sample=gsub('NFlex', 'NF', Sample))


"""
Prevalence: set to 0.2 so 4/20 samples must have taxon
"""
fit_data <- Maaslin2(
  inputdata, metadata, here('XTvFlex/04_maaslin2'), 
  transform = "LOG",
  analysis_method = "LM",
  fixed_effects = c('Method'),
  random_effects = c('Donor'),
  normalization = 'TSS',
  standardize = FALSE,
  min_abundance = 0.05,
  min_prevalence = 0.2,
  plot_scatter = FALSE, 
  plot_heatmap = FALSE)


results <- fit_data$results %>% filter(qval < 0.05) %>% arrange(coef)
results <- mutate(results, DirectionEnriched = ifelse(coef < 0, "Illumina DNA Prep", "Nextera XT"))
results <- mutate(results, feature=gsub(" miscellaneous", "", gsub("\\.", " ", feature)))

ggplot(results, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

ggsave(here("outputs/figures/XTvFlex_genus_maaslin2.pdf"), dpi=300, w=5, h=6)
ggsave(here("outputs/figures/XTvFlex_genus_maaslin2.jpeg"), dpi=300, w=5, h=6)

