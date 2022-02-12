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
library(Maaslin2)

#####Maaslin2 for differential abundance analysis#####
#Installation of Maaslin2
#if(!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install(version = "3.12") 
#BiocManager::install("Maaslin2", version="3.12")

#Other installation method
#install.packages("devtools")
# library(devtools)
#devtools::install_github("biobakery/maaslin2")


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
kraken_genus_og <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)

kraken_genus <- t(read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE))


##Maaslin2 analysis of ZF v ZH
#Filter metadata to only include Zymo kit
benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample
zymoresults <- filter(benchmark_groups, Kit=="Zymo")

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, zymoresults, here('DNA/3.maaslin'), transform = "NONE",
  fixed_effects = c('Temperature'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Temperature,Frozen')

#Filter for significant and large effect sizes
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
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
benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample
omniresults <- filter(benchmark_groups, Kit=="Omnigene")

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, omniresults, here('DNA/3.maaslin'), transform = "NONE",
  fixed_effects = c('Temperature'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Temperature,Frozen')

#Filter for significant and large effect sizes
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
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
benchmark_groups <- read.table(here("outputs/tables/DNA_metadata.tsv"), sep="\t", header=TRUE)
rownames(benchmark_groups) <- benchmark_groups$sample
frozenresults <- filter(benchmark_groups, Temperature=="Frozen")

#Running Maaslin2
fit_data <- Maaslin2(
  kraken_genus, frozenresults, here('DNA/3.maaslin'), transform = "NONE",
  fixed_effects = c('Kit'),
  random_effects = c('Donor'),
  normalization = 'NONE',
  standardize = FALSE,
  min_abundance = 0.01,
  min_prevalence = 0.1, 
  reference = 'Kit,No Preservative')

#Filter for significant and large effect sizes only for NF v ZF
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
nfzfresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Zymo"))

ggplot(nfzfresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

ggsave(here("outputs/figures/DNA_NFZF_maaslin2.pdf"), dpi=300, w=5, h=6)
ggsave(here("outputs/figures/DNA_NFZF_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)

#Filter for significant and large effect sizes only for NF v OF
results <- fit_data$results %>% filter(qval < 0.05) %>% filter(abs(coef)>1) %>% arrange(coef)
nfofresults <- mutate(results, DirectionEnriched = ifelse(coef < 0, "No Preservative", "Omnigene"))

ggplot(nfofresults, aes(x=reorder(feature, -coef), y=coef)) +
  geom_col(aes(fill=DirectionEnriched), alpha = 0.8) + 
  theme_bw() + 
  labs(y = "Effect Size", x = "Genus", fill = "") +
  coord_flip() + 
  #scale_fill_manual(values = flexxt_palette) + 
  theme(axis.title.y = element_blank(), legend.position = "top", legend.justification = "center")

ggsave(here("outputs/figures/DNA_NFOF_maaslin2.pdf"), dpi=300, w=5, h=6)
ggsave(here("outputs/figures/DNA_NFOF_maaslin2.pdf.jpeg"), dpi=300, w=5, h=6)


##Generating heatmap for Kit and Temp comparison
#Read in dataframes
zymoresults <- read.table(here("DNA/3.maaslin/zymo/all_results.tsv"), sep="\t", header=TRUE)
omniresults <- read.table(here("DNA/3.maaslin/omni/all_results.tsv"), sep="\t", header=TRUE)
frozenresults <- read.table(here("DNA/3.maaslin/frozen/all_results.tsv"), sep="\t", header=TRUE)

zymoresults <- filter(zymoresults, value=="Hot") %>% mutate(index=paste(metadata, "Zymo"))
omniresults <- filter(omniresults, value=="Hot") %>% mutate(index=paste(metadata, "Omni"))
frozenresults <- mutate(frozenresults, index=paste(metadata,value))

#Concatenate the dataframes
results <- rbind(zymoresults, frozenresults, omniresults)
write.table(results, here("DNA/3.maaslin/concat_results.tsv"), row.names = FALSE, sep="\t", quote=FALSE)

#Generating the heatmap
results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
results <- filter(results, !grepl("unclassified", feature)) %>% filter(!grepl(".environmental", feature))
sigresults <- filter(results, abs(coef)>0.25)
uniquetaxa <- unique(sigresults$feature)
results <- filter(results, feature %in% uniquetaxa)
results <- mutate(results, coef2=ifelse(coef>3, 3, ifelse(coef < -3, -3, coef)))
limit <- max(abs(results$coef2)) * c(-1, 1)



ggplot(results, aes(x=index, y=feature, fill=coef2)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(~metadata, scales="free_x") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 15)) + 
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

