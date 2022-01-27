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
library(PairedData)


##### HEATMAP #####
results <- read.table(here("test/maaslin_test.txt"), sep="\t", header=TRUE)
results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
limit <- max(abs(results$coef)) * c(-1, 1)


ggplot(results, aes(x=comparison, y=feature, fill=coef)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(~metadata, scales="free_x") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F, shape=8) + 
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
  

heatmap(results, Colv = comparison, Rowv = coef, scale="column")


results2 <- results %>% select(comparison, feature, coef)
as.matrix(results2)


##### TESTING AGGREGATION #####
donor <- c("D1", "D1", "D1", "D1", "D1", "D1", "D2", "D2", "D2", "D2", "D2", "D2", "D3", "D3","D3","D3","D3","D3")
replicate <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
condition <- c("NF", "NF", "NF", "OF", "OF", "OF","NF", "NF", "NF", "OF", "OF", "OF","NF", "NF", "NF", "OF", "OF", "OF")
diversity <- c(1,2,3,4,5,6,7,8,9,2,4,6,7,8,9,1,9,2)
testdf <- data.frame(donor, replicate, condition, diversity)

testdf %>%
  group_by(donor, condition) %>%
  summarise_at(vars(diversity), list(name = mean))


#####Shannon Diversity Plots#####
benchmark_s <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t")
benchmark_groups <- read.csv(here("data/DNAExtraction.tsv"), sep="\t", header=TRUE)
#Separate the SampleID name into Donor, Condition and Replicate columns; remove=FALSE keeps the SampleID column
benchmark_groups <- benchmark_groups %>% separate(SampleID, c("Donor", "Condition", "Replicate"), remove=FALSE)
benchmark_groups <- mutate(benchmark_groups, sample=gsub("_", ".", SampleID))
benchmark_groups <- filter(benchmark_groups, Donor!="NCO" & Donor!="PCO")
benchmark_groups <- mutate(benchmark_groups, Preservation=substr(Condition,1,1))
benchmark_groups <- mutate(benchmark_groups, Temperature=substr(Condition,2,2))


shannon_div_s <- diversity(t(benchmark_s), index = "shannon")
div <- data.frame("shannon_div" = shannon_div_s, "sample" = names(shannon_div_s))
div_meta <- merge(div, benchmark_groups, by = "sample")
div_meta <- filter(div_meta, sample != "D03.ZH.R2")

names(div_meta)
div_grouped_replicates <- div_meta %>% dplyr::select(sample, shannon_div, SampleID, Donor, Condition, Replicate, Preservation, Temperature) %>%
  group_by(Donor, Condition) %>%
  summarise_at(vars(shannon_div), list(shannon_div = mean))

div_grouped_replicates <- mutate(div_grouped_replicates, Preservation=substr(Condition,1,1))
div_grouped_replicates <- mutate(div_grouped_replicates, Temperature=substr(Condition,2,2))

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

div_meta <- div_meta %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))
div_grouped_replicates <- div_grouped_replicates %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))

# DIVERSITY PLOT BY CONDITION, all replicates
# ggplot(div_meta, aes(Condition, shannon_div)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85, aes(fill=Condition), color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
#   labs(x = "",
#        y = "Shannon Diversity", 
#        fill = "") + 
#   scale_color_manual(values = condition_palette) + 
#   scale_fill_manual(values = condition_palette) +
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(size = 12),
#         legend.position = "none",
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) 

# # DIVERSITY PLOT BY CONDITION, mean of replicates
# ggplot(div_grouped_replicates, aes(Condition, shannon_div)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85, aes(fill=Condition), color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
#   labs(x = "",
#        y = "Shannon Diversity", 
#        fill = "") + 
#   scale_color_manual(values = condition_palette) + 
#   scale_fill_manual(values = condition_palette) +
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(size = 12),
#         legend.position = "none",
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) 

# # DIVERSITY PLOT BY CONDITION, mean of replicates
# preserve_effect <- ggplot(div_grouped_replicates %>% filter(Condition %in% c("NF", "OF", "ZF")), aes(Condition, shannon_div)) + 
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#               alpha = 0.85, aes(fill=Condition), color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
#   labs(x = "",
#        y = "Shannon Diversity", 
#        fill = "") + 
#   scale_color_manual(values = condition_palette) + 
#   scale_fill_manual(values = condition_palette) +
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(size = 12),
#         legend.position = "none",
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) 
# preserve_effect


preserve_effect <- ggplot(div_grouped_replicates %>% filter(Condition %in% c("NF", "OF", "ZF")), aes(Condition, shannon_div)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85, aes(fill=Condition), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
  labs(x = "",
       y = "Shannon Diversity", 
       fill = "", 
       title = "Preservative Effect") + 
  ylim(2,4) +
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

temperature_effect <- ggplot(div_grouped_replicates , aes(Condition, shannon_div)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85, aes(fill=Condition), color = "darkgray") +
  geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) + 
  labs(x = "",
       y = "Shannon Diversity", 
       fill = "", 
       title = "Temperature Effect") + 
  ylim(2,4) +
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        title = element_text(size = 13),
        plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm")) +
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("OF", "OR"), c("OF", "OH"), c("ZF", "ZR"), c("ZF", "ZH")), 
                     tip.length = 0, label="p.signif")
temperature_effect

plot_grid(preserve_effect, temperature_effect, rel_widths = c(1,1.9))

ggsave(here("outputs/figures/DNA_alphadiv_paired_by_condition.jpg"), dpi=300, w=8, h=4)
ggsave(here("outputs/figures/DNA_alphadiv_paired_by_condition.pdf"), dpi=300, w=8, h=4)


# ## trying paired plots
# ggpaired(div_grouped_replicates %>% filter(Temperature == "F"), x="Condition", y="shannon_div",
#          fill="Condition", line.color="gray", id="Donor",
#          point.size = 1,
#          line.size = 1) +
#   stat_compare_means(paired=TRUE) +
#   #geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
#   #            alpha = 0.85, aes(fill=Condition), color = "darkgray") +
#   scale_color_manual(values = condition_palette) +
#   scale_fill_manual(values = condition_palette) +
#   #geom_boxplot(outlier.shape = NA, aes(fill = Condition), alpha = 0.7) +
#   labs(x = "",
#        y = "Shannon Diversity",
#        fill = "") +
#   theme(axis.text.x = element_text(size = 12),
#         legend.position = "none",
#         plot.margin = unit(c(0.2, 0.5, 0, 0.2), unit = "cm"))
