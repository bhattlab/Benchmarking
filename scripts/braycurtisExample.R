library(tidyverse)
library(vegan)
library(here)
library(reshape2)
library(ggsignif)
library(ggpubr)


# Read in count data
dist_matrix <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/braycurtis_matrices/braycurtis_distance_species.txt"), sep="\t", header = TRUE)

# convert row names to a main column
dist_matrix <- data.frame(names = row.names(dist_matrix), dist_matrix)
rownames(dist_matrix) <- NULL

# filter rows and columns to exclude controls and contaminated sample
dist_matrix <- dist_matrix %>% select(!starts_with("NCO") & !starts_with("PCO") & !matches("D03.ZH.R2")) # filter columns
dist_matrix <- mutate(dist_matrix, names = gsub('-', '.', names))
dist_matrix <- dist_matrix %>% filter(!stringr::str_detect(names, 'NCO|PCO') & !stringr::str_detect(names, 'D03.ZH.R2'))

# convert data to long format
data_long <- melt(data = dist_matrix, id.vars = "names", variable.name="comparison", value.name = "bc")
names(data_long) <- c("Sample1", "Sample2", "BrayCurtis")

data_long <- data_long %>% separate(Sample1, c("Donor1", "Condition1", "Replicate1"), remove=FALSE)
data_long <- data_long %>% separate(Sample2, c("Donor2", "Condition2", "Replicate2"), remove=FALSE)

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
distance_averages_within_donor <- mutate(distance_averages_within_donor, condition_codes = ifelse(as.character(Condition1) > as.character(Condition2), paste(Condition2, Condition1), paste(Condition1, Condition2)))


condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("NF OH", "NF OR", "NF OF", "NF", "NF ZF", "NF ZR", "NF ZH")

distance_averages_within_donor <- distance_averages_within_donor %>% mutate(condition_codes = fct_relevel(condition_codes, "NF OF", "NF OR", "NF OH", "NF ZF", "NF ZR", "NF ZH"))


ggplot(distance_averages_within_donor %>% filter(condition_codes %in% c("NF OF", "NF OR", "NF OH", "NF ZF", "NF ZR", "NF ZH")), aes(x=condition_codes, y=BrayCurtis, fill = condition_codes)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.8),
              alpha = 0.85,  color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_bw() + 
  xlab("Condition") + 
  ylab("Species-Level Bray Curtis Dissimilarity") + 
  scale_fill_manual(values = condition_palette) + 
  theme(
    legend.position = "none"
  ) + 
  stat_compare_means(method = "wilcox.test", paired=TRUE, comparisons=list(c("NF OF", "NF OR"), c("NF OR", "NF OH"), c("NF OF", "NF OH"), c("NF ZF", "NF ZR"), c("NF ZR", "NF ZH"), c("NF ZF", "NF ZH")), 
                     tip.length = 0, label = "p.signif")
  

