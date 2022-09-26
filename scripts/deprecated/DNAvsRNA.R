library(dplyr)
library(ggplot2)

##### COMPARE DNA AND RNA GENUS LEVEL ABUNDANCES #####

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

# Read in DNA genus
dna_genus <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE)
rna_genus <- read.csv(here("RNA/02_kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_phylum_percentage.txt"), sep="\t", header=TRUE)

# filter DNA abundance data to above abundance threshold, convert to long form, reformat sample names
abundance_threshold <- 0
#abundance_threshold <- sort(rowSums(dna_genus), decreasing = T)[n_taxa]
dna_bracken_plot <- dna_genus[rowSums(dna_genus) >= abundance_threshold,]
dna_bracken_plot <- rbind(dna_bracken_plot, t(data.frame("Other" =  100 - colSums(dna_bracken_plot))))
dna_bracken_plot$Genus <- row.names(dna_bracken_plot)
dna_bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", dna_bracken_plot$Genus)
dna_bracken_long <- melt(dna_bracken_plot, id.vars = "Genus", variable.name = "Sample", value.name = "rel_abundance")
dna_bracken_long <- mutate(dna_bracken_long, Sample=gsub("\\.", "_", Sample))

# filter RNA abundance data to above abundance threshold, convert to long form, reformat sample names
abundance_threshold <- 0
#abundance_threshold <- sort(rowSums(dna_genus), decreasing = T)[n_taxa]
rna_bracken_plot <- rna_genus[rowSums(rna_genus) >= abundance_threshold,]
rna_bracken_plot <- rbind(rna_bracken_plot, t(data.frame("Other" =  100 - colSums(rna_bracken_plot))))
rna_bracken_plot$Genus <- row.names(rna_bracken_plot)
rna_bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", rna_bracken_plot$Genus)
rna_bracken_long <- melt(rna_bracken_plot, id.vars = "Genus", variable.name = "Sample", value.name = "rel_abundance")
rna_bracken_long <- mutate(rna_bracken_long, Sample=gsub("\\.", "_", Sample))

# merge DNA and RNA databases using a full join (include taxa that are unique to either dataset)
bracken_long <- merge(x=dna_bracken_long, y=rna_bracken_long, all=TRUE, by=c("Sample", "Genus"))
names(bracken_long) <- c("Sample", "Genus", "DNA_Abund", "RNA_abund")
bracken_long[is.na(bracken_long)] <- 0  # replace NAs with zero
## I don't think this bit is necessary - turning RNA/DNA columns into long form? 
# bracken_long <- melt(bracken_long, id=c("Genus", "Sample"))
# names(bracken_long) <- c("Sample", "Genus", "SequenceType", "Abundance")


# Merge in the metadata
colnames(metadata)[3]<-"Sample"
metadata <- metadata %>% select(Sample, Donor, Condition, Replicate)
bracken_pheno <- merge(bracken_long, metadata, by = "Sample")
 

condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

bracken_pheno <- bracken_pheno %>% filter(Donor != "NCO" & Donor != "PCO")

ggplot(bracken_pheno %>% filter(Condition != "OH"), aes(x=DNA_Abund, y=RNA_abund, color=Condition)) + 
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  geom_point(alpha=0.8) + 
  scale_color_manual(values = condition_palette) + 
  theme_bw() 

# group to only one point per donor/condition (merge tech reps)
bracken_pheno_grouped_reps <- bracken_pheno %>% group_by(Donor, Condition, Genus) %>%
  summarise_at(vars(DNA_Abund, RNA_abund), mean)

bracken_pheno_grouped_reps <- bracken_pheno_grouped_reps %>% mutate(Sample = paste(Donor, Condition, sep=""))

ggplot(bracken_pheno_grouped_reps %>% filter(Condition != "OH"), aes(x=DNA_Abund, y=RNA_abund, color=Condition)) + 
  geom_abline(intercept = 0, slope = 1, color = "lightgrey") +
  geom_point(alpha=0.8) + 
  scale_color_manual(values = condition_palette) + 
  theme_bw() 


# rough code to replace 0s with pseudocounts and calculate log fold change
bracken_pheno_grouped_reps <- mutate(bracken_pheno_grouped_reps, RNA_abund_pseudo = ifelse(RNA_abund > 0, RNA_abund, 0.001))
bracken_pheno_grouped_reps <- mutate(bracken_pheno_grouped_reps, DNA_abund_pseudo = ifelse(DNA_Abund > 0, DNA_Abund, 0.001))
bracken_pheno_grouped_reps <- mutate(bracken_pheno_grouped_reps, FoldChange = log2(RNA_abund_pseudo)-log2(DNA_abund_pseudo))

# group to only one point per donor/condition (merge tech reps)
bracken_pheno_grouped_donors <- bracken_pheno_grouped_reps %>% group_by(Condition, Genus) %>%
  summarise_at(vars(FoldChange), mean) %>% filter(Condition != "OH")

ggplot(bracken_pheno_grouped_donors, aes(x=reorder(Genus, FoldChange, FUN = median), y=FoldChange)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  xlab("Phylum") + 
  ylab("Log2 Fold Change (RNA/DNA)") + 
  theme(axis.text.x = element_text(angle = 90))



ggplot(bracken_pheno_grouped_reps, aes(x=reorder(Genus, FoldChange, FUN = median), y=FoldChange)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  xlab("Phylum") + 
  ylab("Log2 Fold Change (RNA/DNA)") + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylim(-10, 10) 

ggsave(here("outputs/figures/DNA_RNA_FoldChange_Boxplot.pdf"), dpi=300, w = 10, h=5)


#write.table(bracken_pheno, here("outputs/tables/DNA_RNA_abundance_phyla.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
