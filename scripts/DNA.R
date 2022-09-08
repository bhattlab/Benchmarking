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

# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#read in the QSU relative data as three separate dataframes
raw <- read.csv(here("QSU_Data/raw_data_relative.csv"), header=TRUE) # raw per sample information
model <- read.csv(here("QSU_Data/raw_model_relative.csv"), header=TRUE) # means and confidence intervals for all conditions
sig <- read.csv(here("QSU_Data/raw_significance_relative.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

#Filter QSU relative data to remove outdated absolute data
model <- filter(model, feature!="MicrobesPerGram" & feature!="mean_16s_count" & feature!="DNAConcentration")
raw <- raw %>% select(-mean_16s_count)
raw <-raw %>% select(-DNAConcentration)
raw <-raw %>% select(-Patient, -Sample_Type, -Replication)
sig <- filter(sig, feature!="mean_16s_count")
colnames(model)[colnames(model) == "estimate"] <- "prediction"

# read in the QSU absolute data as three separate dataframes 
rawabs <- read.csv(here("QSU_Data/raw_data_all.csv"), header=TRUE) # raw per sample information
modelabs <- read.csv(here("QSU_Data/raw_data_model_means.csv"), header=TRUE) # means and confidence intervals for all conditions
sigabs <- read.csv(here("QSU_Data/raw_data_significance.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

#merge the absolute and relative QSU data
raw <- merge(raw,rawabs, by="Sample")
model <- rbind(model, modelabs)
sig <- rbind(sig,sigabs)

# format and edit dataframes 
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values
model <- model %>% mutate(Sample_Type = Condition) # make all dataframes have a column called Sample_Type
raw <- raw %>% mutate(Preservative = substr(Sample_Type, 1, 1)) # make preservative column (first character of sample code)
raw <- raw %>% mutate(Temperature = substr(Sample_Type, 2, 2)) # make temperature column (second character of sample code)
raw <- raw %>% mutate(Sample_Type = fct_relevel(Sample_Type, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH")) # force the ordering of the conditions
raw <- raw %>% mutate(TemperatureLong = ifelse(Temperature == "F", "-80C", ifelse(Temperature == "R", "23C", "40C"))) # write out temperature
raw <- raw %>% mutate(PreservativeLong = ifelse(Preservative == "N", "None", ifelse(Preservative == "O", "Omnigene", "Zymo"))) # write out preservative
raw <- raw %>% mutate(Label = paste(TemperatureLong, PreservativeLong, sep="\n")) # make a label of temperature and preserative
raw <- raw %>% mutate(hiddenLabel=ifelse(Sample_Type == "OR", "OMNIgene", ifelse(Sample_Type== "ZR", "Zymo", ifelse(Sample_Type == "NF", "None", "")))) # make a label just for the "R" samples


##### FIGURE 2 #####
##### Microbes per Gram 
# in order, this code: 
# makes the main plot: 
# plots a vertical line to divide the none and Omni and the Omni and Zymo plots
# plots jittered points colored by condition
# colors the points based on the condition_palette
# labels the x axis with our specified condition_labels (the temperatures only)
# adds error bars based on the model dataframe, filtering to only show errors for the MicrobesPerGram rows, and based on the high and low CIs in that table
# adds a mean point based on the model dataframe, filtering to only show means for the MicrobesPerGram rows
# scales the y axis by log10
# plots the p-values based on the sig dataframe, filtering to only show p values for the MicrobesPerGram rows, labeling with asterisks based on the p.signif column, 
#       and placing the p-value at the specified y positions
#       EITHER plotting all p-values (commented out) OR plotting only significant p-values (make sure to edit y.position to only have same number of entries)
# makes the theme bw (a bit cleaner)
# changes the y-axis label
# adjusts the theme to remove the x axis title, remove the vertical gridlines, removes the legend, makes all the text size 12, and remove the margins around the plot
# ---
# makes a second plot: 
# takes the raw dataframe, filters to only look at donor 1, replicate 1, plots such that x is the condition, y is 0
# plots text for each condition where the text is the "hiddenLabel" (where OR and ZR map to Omnigene and Zymo, and OF/OH/ZF/ZH are blank - this makes the labels centered across the three temperatures)
# makes the theme void (no grid lines, axes, labels, etc)
# removes the plot margins
# ---
# combines the plots: 
# uses cowplot to plot both plots, in a single column, with the second plot having a very short height, and aligning the left axis of the plots (not centering relative to each other)

#Figure 2A: Absolute Count across Conditions
p <- ggplot(raw, aes(x = Sample_Type, y=MicrobesPerGram)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=1.5) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "MicrobesPerGram"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "MicrobesPerGram"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  scale_y_log10() +
  # stat_pvalue_manual(sig %>% filter(Feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(14.5, 14.3, 14, 13.8, 13, 12.8), 
  #                    tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sig %>% filter(feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(13.8, 13, 12.8), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Total Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

absolute <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='l')
absolute

#Figure2B: Changing Absolute count of Bacteroidetes and Firmicutes across conditions
abs<- raw[c("Sample_Type","Absolute.Abundance..Firmicutes", "Absolute.Abundance..Bacteroidetes")]
meltabs <- melt(abs, id.vars = "Sample_Type", variable.name = "Phyla", value.name = "MicrobesPerGram")
fbmodel <- filter(model, feature%in%c("Absolute Abundance: Firmicutes", "Absolute Abundance: Bacteroidetes"))
fbsig <- filter(sig, feature%in%c("Absolute Abundance: Firmicutes", "Absolute Abundance: Bacteroidetes"))

#Plot Order
fbmodel <- mutate(fbmodel, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                            ifelse(Sample_Type == "OF", 2, 
                                                   ifelse(Sample_Type == "OR", 3, 
                                                          ifelse(Sample_Type == "OH", 4, 
                                                                 ifelse(Sample_Type == "ZF", 5,
                                                                        ifelse(Sample_Type == "ZR", 6, 7)))))))

#Plot Firmicutes and Bacteroidetes Absolute Abundance by Condition
fbmodel$feature <- gsub("Absolute Abundance: Firmicutes","Firmicutes",fbmodel$feature)
fbmodel$feature <- gsub("Absolute Abundance: Bacteroidetes","Bacteroidetes",fbmodel$feature)

#Plot Firmicutes and Bacteroidetes NFvOF
nfof <- ggplot(meltabs %>% filter(Sample_Type=="NF"|Sample_Type=="OF"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labels) +
  geom_line(data=fbmodel%>% filter(feature=="Bacteroidetes")%>% filter(Sample_Type =="NF" | Sample_Type=="OF"), aes(y=prediction, group=1),linetype="dashed")+
  geom_line(data=fbmodel%>% filter(feature=="Firmicutes") %>% filter(Sample_Type =="NF" | Sample_Type=="OF"), aes(y=prediction, group=1))+
  geom_point(data=fbmodel %>% filter(Sample_Type=="NF" | Sample_Type=="OF"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits = c(2e11, 4e12)) +
  theme_bw() + 
  ylab("OMNIgene\n\nMicrobes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none",
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle("Preservative Effect")

nfof

#Plot Firmicutes and Bacteroidetes OFvORvOH
omni <- ggplot(meltabs %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labels) +
  geom_line(data=fbmodel%>% filter(feature=="Bacteroidetes")%>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(x=reorder(Sample_Type, PlotOrder), y=prediction, group=1), linetype="dashed")+
  geom_line(data=fbmodel%>% filter(feature=="Firmicutes") %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(x=reorder(Sample_Type, PlotOrder), y=prediction, group=1))+
  geom_point(data=fbmodel %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits = c(2e11, 4e12)) +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle("Temperature Effect")

omni

o <- plot_grid(nfof, omni, nrow=1, ncol=2, align="h", rel_widths = c(0.8,1))
o

#Plot Firmicutes and Bacteroidetes NFvZF
nfzf <- ggplot(meltabs %>% filter(Sample_Type=="NF"|Sample_Type=="ZF"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labels) +
  geom_line(data=fbmodel%>% filter(feature=="Bacteroidetes")%>% filter(Sample_Type =="NF" | Sample_Type=="ZF"), aes(y=prediction, group=1),linetype="dashed")+
  geom_line(data=fbmodel%>% filter(feature=="Firmicutes") %>% filter(Sample_Type =="NF" | Sample_Type=="ZF"), aes(y=prediction, group=1))+
  geom_point(data=fbmodel %>% filter(Sample_Type=="NF" | Sample_Type=="ZF"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits=c(5e10,9e11)) +
  theme_bw() + 
  ylab("Zymo\n\nMicrobes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none") 

nfzf

#Plot Firmicutes and Bacteroidetes ZFvZRvZH
zymo <- ggplot(meltabs %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_line(data=fbmodel%>% filter(feature=="Bacteroidetes")%>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(x=reorder(Sample_Type, PlotOrder), y=prediction, group=1), linetype="dashed")+
  geom_line(data=fbmodel%>% filter(feature=="Firmicutes") %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(x=reorder(Sample_Type, PlotOrder), y=prediction, group=1))+
  geom_point(data=fbmodel %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits=c(5e10,9e11)) +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

zymo

z <- plot_grid(nfzf, zymo, nrow=1, ncol=2, align="h", rel_widths = c(0.8,1))
z

oz <- plot_grid(o, z, nrow=2, ncol=1, align="v")
oz

#Make legend
line <- c("solid","dashed")
names(line) <- c("Firmicutes","Bacteroidetes")
l <- ggplot(fbmodel %>% filter(feature=="Firmicutes" | feature=="Bacteroidetes"), aes(x=prediction, y=CI_low, linetype=feature)) +
  geom_line()+
  scale_linetype_manual(values=line)+
  theme_bw()+
  theme(legend.title = element_blank(),legend.direction="horizontal", legend.position = "bottom")
l

le <- get_legend(l)

fig2b <- plot_grid(oz, le, nrow=2, ncol=1, rel_heights=c(1, 0.05))
fig2b

#ggsave(here("QSU_Data/Fig2Bline.jpeg"), dpi=300, h=5, w=8)


#Figure 2C: Firmcutes and Bacteroidetes Ratio
#Read in QSU ratio data as three dataframes
rratio <- read.csv(here("QSU_Data/rawratio.csv"), header=TRUE) # raw per sample information
mratio <- read.csv(here("QSU_Data/modelratio.csv"), header=TRUE) # means and confidence intervals for all conditions
sratio <-read.csv(here("QSU_Data/sigratio.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

rratio <- mutate(rratio, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                          ifelse(Sample_Type == "OF", 2, 
                                                 ifelse(Sample_Type == "OR", 3, 
                                                        ifelse(Sample_Type == "OH", 4, 
                                                               ifelse(Sample_Type == "ZF", 5,
                                                                      ifelse(Sample_Type == "ZR", 6, 7)))))))

r <- ggplot(rratio, aes(x = reorder(Sample_Type, PlotOrder), y=ratio.of.B.F)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  #geom_violin(aes(color=Sample_Type)) + 
  geom_jitter(width=0.2, shape=16, size=1.5, aes(color=Sample_Type))+
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=mratio %>% filter(feature == "Bacteroidetes/Firmicutes"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=mratio %>% filter(feature == "Bacteroidetes/Firmicutes"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  #stat_pvalue_manual(sig %>% filter(Feature == "Richness 0.01%"), y.position=c(127, 124, 118, 121, 118, 121), 
  #                   tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sratio %>% filter(feature == "Bacteroidetes/Firmicutes") %>% filter(p.adj <= 0.05), y.position=c(2.5, 2.4,2.3,2.2), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Bacteroidetes/Firmicutes \n") + 
  ylim(0,2.5) + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
r


#ggsave(here("QSU_Data/Figure2ratio.jpeg"), w=5, h=4, dpi=300)

#Plot Figure 2C together
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.05, 0.05) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b
r<-plot_grid(r,b, nrow = 2, ncol = 1,rel_heights=c(1,0.1), align="v", axis="l")
r

#ggsave(here("QSU_Data/Figure2Ceven.jpeg"), w=6.5, h=6.5, dpi=300)

#Plot Figure 2!
two<-plot_grid(absolute, fig2b, r, nrow=1, ncol=3, scale=0.9, labels=c("A","B", "C"), align = "v",axis="tb")
two
ggsave(here("QSU_Data/Figure2.jpeg"), dpi=300, h=8, w=22)
ggsave(here("QSU_Data/Figure2.pdf"), dpi=300, h=5.5, w=16)


##### FIGURE 3 #####
#Relative abundance analysis
#3A: Stacked bar plot for all donors
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

genus <- read.csv(here("DNA/2.kraken/kraken2_classification/processed_results/taxonomy_matrices_classified_only/bracken_genus_percentage.txt"), sep="\t", header=TRUE)
phylum_plottingorder <- read.csv(here("QSU_Data/genus_to_phylum.csv"), header=TRUE)
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

#Filter out controls
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
barplot_pal <- c(pinks[3], rev(oranges[2:5]), rev(blues[1:10]), "#939393")

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
condition_labels2 <- c("OMNI 40°C","OMNI 23°C","OMNI -80°C","No Preservative -80°C","Zymo -80°C","Zymo 23°C","Zymo 40°C")
names(condition_labels2) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

dummy <- ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder.x), y=rel_abundance)) +
  geom_bar(stat="identity", aes(fill=Condition)) + 
  scale_fill_manual(values=condition_palette, labels=condition_labels2) + 
  theme(legend.direction="horizontal", legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
dummy
b <- get_legend(dummy)
fig3a <- plot_grid(r, b, nrow=2, ncol=1, rel_heights = c(1, 0.08), rel_widths = c(1, 1))
fig3a

#ggsave(here("QSU_Data/Fig3ADNA_stackedbar.pdf"), dpi=300, h=6, w=15)
#ggsave(here("QSU_Data/Fig3ADNA_stackedbar.jpeg"), dpi=300, h=6, w=15)

#3B: Shannon entropy across conditions
#Change column name
colnames(raw)[colnames(raw) == "Shannon.Entropy"] <- "Shannon"

#Shannon entropy
se <- ggplot(raw , aes(Sample_Type, Shannon)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Shannon Entropy"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Shannon Entropy"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  ylim(1,4.5) +
  stat_pvalue_manual(sig %>% filter(feature == "Shannon Entropy") %>% filter(p.adj <= 0.05), y.position=c(4.1, 4.3, 3.7,3.9), 
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
shannon <- plot_grid(se,b, nrow=2, ncol=1, rel_heights=c(1,0.1), align="v", axis='lr')
shannon

#Figure 3E:Phylum forest plotting testing

raw <- mutate(raw, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                                        ifelse(Sample_Type == "OF", 2, 
                                                               ifelse(Sample_Type == "OR", 3, 
                                                                      ifelse(Sample_Type == "OH", 4, 
                                                                             ifelse(Sample_Type == "ZF", 5,
                                                                                    ifelse(Sample_Type == "ZR", 6, 7)))))))
model <- mutate(model, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                    ifelse(Sample_Type == "OF", 2, 
                                           ifelse(Sample_Type == "OR", 3, 
                                                  ifelse(Sample_Type == "OH", 4, 
                                                         ifelse(Sample_Type == "ZF", 5,
                                                                ifelse(Sample_Type == "ZR", 6, 7)))))))

bact <- ggplot(raw , aes(x=reorder(Sample_Type, PlotOrder), Relative.Abundance..Bacteroidetes*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Bacteroidetes"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Bacteroidetes"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Bacteroidetes") %>% filter(p.adj <= 0.05), y.position=c(90, 95, 70, 75), 
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
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Firmicutes"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Firmicutes"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Firmicutes") %>% filter(p.adj <= 0.05), y.position=c(90, 95, 70,75), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylim(0,100) + 
  ylab("Relative Abundance") + 
  xlab("Firmicutes") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())

act <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Actinobacteria*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Actinobacteria"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Actinobacteria"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Actinobacteria") %>% filter(p.adj <= 0.05), y.position=c(25, 29, 17, 21, 17), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Relative Abundance") + 
  xlab("Actinobacteria") + 
  ylim(0,100)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
act
vir <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Viruses*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Viruses"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Viruses"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Viruses") %>% filter(p.adj <= 0.05), y.position=c(0.45, 0.38, 0.41, 0.38, 0.41),
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Relative Abundance") + 
  xlab("Viruses") + 
  ylim(0,0.5)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())

vir

fungi <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Fungi*100)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Fungi"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Fungi"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=3) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Fungi") %>% filter(p.adj <= 0.05), y.position=c(0.15, 0.13),
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Relative Abundance") + 
  xlab("Fungi") + 
  ylim(0,0.5) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())

fungi

pdiff <- plot_grid(bact, firm, act, vir, fungi, nrow=1, ncol=5, rel_widths = c(1.3, 1, 1, 1.1, 1))
pdiff
ggsave(here("QSU_Data/PhylumEnrichmentDraft.pdf"), dpi=300, h=3.2, w=4)

#Figure 3D: Metagenomic Bray Curtis #
bc_raw <- read.csv(here("QSU_Data/raw_data_all_braycurtis.csv"), header=TRUE)
bc_raw <- bc_raw %>% mutate(UniqueComparison = ifelse(Protocol_1 > Protocol_2, paste(Patient, Protocol_1, Protocol_2, sep="_"), paste(Patient, Protocol_2, Protocol_1, sep="_")))
bc_raw <- bc_raw %>% select(Patient, bcdist, UniqueComparison)
bc_unique <- distinct(bc_raw, UniqueComparison, bcdist, .keep_all = TRUE)

bc_model <- read.csv(here("QSU_Data/raw_data_model_means_braycurtis.csv"), header=TRUE)
bc_model <- bc_model %>% separate(Protocols, c("group1", "group2"), remove=FALSE) %>% filter(group1 == "NF") %>% filter(group1 != group2)
bc_model <- bc_model %>% mutate(Condition=group2)

bc_sig <- read.csv(here("QSU_Data/raw_data_significance_braycurtis.csv"), header=TRUE)
bc_sig_toNF <- bc_sig %>% filter(Feature == "BrayCurtisBetweenConditions")
bc_sig_toNF <- bc_sig_toNF %>% mutate(group1=gsub("_NF", "", group1)) %>% mutate(group2=gsub("_NF", "", group2))

# across condition Bray Curtis
bc_across <- bc_unique %>% separate(UniqueComparison, c("Donor", "group1", "group2"), remove=FALSE)
bc_across <- bc_across %>% filter(group1 != group2) %>% filter(group1 == "NF" | group2 == "NF")
bc_across <- bc_across %>% mutate(modelID = paste(group2, group1, sep="_"))
bc_across <- bc_across %>% mutate(group1 = fct_relevel(group1,  "OF", "OR", "OH", "ZF", "ZR", "ZH"))

bc_plot <- ggplot(bc_across, aes(x=group1, y=bcdist)) + 
  geom_vline(aes(xintercept=3.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=group1),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=Mean), size=2.5) +
  ylim(0,1) +
  stat_pvalue_manual(bc_sig_toNF %>% filter(p.adj <= 0.05), y.position=c(1, .9, .7, .8), 
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
bc_full <- plot_grid(bc_plot,b, nrow=2, ncol=1, rel_heights=c(1,0.05), align="v", axis='l')
bc_full

#Figure3C: Relative Richness 0.01%
p <- ggplot(raw, aes(x = Sample_Type, y=Richness.0.01.)) +
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) +
  geom_jitter(width=0.2, aes(color=Sample_Type), shape=16, size=2) +
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  #stat_pvalue_manual(sig %>% filter(Feature == "Richness 0.01%"), y.position=c(127, 124, 118, 121, 118, 121),
  #                   tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sig %>% filter(feature == "Richness 0.01%") %>% filter(p.adj <= 0.05), y.position=c(127, 124, 121),
                     tip.length=0, label = "p.signif") +
  theme_bw() +
  ylab("Number of Genera\n>0.01% Relative Abundance") +
  ylim(50,130) +
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

#Plot Figure 3!
a <- plot_grid(fig3a,shannon,nrow=1,ncol=2,scale=0.9,rel_widths=c(1, 0.3),labels=c("a","b"), align = "v")
a
b <- plot_grid(richness,bc_full,pdiff,nrow=1,ncol=3, scale=0.9, rel_widths=c(0.7, 0.6,1.1),labels=c("c","d","e"), align="v", axis="t")
b
three<-plot_grid(a, b, nrow=2, ncol=1, rel_widths=c(1, 1), align="h", axis="lr")
three

ggsave(here("QSU_Data/Figure3.pdf"), dpi=300, h=9, w=16)
ggsave(here("QSU_Data/Figure3.jpeg"), dpi=300, h=9, w=16)

##### GENUS HEATMAP #####
genus_sig <- read.csv(here("QSU_Data/genus_significance.csv"), header=TRUE)
genus_sig <- genus_sig %>% mutate(PercentFormatted = as.numeric(gsub("%.*", "", Percent)))
genus_sig <- genus_sig %>% mutate(PFormatted = as.numeric(ifelse(p == "< 0.001", "0.001", p)))
genus_sig <- genus_sig %>% mutate(Condition = gsub("_.*", "", Comparison))
genus_sig <- genus_sig %>% mutate(y.position = ifelse(Condition == "OF", 2, 
                                                      ifelse(Condition=="ZF", 1, 
                                                             ifelse(Condition=="OR", 4, 
                                                                    ifelse(Condition == "OH", 3, 
                                                                           ifelse(Condition == "ZR", 2, 1))))))

raw_abundance <- read.csv(here("QSU_Data/genus_model_means.csv"), header=TRUE)
raw_abundance <- raw_abundance %>% filter(Condition == "NF") %>% select(feature, Mean)
names(raw_abundance) <- c("feature", "NFMean")
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: unclassified Clostridiales.*", "Firmicutes: unclassified Clostridiales (miscel...", feature))
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: unclassified Firmicutes sensu.*", "Firmicutes: unclassified Firmicutes sensu stri...", feature))
raw_abundance <- raw_abundance %>% mutate(Phylum = gsub(":.*", "", feature))
raw_abundance <- raw_abundance %>% 
  group_by(Phylum) %>% 
  arrange(Phylum, desc(NFMean))
raw_abundance$x.position <- as.numeric(row.names(raw_abundance))


genus_sig <- merge(genus_sig, raw_abundance, by="feature", all.x = TRUE)
genus_sig <- genus_sig %>% mutate(Phylum = gsub(":.*", "", feature))
genus_sig <- genus_sig %>% mutate(Phylum = ifelse(Phylum == "null", "Virus", Phylum))
genus_sig <- genus_sig %>% mutate(Plabel = ifelse(PFormatted <= 0.05, "*", ""))
genus_sig <- genus_sig %>% mutate(Genus = gsub(" sensu.*", "", gsub("unclassified ", "other ", gsub(".*: ", "", gsub(" \\(misc.*", "", gsub(": environmental.*", "", feature))))))

preservative <- ggplot(genus_sig %>% filter(Condition == "OF" | Condition == "ZF"), aes(x=x.position, y=y.position, fill=PercentFormatted)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1, 2), labels=c("Zymo -80°C", "OMNI -80°C")) + 
  labs(title = "Preservative Effect (Relative to no preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-81, 81)) +
  labs(fill="Percent\nEnrichment")+
  geom_point(data = genus_sig %>% filter(Condition == "OF" | Condition == "ZF") %>% filter(Plabel == "*"), size=1, color="white") + 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  geom_vline(xintercept=48.5) + 
  geom_vline(xintercept=49.5) + 
  geom_vline(xintercept=50.5) + 
  geom_vline(xintercept=53.5) + 
  theme(panel.border = element_rect(size=1))

preservative

z <- get_legend(preservative)

temperature <- ggplot(genus_sig %>% filter(Condition == "OH" | Condition == "ZH"), aes(x=x.position, y=y.position, fill=PercentFormatted)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0), breaks=genus_sig$x.position, labels=genus_sig$Genus) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1,2), labels=c("Zymo 40°C", "OMNI 40°C")) + 
  labs(title = "Temperature Effect (Relative to -80°C for each preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(),
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-81, 81)) +
  geom_point(data = genus_sig %>% filter(Condition == "OH") %>% filter(Plabel == "*"), size=1, color="white") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  geom_vline(xintercept=48.5) + 
  geom_vline(xintercept=49.5) + 
  geom_vline(xintercept=50.5) + 
  geom_vline(xintercept=53.5) + 
  theme(panel.border = element_rect(size=1))


temperature

wolegend <- plot_grid(preservative + theme(legend.position="none"), temperature + theme(legend.position="none"), nrow=2, ncol=1, align="v", axis="lr")
main <- plot_grid(wolegend, z, ncol=2, nrow=1, rel_widths=c(1, 0.3), rel_heights = c(1, 0.8))

raw_abund <- ggplot(genus_sig %>% filter(Condition == "ZF"), aes(x=x.position, y=y.position, fill=NFMean)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1), labels=c("Abundance")) + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("#F3F3F3", "black")) +
  labs(fill="Baseline\nAbundance")+ 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  geom_vline(xintercept=48.5) + 
  geom_vline(xintercept=49.5) + 
  geom_vline(xintercept=50.5) + 
  geom_vline(xintercept=53.5) + 
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
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  geom_vline(xintercept=48.5) + 
  geom_vline(xintercept=49.5) + 
  geom_vline(xintercept=50.5) + 
  geom_vline(xintercept=53.5) + 
  theme(panel.border = element_rect(size=1))
phylum_scale 

z3 <- get_legend(phylum_scale)


main_andabund <- plot_grid(phylum_scale + theme(legend.position="none"),
                           raw_abund + theme(legend.position="none"), 
                           preservative + theme(legend.position="none"), 
                           temperature + theme(legend.position="none"), 
          nrow=4, ncol=1, align="v", axis="lr", rel_heights = c(0.5, 0.5, 2.5, 2.5))

legend_plot <- plot_grid(z3,z2, z, ncol=3, nrow=1, align = "v")
legend_plot

plot_grid(main_andabund, legend_plot, nrow=1, ncol=2, rel_widths=c(1, 0.2), scale=0.9)
ggsave(here("QSU_Data/Figure3_heatmap_phylumscale.pdf"), dpi=300, w=14, h=5)
ggsave(here("QSU_Data/Figure3_heatmap_phylumscale.jpeg"), dpi=300, w=14, h=5)

##### heatmap with phylum labels #####
phylum_positions <- genus_sig %>% select(Phylum, x.position) %>%
  group_by(Phylum) %>%
  summarise_at(vars(x.position), list(Position = mean))

test2 <- ggplot(phylum_positions, aes(x=Position, y=0)) + 
  geom_text(aes(label=Phylum), size=3, hjust=1) + 
  theme_void() + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

plot_grid(test2, 
          raw_abund + theme(legend.position="none"), 
          preservative + theme(legend.position="none"), 
          temperature + theme(legend.position="none"), 
          nrow=4, ncol=1, align="v", axis="lr", rel_heights = c(1, 0.5, 2.5, 2.5))


##### heatmap without phylum panel
main_andabund <- plot_grid(raw_abund + theme(legend.position="none"), 
                           preservative + theme(legend.position="none"), 
                           temperature + theme(legend.position="none"), 
                           nrow=3, ncol=1, align="v", axis="lr", rel_heights = c(0.5, 2, 2.5))

main_andabund
legend_plot <- plot_grid(z2, z, ncol=2, nrow=1, align = "v")
legend_plot

plot_grid(main_andabund, legend_plot, nrow=1, ncol=2, rel_widths=c(1, 0.2), scale=0.9)
ggsave(here("QSU_Data/Figure3_heatmap.pdf"), dpi=300, w=12, h=4)
ggsave(here("QSU_Data/Figure3_heatmap.jpeg"), dpi=300, w=12, h=4)


##### heatmap with only top phyla
genus_top <- genus_sig %>% filter(Phylum %in% c("Actinobacteria", "Bacteroidetes", "Firmicutes"))

preservative <- ggplot(genus_top %>% filter(Condition == "OF" | Condition == "ZF"), aes(x=x.position, y=y.position, fill=PercentFormatted)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1, 2), labels=c("Zymo -80°C", "OMNI -80°C")) + 
  labs(title = "Preservative Effect (Relative to no preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-81, 81)) +
  labs(fill="Percent\nEnrichment")+
  geom_point(data = genus_top %>% filter(Condition == "OF" | Condition == "ZF") %>% filter(Plabel == "*"), size=1, color="white") + 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  theme(panel.border = element_rect(size=1))

preservative

z <- get_legend(preservative)

temperature <- ggplot(genus_top %>% filter(Condition != "OF" | Condition != "ZF"), aes(x=x.position, y=y.position, fill=PercentFormatted)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0), breaks=genus_sig$x.position, labels=genus_sig$Genus) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1,2,3,4), labels=c("Zymo 40°C", "Zymo 23°C", "OMNI 40°C", "OMNI 23°C")) + 
  labs(title = "Temperature Effect (Relative to -80°C for each preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(),
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-81, 81)) +
  geom_point(data = genus_top %>% filter(Condition != "ZF" & Condition != "OF") %>% filter(Plabel == "*"), size=1, color="white") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  geom_hline(yintercept=2.5) +
  theme(panel.border = element_rect(size=1))


temperature

wolegend <- plot_grid(preservative + theme(legend.position="none"), temperature + theme(legend.position="none"), nrow=2, ncol=1, align="v", axis="lr")
main <- plot_grid(wolegend, z, ncol=2, nrow=1, rel_widths=c(1, 0.3), rel_heights = c(1, 0.8))

raw_abund <- ggplot(genus_top %>% filter(Condition == "ZF"), aes(x=x.position, y=y.position, fill=NFMean)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1), labels=c("Abundance")) + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("#F3F3F3", "black")) +
  labs(fill="Baseline\nAbundance")+ 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  theme(panel.border = element_rect(size=1))
raw_abund

z2 <- get_legend(raw_abund)

phylum_scale <- ggplot(genus_top %>% filter(Condition == "ZF"), aes(x=x.position, y=y.position, fill=Phylum)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1), labels=c("Phylum")) + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(fill="Phylum") + 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=11.5) + 
  theme(panel.border = element_rect(size=1))
phylum_scale 

z3 <- get_legend(phylum_scale)


main_andabund <- plot_grid(phylum_scale + theme(legend.position="none"),
                           raw_abund + theme(legend.position="none"), 
                           preservative + theme(legend.position="none"), 
                           temperature + theme(legend.position="none"), 
                           nrow=4, ncol=1, align="v", axis="lr", rel_heights = c(0.5, 0.5, 2.5, 2.5))

legend_plot <- plot_grid(z3,z2, z, ncol=3, nrow=1, align = "v")
legend_plot

plot_grid(main_andabund, legend_plot, nrow=1, ncol=2, rel_widths=c(1, 0.2), scale=0.9)
ggsave(here("QSU_Data/Figure3_heatmap_phylumscale_withRT.pdf"), dpi=300, w=14, h=6)
ggsave(here("QSU_Data/Figure3_heatmap_phylumscale_withRT.jpeg"), dpi=300, w=14, h=6)

##### heatmap with phylum labels #####
phylum_positions <- genus_top %>% select(Phylum, x.position) %>%
  group_by(Phylum) %>%
  summarise_at(vars(x.position), list(Position = mean))

test2 <- ggplot(phylum_positions, aes(x=Position, y=0)) + 
  xlim(c(0,48)) +
  geom_text(aes(label=Phylum), size=3.5, hjust=1.3) + 
  theme_void() + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
test2

test3 <- plot_grid(test2, 
          raw_abund + theme(legend.position="none"), 
          preservative + theme(legend.position="none"), 
          temperature + theme(legend.position="none"), 
          nrow=4, ncol=1, align="v", axis="lr", rel_heights = c(1, 0.5, 2.5, 2.5))

legend_plot <- plot_grid(z2, z, ncol=2, nrow=1, align = "v")
legend_plot

plot_grid(test3, legend_plot, nrow=1, ncol=2, rel_widths=c(1, 0.2), scale=0.9)
ggsave(here("QSU_Data/Figure3_heatmap_phylumlabels.pdf"), dpi=300, w=17, h=6)
ggsave(here("QSU_Data/Figure3_heatmap_phylumlabels.jpeg"), dpi=300, w=17, h=6)


##### SUPPLEMENT
# DNA BC
bc_raw <- read.csv(here("QSU_Data/raw_data_all_braycurtis.csv"), header=TRUE)
bc_raw <- bc_raw %>% mutate(UniqueComparison = ifelse(Protocol_1 > Protocol_2, paste(Patient, Protocol_1, Protocol_2, sep="_"), paste(Patient, Protocol_2, Protocol_1, sep="_")))
bc_raw <- bc_raw %>% select(Patient, bcdist, UniqueComparison)
bc_unique <- distinct(bc_raw, UniqueComparison, bcdist, .keep_all = TRUE)
bc_unique <- bc_unique %>% separate(UniqueComparison, into=c("Donor", "Group1", "Group2"), sep="_", remove=FALSE)
bc_unique <- bc_unique %>% filter(Group1 == Group2)
bc_unique <- bc_unique %>% mutate(Condition = Group1)


bc_model <- read.csv(here("QSU_Data/raw_data_model_means_braycurtis.csv"), header=TRUE)
bc_model <- bc_model %>% separate(Protocols, c("group1", "group2"), remove=FALSE) 
bc_model <- bc_model %>% filter(group1 == group2)
bc_model <- bc_model %>% mutate(Condition=group2)

bc_sig <- read.csv(here("QSU_Data/raw_data_significance_braycurtis.csv"), header=TRUE)
bc_sig <- bc_sig %>% filter(Feature == "BrayCurtisWithinCondition")
bc_sig <- bc_sig %>% mutate(y.position = 0.2)
bc_unique <- bc_unique %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH"))


bc_plot <- ggplot(bc_unique, aes(x=Condition, y=bcdist)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) + 
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Condition),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=Mean), size=2.5) +
  ylim(0,1) +
  stat_pvalue_manual(bc_sig %>% filter(p.adj <= 0.05), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Bray-Curtis Dissimilarity") + 
  ggtitle("Metagenomic") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
bc_plot

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1") , aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") +
  ylim(-0.05, 0.05) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"))
b
bc_full_dna <- plot_grid(bc_plot,b, nrow=2, ncol=1, rel_heights=c(1,0.05), align="v", axis='l')
bc_full_dna

# RNA BC
bc_raw <- read.csv(here("QSU_Data/raw_data_rna_all_braycurtis.csv"), header=TRUE)
bc_raw <- bc_raw %>% mutate(UniqueComparison = ifelse(Protocol_1 > Protocol_2, paste(Patient, Protocol_1, Protocol_2, sep="_"), paste(Patient, Protocol_2, Protocol_1, sep="_")))
bc_raw <- bc_raw %>% select(Patient, bcdist, UniqueComparison)
bc_unique <- distinct(bc_raw, UniqueComparison, bcdist, .keep_all = TRUE)
bc_unique <- bc_unique %>% separate(UniqueComparison, into=c("Donor", "Group1", "Group2"), sep="_", remove=FALSE)
bc_unique <- bc_unique %>% filter(Group1 == Group2)
bc_unique <- bc_unique %>% mutate(Condition = Group1)


bc_model <- read.csv(here("QSU_Data/raw_data_rna_model_means_braycurtis.csv"), header=TRUE)
bc_model <- bc_model %>% separate(Protocols, c("group1", "group2"), remove=FALSE) 
bc_model <- bc_model %>% filter(group1 == group2)
bc_model <- bc_model %>% mutate(Condition=group2)

bc_sig <- read.csv(here("QSU_Data/raw_data_rna_significance_braycurtis.csv"), header=TRUE)
bc_sig <- bc_sig %>% filter(feature == "BrayCurtisWithinCondition")
bc_sig <- bc_sig %>% mutate(y.position = 0.2)
bc_unique <- bc_unique %>% mutate(Condition = fct_relevel(Condition, "NF", "OF", "OR", "ZF", "ZR", "ZH"))

bc_plot <- ggplot(bc_unique, aes(x=Condition, y=bcdist)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) + 
  geom_vline(aes(xintercept=3.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Condition),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=estimate), size=2.5) +
  ylim(0,1) +
  stat_pvalue_manual(bc_sig %>% filter(p.adj <= 0.05), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Bray-Curtis Dissimilarity") + 
  ggtitle("Metatranscriptomic") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
bc_plot

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1") , aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") +
  ylim(-0.05, 0.05) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"))
b
bc_full_rna <- plot_grid(bc_plot,b, nrow=2, ncol=1, rel_heights=c(1,0.05), align="v", axis='l')
bc_full_rna


plot_grid(bc_full_dna, bc_full_rna, nrow=1, ncol=2, scale=0.9, labels=c("a", "b"))
ggsave(here("QSU_Data/Supplement_BC.jpg"), dpi=300, w=12, h=5)
ggsave(here("QSU_Data/Supplement_BC.pdf"), dpi=300, w=12, h=5)

