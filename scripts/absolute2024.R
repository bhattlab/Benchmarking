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
library(ggtext)
library(stringr)
library(ggpmisc)
library(ggbeeswarm)

# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# #read in the QSU relative data as three separate dataframes
# raw <- read.csv(here("QSU_Data/raw_data_relative.csv"), header=TRUE) # raw per sample information
# model <- read.csv(here("QSU_Data/raw_model_relative.csv"), header=TRUE) # means and confidence intervals for all conditions
# sig <- read.csv(here("QSU_Data/raw_significance_relative.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests
# 
# #Filter QSU relative data to remove outdated absolute data
# model <- filter(model, feature!="MicrobesPerGram" & feature!="mean_16s_count" & feature!="DNAConcentration")
# raw <- raw %>% select(-mean_16s_count)
# raw <-raw %>% select(-DNAConcentration)
# raw <-raw %>% select(-Patient, -Sample_Type, -Replication)
# sig <- filter(sig, feature!="mean_16s_count")
# colnames(model)[colnames(model) == "estimate"] <- "prediction"
# 
# # read in the QSU absolute data as three separate dataframes 
# rawabs <- read.csv(here("QSU_Data/raw_data_all.csv"), header=TRUE) # raw per sample information
# modelabs <- read.csv(here("QSU_Data/raw_data_model_means.csv"), header=TRUE) # means and confidence intervals for all conditions
# sigabs <- read.csv(here("QSU_Data/raw_data_significance.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests
# 
# #merge the absolute and relative QSU data
# raw <- merge(raw,rawabs, by="Sample")
# model <- rbind(model, modelabs)
# sig <- rbind(sig,sigabs)

#2024 absolute abundance data with corrected standard curves
raw <- read.csv(here("QSU_Data/2024/raw_data_dna_2024.csv"), header=TRUE) # raw per sample information
model <- read.csv(here("QSU_Data/2024/model_data_dna_2024.csv"), header=TRUE) # means and confidence intervals for all conditions
sig <- read.csv(here("QSU_Data/2024/sig_data_dna_2024.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

# format and edit dataframes 
sig <- sig[, -c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15)]
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig$p.signif <- str_extract(sig$percentchange, "\\]\\s*(.*)$") # Extract psignif after ']'
sig$p.signif <- str_replace(sig$p.signif, "\\]\\s*", "") # Remove the leading whitespace and closing bracket
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values
model <- model %>% mutate(Condition = Sample_Type) # make all dataframes have a column called Sample_Type
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
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=1) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels, guide = guide_axis(angle = 45)) +
  geom_errorbar(data=model %>% filter(feature == "MicrobesPerGram"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0, size=1) +
  geom_point(data=model %>% filter(feature == "MicrobesPerGram"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2) +
  scale_y_log10() +
  stat_pvalue_manual(sig %>% filter(feature == "MicrobesPerGram") %>% filter(p <= 0.05), y.position=c(13.8, 13.4, 12.8, 12.2), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Total Prokaryotes per Dry Gram") + 
  ggtitle(" ") +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.title.y = element_text(size=10), plot.title=element_text(size=8.5))
p
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3) + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

absolute <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.1 ), align="v", axis='l')
absolute
#ggsave(here("QSU_Data/deprecated/2024fig5a.jpeg"), dpi=300, h=5, w=5)

#Figure2C: Changing Absolute count of Bacteroidetes and Firmicutes across conditions
abs<- raw[c("Sample_Type","Absolute.Abundance..Firmicutes", "Absolute.Abundance..Bacteroidetes")]
meltabs <- melt(abs, id.vars = "Sample_Type", variable.name = "Phyla", value.name = "MicrobesPerGram")
fbmodel <- filter(model, feature%in%c("Absolute Abundance: Firmicutes", "Absolute Abundance: Bacteroidetes"))
fbsig <- filter(sig, feature%in%c("Absolute Abundance: Firmicutes", "Absolute Abundance: Bacteroidetes"))

condition_labelss <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labelss) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

condition_labels2 <- c("OMNIgene 40°C","OMNIgene 23°C","OMNIgene -80°C","No Preservative -80°C","Zymo -80°C","Zymo 23°C","Zymo 40°C")
names(condition_labels2) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plot Order
fbmodel <- mutate(fbmodel, PhyCondition=gsub("Absolute Abundance: ", "", paste(Condition,feature,sep="")))
fbmodel <- mutate(fbmodel, PlotOrder=ifelse(PhyCondition == "NFFirmicutes", 1, 
                                            ifelse(PhyCondition == "NFBacteroidetes", 2, 
                                                   ifelse(PhyCondition == "OFFirmicutes", 3, 
                                                          ifelse(PhyCondition == "OFBacteroidetes", 4, 
                                                                 ifelse(PhyCondition == "ORFirmicutes", 5,
                                                                        ifelse(PhyCondition == "ORBacteroidetes", 6,
                                                                               ifelse(PhyCondition == "OHFirmicutes", 7,
                                                                                      ifelse(PhyCondition == "OHBacteroidetes", 8,
                                                                                             ifelse(PhyCondition == "ZFFirmicutes", 9,
                                                                                                    ifelse(PhyCondition == "ZFBacteroidetes", 10,
                                                                                                           ifelse(PhyCondition == "ZRFirmicutes", 11,
                                                                                                                  ifelse(PhyCondition == "ZRBacteroidetes", 12,
                                                                                                                         ifelse(PhyCondition == "ZHFirmcutes", 13,
                                                                                                                                ifelse(PhyCondition == "ZHBacteroidetes", 14,15)))))))))))))))

#Plot Firmicutes and Bacteroidetes Absolute Abundance by Condition
fbmodel$feature <- gsub("Absolute Abundance: Firmicutes","Firmicutes",fbmodel$feature)
fbmodel$feature <- gsub("Absolute Abundance: Bacteroidetes","Bacteroidetes",fbmodel$feature)
fbmodel <- mutate(fbmodel, PlotPhyOrder = ifelse(feature=="Firmicutes", PlotOrder*10,PlotOrder))

### FIXME: need to figure out how to get y.position for significance bars in range
#Plot Firmicutes and Bacteroidetes NFvOF
nfof <- ggplot(meltabs %>% filter(Sample_Type=="NF"|Sample_Type=="OF"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labels) +
  geom_line(data=fbmodel%>% filter(Sample_Type =="NF" | Sample_Type=="OF"), aes(y=prediction, group=feature, linetype=feature),position=position_dodge(0.2), color="grey", size=0.5)+
  geom_errorbar(data=fbmodel %>% filter(Sample_Type=="NF" | Sample_Type=="OF"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low, ymax=CI_high, group=feature),position=position_dodge(0.2), size=0.5, width=0) +
  geom_point(data=fbmodel %>% filter(Sample_Type=="NF" | Sample_Type=="OF"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction, group=feature, color=Sample_Type), size=1, show.legend = FALSE, position=position_dodge(width=0.2)) +
  scale_y_log10(limits = c(1.3e10, 7e11)) +
  theme_bw() + 
  ylab("**OMNIgene**<br/>Prokaryotes/Dry Gram") + 
  ggtitle("Preservative") +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none",
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title=element_text(size=8.5),
        axis.title.y = element_markdown(size=9.5)) 
# stat_pvalue_manual(sig %>% filter(feature == "Absolute Abundance: Firmicutes") %>% filter(p.adj <= 0.05) %>% filter(group1 == "NF") %>% filter(group2 == "OF"), 
#                    tip.length=0, label = "p.signif") +
# stat_pvalue_manual(sig %>% filter(feature == "Absolute Abundance: Bacteroidetes") %>% filter(p.adj <= 0.05) %>% filter(group1 == "NF") %>% filter(group2 == "OF"), 
#                    label = "p.signif", tip.length =0, remove.bracket = TRUE, y.position=c(200000000000), position = position_dodge(width = 0.2)) 


nfof

#Plot Firmicutes and Bacteroidetes OFvORvOH
omni <- ggplot(meltabs %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(x=reorder(Sample_Type, PlotPhyOrder), y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labels) +
  geom_line(data=fbmodel%>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(y=prediction, group=feature, linetype=feature),position=position_dodge(0.2), color="grey", size=0.5)+
  geom_errorbar(data=fbmodel %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotPhyOrder), ymin=CI_low, ymax=CI_high, group=feature),position=position_dodge(0.2), size=0.5, width=0) +
  geom_point(data=fbmodel %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotPhyOrder), y=prediction, group=feature, color=Sample_Type), size=1, show.legend = FALSE, position=position_dodge(width=0.2)) +
  scale_y_log10(limits = c(1.3e10, 7e11)) +
  theme_bw() + 
  ylab("Prokaryotes/Dry Gram") + 
  ggtitle("Temperature") +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.title=element_text(size=8.5)) 


omni

o <- plot_grid(nfof, omni, nrow=1, ncol=2, align="h", rel_widths = c(1.2,1))
o

#Plot Firmicutes and Bacteroidetes NFvZF
nfzf <- ggplot(meltabs %>% filter(Sample_Type=="NF"|Sample_Type=="ZF"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labelss, guide = guide_axis(angle = 45)) +
  geom_line(data=fbmodel%>% filter(Sample_Type =="NF" | Sample_Type=="ZF"), aes(y=prediction, group=feature, linetype=feature),position=position_dodge(0.2), color="grey", size=0.5)+
  geom_errorbar(data=fbmodel %>% filter(Sample_Type=="NF" | Sample_Type=="ZF"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low, ymax=CI_high, group=feature),position=position_dodge(0.2), size=0.5, width=0) +
  geom_point(data=fbmodel %>% filter(Sample_Type=="NF" | Sample_Type=="ZF"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction, group=feature, color=Sample_Type), size=1, show.legend = FALSE, position=position_dodge(width=0.2)) +
  scale_y_log10(limits=c(3e09,3e12)) +
  theme_bw() + 
  ylab("**Zymo**<br/>Prokaryotes/Dry Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none",
        axis.title.y = element_markdown(size=9.5)) 

nfzf

#Plot Firmicutes and Bacteroidetes ZFvZRvZH
zymo <- ggplot(meltabs %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(x=reorder(Sample_Type, PlotPhyOrder), y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labelss, guide = guide_axis(angle = 45)) +
  geom_line(data=fbmodel%>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(y=prediction, group=feature, linetype=feature),position=position_dodge(0.2), color="grey", size=0.5)+
  geom_errorbar(data=fbmodel %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high, group=feature),position=position_dodge(0.2), size=0.5, width=0) +
  geom_point(data=fbmodel %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction, group=feature, color=Sample_Type), size=1, show.legend = FALSE, position=position_dodge(width=0.2)) +
  scale_y_log10(limits=c(3e09,3e12)) +
  theme_bw() + 
  ylab("Microbes/Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

zymo

z <- plot_grid(nfzf, zymo, nrow=1, ncol=2, align="h", rel_widths = c(1.2,1))
z

oz <- plot_grid(o, z, nrow=2, ncol=1, rel_heights=c(1,1.2), align="v")
oz

#Make legend
line <- c("solid","dashed")
names(line) <- c("Bacteroidetes","Firmicutes")
l <- ggplot(fbmodel %>% filter(feature=="Firmicutes" | feature=="Bacteroidetes"), aes(x=prediction, y=CI_low, linetype=feature)) +
  geom_line()+
  scale_linetype_manual(values=line)+
  theme_bw()+
  theme(legend.title = element_blank(),legend.direction="horizontal", legend.position = "bottom", legend.text = element_text(size=8))
l

le <- get_legend(l)

fig2c <- plot_grid(oz, le, nrow=2, ncol=1, rel_heights=c(1, 0.1))
fig2c

#ggsave(here("QSU_Data/deprecated/2024fig5c.jpeg"), dpi=300, h=5, w=7)


#Figure 2D: Firmcutes and Bacteroidetes Ratio
#Read in 2023 data from QSU
raw <- read.csv(here("QSU_Data/raw_data_dna_full.csv"), header=TRUE) # raw per sample information
model <- read.csv(here("QSU_Data/model_data_dna_full.csv"), header=TRUE) # means and confidence intervals for all conditions
sig <- read.csv(here("QSU_Data/sig_data_dna_full.tsv"), sep="\t", header=TRUE) # p-values and percent enrichment/depletion for all tests


# format and edit dataframes 
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values
model <- model %>% mutate(Condition = Sample_Type) # make all dataframes have a column called Sample_Type
raw <- raw %>% mutate(Preservative = substr(Sample_Type, 1, 1)) # make preservative column (first character of sample code)
raw <- raw %>% mutate(Temperature = substr(Sample_Type, 2, 2)) # make temperature column (second character of sample code)
raw <- raw %>% mutate(Sample_Type = fct_relevel(Sample_Type, "NF", "OF", "OR", "OH", "ZF", "ZR", "ZH")) # force the ordering of the conditions
raw <- raw %>% mutate(TemperatureLong = ifelse(Temperature == "F", "-80C", ifelse(Temperature == "R", "23C", "40C"))) # write out temperature
raw <- raw %>% mutate(PreservativeLong = ifelse(Preservative == "N", "None", ifelse(Preservative == "O", "Omnigene", "Zymo"))) # write out preservative
raw <- raw %>% mutate(Label = paste(TemperatureLong, PreservativeLong, sep="\n")) # make a label of temperature and preserative
raw <- raw %>% mutate(hiddenLabel=ifelse(Sample_Type == "OR", "OMNIgene", ifelse(Sample_Type== "ZR", "Zymo", ifelse(Sample_Type == "NF", "None", "")))) # make a label just for the "R" samples

raw <- mutate(raw, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                    ifelse(Sample_Type == "OF", 2, 
                                           ifelse(Sample_Type == "OR", 3, 
                                                  ifelse(Sample_Type == "OH", 4, 
                                                         ifelse(Sample_Type == "ZF", 5,
                                                                ifelse(Sample_Type == "ZR", 6, 7)))))))

r <- ggplot(raw, aes(x = reorder(Sample_Type, PlotOrder), y=ratio.of.B.F)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  #geom_violin(aes(color=Sample_Type)) + 
  geom_jitter(width=0.2, shape=16, size=1, aes(color=Sample_Type))+
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels, guide = guide_axis(angle = 45)) +
  geom_errorbar(data=model %>% filter(feature == "Bacteroidetes/Firmicutes"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0, size=1) +
  geom_point(data=model %>% filter(feature == "Bacteroidetes/Firmicutes"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2) +
  stat_pvalue_manual(sig %>% filter(feature == "Bacteroidetes/Firmicutes") %>% filter(p.adj <= 0.05), y.position=c(2.38, 2.5, 2.14, 2.26), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Bacteroidetes/Firmicutes \n") + 
  ggtitle(" ") +
  ylim(0,2.5) + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=10), axis.title.y = element_text(size=10), 
        plot.margin = unit(c(0,0,0,0), "cm"), plot.title=element_text(size=8.5))
r

#Plot Figure 2C together
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3) + 
  ylim(-0.05, 0.05) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b
r<-plot_grid(r,b, nrow = 2, ncol = 1,rel_heights=c(1,0.1), align="v", axis="l")
r

cfu_palette <- c("#acaaaf","#9a6faa","#4d9222") 
names(cfu_palette) <- c("Media", "OMNIgene", "Zymo")

cfu_styling <-   theme(panel.grid.major.x = element_blank(), text = element_text(size=10), plot.title = element_text(size=8.5))

efa <- read.csv(here("CFU/OmniZymoEfaecalisGrowthReps.csv"), header=TRUE) 
eco <- read.csv(here("CFU/OmniZymoEcoliGrowthReps.csv"), header=TRUE) 
bth <- read.csv(here("CFU/OmniZymoBthetaGrowthReps.csv"), header=TRUE) 

efa <- mutate(efa, taxon="E. faecalis")
eco <- mutate(eco, taxon="E. coli")
bth <- mutate(bth, taxon="B. theta")
efa <- mutate(efa, Condition=gsub("E faecalis", "Media", Condition))
eco <- mutate(eco, Condition = gsub("E coli", "Media", Condition))
bth <- mutate(bth, Condition = gsub("B theta", "Media", Condition))

results <- rbind(efa, eco, bth)

results <- mutate(results, Time=as.factor(Time))
bt <- ggplot(results %>% filter(taxon == "B. theta"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=1, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("B. theta"))) + 
  scale_color_manual(values = cfu_palette) + 
  scale_fill_manual(values = cfu_palette) +
  cfu_styling + 
  theme(legend.position = "bottom") +
  geom_vline(xintercept=1.5, color="grey")
bt

ec <- ggplot(results %>% filter(taxon == "E. coli"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=1, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("E. coli"))) + 
  scale_color_manual(values = cfu_palette) + 
  scale_fill_manual(values = cfu_palette) +
  cfu_styling +
  geom_vline(xintercept=1.5, color="grey")

ef <- ggplot(results %>% filter(taxon == "E. faecalis"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=1, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("E. faecalis"))) + 
  scale_color_manual(values = cfu_palette) + 
  scale_fill_manual(values = cfu_palette) +
  cfu_styling +
  geom_vline(xintercept=1.5, color="grey")

cfuleg <- get_legend(bt)

cfu_midplot <- plot_grid(bt + theme(legend.position = "none"), 
                         ec + theme(axis.title.y = element_blank(), legend.position = "none"), 
                         ef + theme(axis.title.y = element_blank(), legend.position = "none"), 
                         nrow=1, ncol=3, rel_widths = c(1.1,1,1))
cfu_midplot
cfu_plot <- plot_grid(cfu_midplot, cfuleg, nrow=2, ncol=1, rel_heights= c(1, 0.1))
cfu_plot

#ggsave(here("qPCR/mai/CFU.jpeg"), dpi=300, h=3, w=6)

#Plot Figure 2!
two<-plot_grid(absolute, cfu_plot, fig2c, r, nrow=2, ncol=2, scale=0.9, labels=c("a","b", "c", "d"), align = "vh",axis="tb", rel_widths = c(1, 1.2, 1, 1.2))
two

#Plot briefing figure
figone <-plot_grid(absolute, fig2c, nrow=1, ncol=2, scale=0.9, labels=c("a","b"), align = "vh",axis="tb", rel_widths = c(1, 1))
figone
one <- plot_grid(figone, miniplot_legend, nrow=2, ncol=1, rel_heights=c(1, 0.07), rel_widths = c(1,0.5))
one
# get legend: 
miniplot <- ggplot(raw, aes(x=Sample_Type, y=0)) + 
  geom_point(aes(color=Sample_Type)) + 
  theme_void() +
  theme(legend.position="bottom", legend.title=element_blank(), 
        legend.spacing.x = unit(0.08, "cm")) + 
  scale_color_manual(values=condition_palette, labels=condition_labels2, breaks=c("NF","OF","OR","OH","ZF","ZR","ZH")) + 
  guides(color=guide_legend(nrow=1))
miniplot

miniplot_legend <- get_legend(miniplot)

plot_grid(two, miniplot_legend, nrow=2, ncol=1, rel_heights=c(1, 0.07), rel_widths = c(1,0.9))

ggsave(here("outputs/figures/Figure5.jpeg"), dpi=300, h=8.3, w=8.5)
ggsave(here("outputs/figures/Figure5.pdf"), dpi=300, h=8.3, w=8.5)

ggsave(here("outputs/figures/BriefingFigure.jpeg"), dpi=300, h=5, w=8.5)
ggsave(here("outputs/figures/BriefingFigure.pdf"), dpi=300, h=5, w=8.5)
#####Supplementary Figures: Figure 6, 7, 8 #####
#Supplementary Figure 6
#qPCR Standard Curves
#Plate 1: Standard Curve and Absolute Count
#Create and edit unified sample datafile
qPCR_samplewellp1 <- read.csv(here("qPCR/2023/Plate1/plate1_384layout.csv"), sep=",", header=TRUE)
qPCR_datap1 <- read.csv(here("qPCR/2023/Plate1/2023qpcrcq_plate1.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp1)[colnames(qPCR_samplewellp1) == "Well384"] <- "Well"
qPCRplate1 <- merge(qPCR_samplewellp1,qPCR_datap1, by="Well")

#Make standard curve for plate 1
standard1 <- qPCRplate1 %>% filter(grepl("standard", SampleName))
standard1 <- standard1 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                  "DilutionFactor"), remove=FALSE, sep='_')
standard1 <- mutate(standard1, DilutionFactor=as.numeric(DilutionFactor))
standard1 <- mutate(standard1, Cq=as.numeric(Cq))
standard1 <- mutate(standard1, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.692*10^11
copyb <- 1.972*10^11
standard1 <- mutate(standard1, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard1 <- standard1 %>% mutate(logCopyNumber=log10(CopyNumber))
standard1 <- standard1 %>% filter(Cq<22.11) #filter out samples above limit of blank
standard1 <- standard1 %>% filter(DilutionFactor!=1e+02 & DilutionFactor!=1e+08)

stdp1model <- lm(Cq~logCopyNumber, data=standard1)
stdplate1 <- ggplot(data = standard1, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 22.11, linetype = "dotted", color = "black") +
  #annotate("text", x = max(standard1$logCopyNumber), y = 24, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7, y = 20, label = "y = 36.3 - 3.56x \nR^2=0.99", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 1") +
  xlab("log(Copy Number)") +
  ylim(5, 30) +
  xlim(4,9)
  

stdplate1

#Plate 2: Standard Curve and Absolute Count
#Create and edit unified sample datafile
qPCR_samplewellp2 <- read.csv(here("qPCR/2023/Plate2/plate2_384layout.csv"), sep=",", header=TRUE)
qPCR_datap2 <- read.csv(here("qPCR/2023/Plate2/2023qpcrcq_plate2.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp2)[colnames(qPCR_samplewellp2) == "Well384"] <- "Well"
qPCRplate2 <- merge(qPCR_samplewellp2,qPCR_datap2, by="Well")

#Make standard curve for plate 2
standard2 <- qPCRplate2 %>% filter(grepl("standard", SampleName))
standard2 <- standard2 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                  "DilutionFactor"), remove=FALSE, sep='_')
standard2 <- mutate(standard2, DilutionFactor=as.numeric(DilutionFactor))
standard2 <- mutate(standard2, Cq=as.numeric(Cq))
standard2 <- mutate(standard2, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.692*10^11
copyb <- 1.972*10^11
standard2 <- mutate(standard2, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard2 <- standard2 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
standard2 <- standard2 %>% filter(Cq<22.83) #filter out samples above limit of blank
standard2 <- standard2 %>% filter(DilutionFactor!=1e+02 & DilutionFactor!=1e+08)


stdp2model <- lm(Cq~logCopyNumber, data=standard2)
stdplate2 <- ggplot(data = standard2, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 22.83, linetype = "dotted", color = "black") +
  #annotate("text", x = max(standard2$logCopyNumber), y = 25, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7.5, y = 20, label = "y = 36.4 - 3.59x \nR^2=0.98", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 2")+
  xlab("log(Copy Number)") +
  ylim(5, 30) +
  xlim(4,9)

stdplate2

# Plate 3: Standard Curve and Absolute Count
#Create and edit unified sample datafile
qPCR_samplewellp3 <- read.csv(here("qPCR/2023/Plate3/plate3_384layout.csv"), sep=",", header=TRUE)
qPCR_datap3 <- read.csv(here("qPCR/2023/Plate3/2024qpcrcq_plate3_r3.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp3)[colnames(qPCR_samplewellp3) == "Well384"] <- "Well"
qPCRplate3 <- merge(qPCR_samplewellp3,qPCR_datap3, by="Well")

#Make standard curve for plate 3
standard3 <- qPCRplate3 %>% filter(grepl("standard", SampleName))
standard3 <- standard3 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                  "DilutionFactor"), remove=FALSE, sep='_')
standard3 <- mutate(standard3, DilutionFactor=as.numeric(DilutionFactor))
standard3 <- mutate(standard3, Cq=as.numeric(Cq))
standard3 <- mutate(standard3, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.692*10^11
copyb <- 1.972*10^11
standard3 <- mutate(standard3, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard3 <- standard3 %>% mutate(logCopyNumber=log10(CopyNumber))

#Filter out standards above limit of blank
standard3 <- standard3 %>% filter(Cq<25.16)
standard3 <- filter(standard3, DilutionFactor!=1e+02 & DilutionFactor!=1e+03 & DilutionFactor!=1e+08)

stdp3model <- lm(Cq~logCopyNumber, data=standard3)
stdplate3 <- ggplot(data = standard3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 25.16, linetype = "dotted", color = "black") +
  #annotate("text", x = max(standard3$logCopyNumber), y = 26, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7, y = 23, label = "y = 37.8 - 3.69x \nR^2=0.93", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 3") +
  xlab("log(Copy Number)") +
  ylim(5, 30) +
  xlim(4,9)

stdplate3

#Plate 4: Standard Curve and Absolute Count of missing Cq samples
#Create and edit unified sample datafile
qPCR_samplewellp4 <- read.csv(here("qPCR/2023/Plate4/plate4_384layout.csv"), sep=",", header=TRUE)
qPCR_datap4 <- read.csv(here("qPCR/2023/Plate4/2024_qpcrcq_plate4.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp4)[colnames(qPCR_samplewellp4) == "Well384"] <- "Well"
qPCRplate4 <- merge(qPCR_samplewellp4,qPCR_datap4, by="Well")

#Make standard curve for plate 4
standard4 <- qPCRplate4 %>% filter(grepl("standard", SampleName))
standard4 <- standard4 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                  "DilutionFactor"), remove=FALSE, sep='_')
standard4 <- mutate(standard4, DilutionFactor=as.numeric(DilutionFactor))
standard4 <- mutate(standard4, Cq=as.numeric(Cq))
standard4 <- mutate(standard4, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
#Note: using a new stock of standards
copyf <- 1.077*10^11
copyb <- 1.071*10^11
standard4 <- mutate(standard4, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard4 <- standard4 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
standard4 <- standard4 %>% filter(Cq<27.5) #filter out samples above limit of blank
standard4 <- standard4 %>% filter(DilutionFactor!=1e+02
                                  & DilutionFactor!=1e+08)

stdp4model <- lm(Cq~logCopyNumber, data=standard4)
stdplate4 <- ggplot(data = standard4, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 27.5, linetype = "dotted", color = "black") +
  #annotate("text", x = max(standard4$logCopyNumber), y = 29, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 6.5, y = 24, label = "y = 35.9 - 3.63x \nR^2=0.99", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 4") +
  xlab("log(Copy Number)") +
  ylim(5, 30) +
  xlim(4,9)

stdplate4

#Plate 5: Standard Curve
#Create and edit unified sample datafile
qPCR_samplewellp5 <- read.csv(here("qPCR/2023/Plate5/plate5_384layout.csv"), sep=",", header=TRUE)
qPCR_datap5 <- read.csv(here("qPCR/2023/Plate5/2024_qpcrcq_plate5.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp5)[colnames(qPCR_samplewellp5) == "Well384"] <- "Well"
qPCRplate5 <- merge(qPCR_samplewellp5,qPCR_datap5, by="Well")

#Make standard curve for plate 5
standard5 <- qPCRplate5 %>% filter(grepl("standard", SampleName))
standard5 <- standard5 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                  "DilutionFactor"), remove=FALSE, sep='_')
standard5 <- mutate(standard5, DilutionFactor=as.numeric(DilutionFactor))
standard5 <- mutate(standard5, Cq=as.numeric(Cq))
standard5 <- mutate(standard5, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
#NOTE: new stock of standards used
copyf <- 1.077*10^11
copyb <- 1.071*10^11
standard5 <- mutate(standard5, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard5 <- standard5 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
standard5 <- standard5 %>% filter(Cq<28.98) #filter out samples above limit of blank
standard5 <- standard5 %>% filter(DilutionFactor!=1e+08 & DilutionFactor!=1e+02)

stdp5model <- lm(Cq~logCopyNumber, data=standard5)
stdplate5 <- ggplot(data = standard5, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 28.98, linetype = "dotted", color = "black") +
  #annotate("text", x = max(standard5$logCopyNumber), y = 30, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7.5, y = 26, label = "y = 35.6 - 3.59x \nR^2=0.99", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 5") +
  xlab("log(Copy Number)") +
  ylim(5, 30) +
  xlim(4,9)

stdplate5

stdlegend <- ggplot(standard5, aes(x=logCopyNumber, y=Cq)) +
  geom_jitter(width=0.1, aes(color=Standard, shape=StandardReplicate)) +
  labs(
    x = "Sample",
    y = "Cq"
  ) +
  theme_bw() +
  scale_color_manual(values=standard_palette)
stdlegend
legend <- get_legend(stdlegend)
legend

allstd <- plot_grid(stdplate1, stdplate2, stdplate3, 
                    stdplate4, stdplate5,legend, nrow=2, ncol=3, 
                    rel_widths=c(1, 1, 1, 1, 1, 0.4))
allstd
#ggsave(here("qPCR/2023/outputs/allstdcurves.jpg"), dpi=300, w=10, h=7)
ggsave(here("outputs/figures/SupplementaryFigure6_allstdcurves.pdf"), dpi=300, w=10, h=7)
#Supplementary Figure 7 
#Absolute Count across Conditions, By Donor
p <- ggplot(raw, aes(x = Sample_Type, y=MicrobesPerGram)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_point(aes(color=Sample_Type),  shape=16, size=1.5) + 
  scale_color_manual(values=condition_palette, labels = c("OMNI 40°C", "OMNI 23°C", "OMNI -80°C", "No Preservative -80°C", "Zymo -80°C", "Zymo 23°C", "Zymo 40°C")) +
  scale_x_discrete(labels=condition_labels) +
  scale_y_log10() +
  theme_bw() + 
  ylab("Total Microbes per Dry Gram") + 
  xlab("Donor and Condition") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), 
        legend.position = "bottom", legend.direction="horizontal") + 
  facet_wrap(vars(Patient), nrow=2) 
p

ggsave(here("outputs/figures/SupplementaryFigure7_AbsoluteClustering.pdf"), dpi=300, w=6, h=3.5)
ggsave(here("outputs/figures/SupplementaryFigure7_AbsoluteClustering.jpeg"), dpi=300, w=6, h=3.5)

#Supplementary Figure 8
#Actino Absolute Abundance
abs<- raw[c("Sample_Type","Absolute.Abundance..Actinobacteria")]
meltabs <- melt(abs, id.vars = "Sample_Type", variable.name = "Phyla", value.name = "MicrobesPerGram")
amodel <- filter(model, feature%in%c("Absolute Abundance: Actinobacteria"))
asig <- filter(sig, feature%in%c("Absolute Abundance: Actinobacteria"))

condition_labelss <- c("40°C","23°C","-80°C","No Preservative\n-80°C","-80°C","23°C","40°C")
names(condition_labelss) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plot Order
amodel <- mutate(amodel, PlotOrder=ifelse(Sample_Type == "NF", 1, 
                                          ifelse(Sample_Type == "OF", 2, 
                                                 ifelse(Sample_Type == "OR", 3, 
                                                        ifelse(Sample_Type == "OH", 4, 
                                                               ifelse(Sample_Type == "ZF", 5,
                                                                      ifelse(Sample_Type == "ZR", 6, 7)))))))

#Plot Actinobacteria Absolute Abundance by Condition
amodel$feature <- gsub("Absolute Abundance: Actinobacteria","Actinobacteria",amodel$feature)

#Plot Actinobacteria NFvOF
nfof <- ggplot(meltabs %>% filter(Sample_Type=="NF"|Sample_Type=="OF"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labelss) +
  geom_line(data=amodel%>%filter(Sample_Type =="NF" | Sample_Type=="OF"), aes(y=prediction, group=1), color="grey")+
  geom_errorbar(data=amodel %>% filter(Sample_Type =="NF" | Sample_Type=="OF"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=0.5) +
  geom_point(data=amodel %>% filter(Sample_Type=="NF" | Sample_Type=="OF"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits = c(2e09, 3.5e10)) +
  theme_bw() + 
  ylab("**OMNIgene**<br/>Microbes/Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none",
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_markdown(size=12)) +
  ggtitle("Preservative")

nfof

#Plot Actiobacteria OFvORvOH
omni <- ggplot(meltabs %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labelss) +
  geom_line(data=amodel%>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), aes(x=reorder(Sample_Type, PlotOrder), y=prediction, group=1), color="grey")+
  geom_errorbar(data=amodel %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=0.5) +
  geom_point(data=amodel %>% filter(Sample_Type=="OF"|Sample_Type=="OR"|Sample_Type=="OH"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits = c(2e09, 3.5e10)) +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0.1,0,0), "cm"), legend.position = "none", axis.title.y = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle("Temperature")

omni

o <- plot_grid(nfof, omni, nrow=1, ncol=2, align="h", rel_widths = c(1,1))
o

#Plot Actinobacteria NFvZF
nfzf <- ggplot(meltabs %>% filter(Sample_Type=="NF"|Sample_Type=="ZF"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labelss) +
  geom_line(data=amodel%>%filter(Sample_Type =="NF" | Sample_Type=="ZF"), aes(y=prediction, group=1), color="grey")+
  geom_errorbar(data=amodel %>% filter(Sample_Type =="NF" | Sample_Type=="ZF"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=0.5) +
  geom_point(data=amodel %>% filter(Sample_Type=="NF" | Sample_Type=="ZF"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits = c(1e8,2e10)) +
  theme_bw() + 
  ylab("**Zymo**<br/>Microbes/Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", 
        axis.title.y = element_markdown(size=12))

nfzf

#Plot Actinobacteria ZFvZRvZH
zymo <- ggplot(meltabs %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(x=Sample_Type, y=MicrobesPerGram)) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labelss) +
  geom_line(data=amodel%>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), aes(x=reorder(Sample_Type, PlotOrder), y=prediction, group=1), color="grey")+
  geom_errorbar(data=amodel %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=0.5) +
  geom_point(data=amodel %>% filter(Sample_Type=="ZF"|Sample_Type=="ZR"|Sample_Type=="ZH"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction, fill=feature, color=Sample_Type), size=2.5, show.legend = FALSE) +
  scale_y_log10(limits = c(1e8,2e10)) +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), plot.margin = unit(c(0,0.1,0,0), "cm"), legend.position = "none", axis.title.y = element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())

zymo

z <- plot_grid(nfzf, zymo, nrow=1, ncol=2, align="h", rel_widths = c(1,1))
z

oz <- plot_grid(o, z, nrow=2, ncol=1, align="v")
oz

miniplot <- ggplot(raw, aes(x=Sample_Type, y=0)) + 
  geom_point(aes(color=Sample_Type)) + 
  theme_void() +
  theme(legend.title=element_blank()) + 
  scale_color_manual(values=condition_palette, labels=condition_labelss, breaks=c("NF","OF","OR","OH","ZF","ZR","ZH")) + 
  guides(color=guide_legend())
miniplot

miniplot_legend <- get_legend(miniplot)

plot_grid(oz, miniplot_legend, nrow=1, ncol=2, rel_widths = c(1,.3))
ggsave(here("outputs/figures/SupplementaryFigure8_ActinoAbsoluteAbundance.jpeg"), dpi=300, w=7, h=5)
ggsave(here("outputs/figures/SupplementaryFigure8_ActinoAbsoluteAbundance.pdf"), dpi=300, w=7, h=5)

