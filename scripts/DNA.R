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
  stat_pvalue_manual(sig %>% filter(feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(13.8, 13.2, 12.8), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Total Microbes per Gram") + 
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


#Figure2B: Changing Absolute count of Bacteroidetes and Firmicutes across conditions
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
  scale_y_log10(limits = c(1.3e11, 7e12)) +
  theme_bw() + 
  ylab("**OMNIgene**<br/>Microbes/Gram") + 
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
  scale_y_log10(limits = c(1.3e11, 7e12)) +
  theme_bw() + 
  ylab("Microbes/Gram") + 
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
  scale_y_log10(limits=c(3e10,4e12)) +
  theme_bw() + 
  ylab("**Zymo**<br/>Microbes/Gram") + 
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
  scale_y_log10(limits=c(3e10,4e12)) +
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

fig2b <- plot_grid(oz, le, nrow=2, ncol=1, rel_heights=c(1, 0.1))
fig2b

#ggsave(here("QSU_Data/Fig2Bline.jpeg"), dpi=300, h=5, w=8)


#Figure 2C: Firmcutes and Bacteroidetes Ratio
#Read in QSU ratio data as three dataframes
# rratio <- read.csv(here("QSU_Data/rawratio.csv"), header=TRUE) # raw per sample information
# mratio <- read.csv(here("QSU_Data/modelratio.csv"), header=TRUE) # means and confidence intervals for all conditions
# sratio <-read.csv(here("QSU_Data/sigratio.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests
# 
# rratio <- mutate(rratio, PlotOrder=ifelse(Sample_Type == "NF", 1, 
#                                           ifelse(Sample_Type == "OF", 2, 
#                                                  ifelse(Sample_Type == "OR", 3, 
#                                                         ifelse(Sample_Type == "OH", 4, 
#                                                                ifelse(Sample_Type == "ZF", 5,
#                                                                       ifelse(Sample_Type == "ZR", 6, 7)))))))

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
  theme(legend.position = "bottom")

ec <- ggplot(results %>% filter(taxon == "E. coli"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=1, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("E. coli"))) + 
  scale_color_manual(values = cfu_palette) + 
  scale_fill_manual(values = cfu_palette) +
  cfu_styling

ef <- ggplot(results %>% filter(taxon == "E. faecalis"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=1, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("E. faecalis"))) + 
  scale_color_manual(values = cfu_palette) + 
  scale_fill_manual(values = cfu_palette) +
  cfu_styling

cfuleg <- get_legend(bt)

cfu_midplot <- plot_grid(bt + theme(legend.position = "none"), 
          ec + theme(axis.title.y = element_blank(), legend.position = "none"), 
          ef + theme(axis.title.y = element_blank(), legend.position = "none"), 
          nrow=1, ncol=3, rel_widths = c(1.1,1,1))
cfu_midplot
cfu_plot <- plot_grid(cfu_midplot, cfuleg, nrow=2, ncol=1, rel_heights= c(1, 0.1))
cfu_plot

#Plot Figure 2!
two<-plot_grid(absolute, cfu_plot, fig2b, r, nrow=2, ncol=2, scale=0.9, labels=c("a","b", "c", "d"), align = "vh",axis="tb", rel_widths = c(1, 1.2, 1, 1.2))
two

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
    plot.title = element_text(face = "plain", size = 12),
    legend.text = element_text(size = 9),
    legend.title = element_text(size=9),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
    strip.text = element_text(color = "black", size = 12), 
    axis.text.y = element_text(size = 10), 
    axis.title.y = element_text(size=10)) + 
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
condition_labels2 <- c("OMNIgene 40°C","OMNIgene 23°C","OMNIgene -80°C","No Preservative -80°C","Zymo -80°C","Zymo 23°C","Zymo 40°C")
names(condition_labels2) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

dummy <- ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder.x), y=rel_abundance)) +
  geom_bar(stat="identity", aes(fill=Condition)) + 
  scale_fill_manual(values=condition_palette, labels=condition_labels2, breaks=c("NF","OF","OR","OH","ZF","ZR","ZH")) + 
  theme(legend.direction="horizontal", legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(legend.text = element_text(size=7.5), legend.title= element_blank())  
dummy
b <- get_legend(dummy)
fig3a <- plot_grid(r, b, nrow=2, ncol=1, rel_heights = c(1, 0.08), rel_widths = c(1, 0.7))
fig3a


#3B: Shannon entropy across conditions
#Change column name
colnames(raw)[colnames(raw) == "Shannon.Entropy"] <- "Shannon"

#Shannon entropy
se <- ggplot(raw , aes(Sample_Type, Shannon)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=1) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels, guide = guide_axis(angle = 45)) +
  geom_errorbar(data=model %>% filter(feature == "Shannon Entropy"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Shannon Entropy"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  ylim(1,4.5) +
  stat_pvalue_manual(sig %>% filter(feature == "Shannon Entropy") %>% filter(p.adj <= 0.05), y.position=c(4.1, 4.3, 3.7,3.9), 
                     tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylab("Shannon Entropy") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), axis.title.y = element_text(size=10), 
        plot.margin = unit(c(0,0,0,0), "cm"))
se

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3.5) +
  ylim(-0.1, 0.1) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"), text=element_text(size=16))
b
shannon <- plot_grid(se,b, nrow=2, ncol=1, rel_heights=c(1,0.1), align="v", axis='l')
shannon


#inverse simpson 
simpson <- ggplot(raw , aes(Sample_Type, Inv.Simpson)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=1) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Inv Simpson"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Inv Simpson"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  ylim(0,20) +
  stat_pvalue_manual(sig %>% filter(feature == "Inv Simpson") %>% filter(p.adj <= 0.05), y.position=c(19, 20, 18,19), 
                      tip.length=0) + #HERE
  theme_bw() + 
  ylab("Inverse Simpson Index") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
simpson

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size = 3.5) +
  ylim(-0.1, 0.1) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"))
b
simpsonplot <- plot_grid(simpson,b, nrow=2, ncol=1, rel_heights=c(1,0.1), scale=0.95, align="v", axis='l')
simpsonplot

miniplot_simpson <- ggplot(raw, aes(x=Sample_Type, y=0)) + 
  geom_point(aes(color=Sample_Type)) + 
  theme_void() +
  theme(legend.title=element_blank()) + 
  scale_color_manual(values=condition_palette, labels=condition_labels2, breaks=c("NF","OF","OR","OH","ZF","ZR","ZH")) + 
  guides(color=guide_legend())
miniplot_simpson

miniplot_legend_simpson <- get_legend(miniplot_simpson)
plot_grid(simpsonplot, miniplot_legend_simpson, ncol=2, nrow=1, rel_widths=c(1, 0.4))

ggsave(here("outputs/figures/SupplementaryFigure2_InverseSimpson.pdf"), width=5.7, h=3.5, dpi=300)
ggsave(here("outputs/figures/SupplementaryFigure2_InverseSimpson.jpeg"), width=5.7, h=3.5, dpi=300)


#Figure 3E:Phylum forest plotting 

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
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Bacteroidetes"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0, size=0.8) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Bacteroidetes"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction*100, color=Sample_Type), size=2) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Bacteroidetes") %>% filter(p.adj <= 0.05), y.position=c(88, 95, 68, 75), 
                     tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylim(0,100) + 
  ylab("Relative Abundance (%)  ") + 
  xlab("Bacteroidetes") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), text = element_text(size=12),
        legend.position = "none", axis.title.y=element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(size=10))
bact

firm <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Firmicutes*100)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Firmicutes"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0, size=0.8) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Firmicutes"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=2) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Firmicutes") %>% filter(p.adj <= 0.05), y.position=c(90, 95, 70,75), 
                     tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylim(0,100) + 
  ylab("Relative Abundance") + 
  xlab("Firmicutes") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_text(size=10))

act <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Actinobacteria*100)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Actinobacteria"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0, size=0.8) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Actinobacteria"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=2) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Actinobacteria") %>% filter(p.adj <= 0.05), y.position=c(15, 17, 11, 13, 13), 
                     tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylab("Relative Abundance") + 
  xlab("Actinobacteria") + 
  ylim(0,20)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.title.x= element_text(size=10))
act
vir <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Viruses*100)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Viruses"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0, size=0.8) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Viruses"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=2) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Viruses") %>% filter(p.adj <= 0.05), y.position=c(0.45, 0.38, 0.41, 0.38, 0.41),
                     tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylab("Relative Abundance") + 
  xlab("Viruses") + 
  ylim(0,0.5)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_text(size=10))

vir

fungi <- ggplot(raw , aes(reorder(Sample_Type, PlotOrder), Relative.Abundance..Fungi*100)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Relative Abundance: Fungi"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low*100, ymax=CI_high*100), width=0, size=0.8) +
  geom_point(data=model %>% filter(feature == "Relative Abundance: Fungi"), inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=prediction*100, color=Sample_Type), size=2) +
  stat_pvalue_manual(sig %>% filter(feature == "Relative Abundance: Fungi") %>% filter(p.adj <= 0.05), y.position=c(0.2, 0.15),
                     tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylab("Relative Abundance") + 
  xlab("Fungi") + 
  ylim(0,0.5) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", 
        text = element_text(size=12),  axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_text(size=10))

fungi

pdiff <- plot_grid(bact, firm, act, vir, fungi, nrow=1, ncol=5, rel_widths = c(1.15, 1, 1, 1, 1))
pdiff


dummy <- ggplot(bracken_pheno, aes(x=reorder(Sample, PlotOrder.x), y=rel_abundance)) +
  geom_point(stat="identity", aes(color=Condition)) + 
  theme_bw() + 
  scale_color_manual(values=condition_palette, labels=condition_labels2) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(text = element_text(size=16)) 
dummy
b <- get_legend(dummy)

pdiff_all <- plot_grid(pdiff, b, nrow=1, ncol=2, rel_widths = c(1, 0.24))
pdiff_all

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
  geom_jitter(width=0.2, aes(color=group1),  shape=16, size=1) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels, guide = guide_axis(angle = 45)) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=Mean), size=2.5) +
  ylim(0,0.75) +
  stat_pvalue_manual(bc_sig_toNF %>% filter(p.adj <= 0.05), y.position=c(0.75, .65, .55, .65), 
                      tip.length=0, size=3, label = "p.signif") + # HERE
  theme_bw() + 
  ylab("Bray-Curtis Dissimilarity") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), axis.title.y = element_text(size=10), 
        plot.margin = unit(c(0,0,0,0), "cm"))
bc_plot

#Legend formatting
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1") %>% filter(Sample_Type != "NF"), aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3.5) +
  ylim(-0.05, 0.05) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0, 0.0, 0), "cm"), text = element_text(size=14))
b
bc_full <- plot_grid(bc_plot,b, nrow=2, ncol=1, rel_heights=c(1,0.1), align="v", axis='l')
bc_full

#Figure3C: Relative Richness 0.01%
p <- ggplot(raw, aes(x = Sample_Type, y=Richness.0.01.)) +
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) +
  geom_jitter(width=0.2, aes(color=Sample_Type), shape=16, size=1) +
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels, guide = guide_axis(angle = 45)) +
  geom_errorbar(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  #stat_pvalue_manual(sig %>% filter(Feature == "Richness 0.01%"), y.position=c(127, 124, 118, 121, 118, 121),
  #                   tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sig %>% filter(feature == "Richness 0.01%") %>% filter(p.adj <= 0.05), y.position=c(127, 124, 121),
                     tip.length=0, size=3, size=3, label = "p.signif") + # HERE
  theme_bw() +
  ylab("Number of Genera\n>0.01% Relative Abundance") +
  ylim(50,130) +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), axis.title.y = element_text(size = 10), plot.margin = unit(c(0,0,0,0), "cm"))

p
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) +
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3.5) +
  ylim(-0.5, 0.5) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

richness <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.1 ), align="v", axis='l')
richness

#Plot Figure 3!
a <- plot_grid(fig3a, nrow=1,ncol=1,scale=1, labels=c("a"))
a
b <- plot_grid(pdiff, nrow=1,ncol=1, scale=1, rel_widths=c(1),labels=c("b"), align="v", axis="tb")
b
c <- plot_grid(shannon, richness, bc_full, nrow=1,ncol=3, scale=0.9, rel_widths=c(1.1, 1.2, 1),labels=c("c","d", "e"), align="v", axis="tb")
c
three<-plot_grid(a, NULL, b, NULL, c, nrow=5, ncol=1, rel_widths=c(1,1, 1, 1, 1), rel_heights=c(1,0.1, 0.6,0.05, 0.8), align="h", axis="lr")
three

ggsave(here("outputs/figures/Figure2_test.pdf"), dpi=300, h=9, w=8.5)
ggsave(here("outputs/figures/Figure2_test.jpeg"), dpi=300, h=9, w=8.5)


##### SUPPLEMENT #####
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
  geom_jitter(width=0.2, aes(color=Condition),  shape=16, size=1) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=Mean), size=2.5) +
  ylim(0,0.75) +
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
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3.5) +
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
  geom_jitter(width=0.2, aes(color=Condition),  shape=16, size=1) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=bc_model, inherit.aes=FALSE, aes(x=Condition, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=bc_model, inherit.aes=FALSE, aes(x=Condition, y=estimate), size=2.5) +
  ylim(0,0.75) +
  stat_pvalue_manual(bc_sig %>% filter(p.adj <= 0.05), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Bray-Curtis Dissimilarity") + 
  ggtitle("Metatranscriptomic") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
bc_plot

foodf <- data.frame(xvals = c(0.3, 2, 4.6), labels=c("None", "OMNIgene", "Zymo"))

test <- ggplot(foodf, aes(x=xvals, y=0)) + 
  geom_text(aes(y=0, label=labels), fontface = "bold", size=3.5) + 
  ylim(-0.5, 0.5) +
  xlim(0,6) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme_void()
test
bc_full_rna <- plot_grid(bc_plot,test, nrow=2, ncol=1, rel_heights=c(1,0.05), align="v", axis='l')
bc_full_rna


bcplot <- plot_grid(bc_full_dna, bc_full_rna, nrow=1, ncol=2, scale=0.9, labels=c("a", "b"))
plot_grid(bcplot, miniplot_legend, ncol=1, nrow=2, rel_heights=c(1,0.07))
ggsave(here("outputs/figures/SupplementaryFigure4_BrayCurtis.pdf"), dpi=300, w=8.5, h=4)
ggsave(here("outputs/figures/SupplementaryFigure4_BrayCurtis.jpeg"), dpi=300, w=8.5, h=4)

##### GENUS HEATMAP
#genus_sig <- read.csv(here("QSU_Data/sig_data_dna_full.tsv"), sep="\t", header=TRUE)
genus_sig <- sig
genus_sig <- genus_sig %>% mutate(PercentFormatted = as.numeric(gsub("%.*", "", percentchange)))
genus_sig <- genus_sig %>% mutate(PFormatted = as.numeric(ifelse(p == "< 0.001", "0.001", p)))
genus_sig <- genus_sig %>% filter(grepl(":", feature)) %>% filter(!grepl("Absolute ", feature)) %>% 
  filter(!grepl("Relative ", feature))
genus_sig <- genus_sig %>% mutate(feature = ifelse(feature == "Firmicutes: unclassified Firmicutes sensu stri...", "Firmicutes: Firmicutes s.s.", feature))


genus_sig <- genus_sig %>% mutate(Condition = group2)
genus_sig <- genus_sig %>% mutate(y.position = ifelse(Condition == "OF", 2, 
                                                      ifelse(Condition=="ZF", 1, 
                                                             ifelse(Condition=="OR", 4, 
                                                                    ifelse(Condition == "OH", 3, 
                                                                           ifelse(Condition == "ZR", 2, 1))))))

#raw_abundance <- read.csv(here("QSU_Data/model_data_dna_full.csv"), header=TRUE)
raw_abundance <- model
raw_abundance <- raw_abundance %>% mutate(Condition=Sample_Type)
raw_abundance <- raw_abundance %>% filter(Condition == "NF") %>% select(feature, prediction)
names(raw_abundance) <- c("feature", "NFMean")
raw_abundance <- raw_abundance %>% filter(grepl(":", feature)) %>% filter(!grepl("Absolute ", feature)) %>% 
  filter(!grepl("Relative ", feature))
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: unclassified Clostridiales.*", "Firmicutes: unclassified Clostridiales (miscel...", feature))
raw_abundance <- raw_abundance %>% mutate(feature = gsub("Firmicutes: unclassified Firmicutes sensu.*", "Firmicutes: unclassified Firmicutes sensu stri...", feature))
raw_abundance <- raw_abundance %>% mutate(feature = ifelse(feature == "Firmicutes: unclassified Firmicutes sensu stri...", "Firmicutes: Firmicutes s.s.", feature))
raw_abundance <- raw_abundance %>% mutate(Phylum = gsub(":.*", "", feature))
raw_abundance <- raw_abundance %>% 
  group_by(Phylum) %>% 
  arrange(Phylum, desc(NFMean))
raw_abundance$x.position <- as.numeric(row.names(raw_abundance))


genus_sig <- merge(genus_sig, raw_abundance, by="feature", all.x = TRUE)
genus_sig <- genus_sig %>% mutate(Phylum = gsub(":.*", "", feature))
genus_sig <- genus_sig %>% mutate(Phylum = ifelse(Phylum == "null", "Virus", Phylum))
genus_sig <- genus_sig %>% mutate(Plabel = ifelse(PFormatted <= 0.05, "*", ""))
genus_sig <- genus_sig %>% mutate(GenusStrange = ifelse((str_count(feature, pattern=" ") > 1) & feature != "null: crAss-like viruses", "*", ""))
genus_sig <- genus_sig %>% mutate(Genus = gsub(" sensu.*", "", gsub("unclassified ", "", gsub(".*: ", "", gsub(" \\(misc.*", "", gsub(": environmental.*", "", feature))))))
genus_sig <- genus_sig %>% mutate(Genus = paste(Genus, GenusStrange, sep=""))

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
  geom_vline(xintercept=8.5) + 
  geom_hline(yintercept=1.5) +
  theme(panel.border = element_rect(size=1))

preservative

z <- get_legend(preservative)

temperature <- ggplot(genus_sig %>% filter(Condition != "OF" & Condition != "ZF"), aes(x=x.position, y=y.position, fill=PercentFormatted)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0), breaks=genus_sig$x.position, labels=genus_sig$Genus) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1,2, 3, 4), labels=c("Zymo 40°C", "Zymo 23°C", "OMNI 40°C", "OMNI 23°C")) + 
  labs(title = "Temperature Effect (Relative to -80°C for each preservative)") +
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(),
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-81, 81)) +
  geom_point(data = genus_sig %>% filter(Condition != "OF" & Condition !="ZF") %>% filter(Plabel == "*"), size=1, color="white") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=8.5) + 
  geom_hline(yintercept=2.5) +
  theme(panel.border = element_rect(size=1))


temperature

raw_abund <- ggplot(genus_sig %>% filter(Condition == "ZF"), aes(x=x.position, y=y.position, fill=NFMean*100)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), breaks=c(1), labels=c("Abundance")) + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), 
                     title = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradientn(colours = c("#F3F3F3", "black")) +
  labs(fill="Baseline %\nAbundance")+ 
  geom_vline(xintercept=3.5) + 
  geom_vline(xintercept=8.5) + 
  theme(panel.border = element_rect(size=1))
raw_abund

z2 <- get_legend(raw_abund)


main_andabund <- plot_grid(raw_abund + theme(legend.position="none"), 
                           preservative + theme(legend.position="none"), 
                           temperature + theme(legend.position="none"), 
                           nrow=3, ncol=1, align="v", axis="lr", rel_heights=c(1,3,8))

main_andabund # works well at 11.5 x 3

legend_plot <- plot_grid(z2, z, ncol=1, nrow=2, align = "v")
legend_plot

plot_grid(main_andabund, legend_plot, nrow=1, ncol=2, rel_widths=c(1, 0.1)) # seems to work well at 12x3
ggsave(here("outputs/figures/SupplementaryFigure6_heatmap_raw.pdf"), dpi=300, w=12, h=4)
ggsave(here("outputs/figures/SupplementaryFigure6_heatmap_raw.jpeg"), dpi=300, w=12, h=4)


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
  scale_y_log10(limits = c(2e10, 3.5e11)) +
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
  scale_y_log10(limits = c(2e10, 3.5e11)) +
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
  scale_y_log10(limits = c(1e9,2e11)) +
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
  scale_y_log10(limits = c(1e9,2e11)) +
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

plot_grid(oz, miniplot_legend_simpson, nrow=1, ncol=2, rel_widths = c(1,.3))
ggsave(here("outputs/figures/SupplementaryFigure2_ActinoAbsoluteAbundance.jpeg"), dpi=300, w=7, h=5)
ggsave(here("outputs/figures/SupplementaryFigure2_ActinoAbsoluteAbundance.pdf"), dpi=300, w=7, h=5)

#### DNA and RNA concentration ####
# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#dna <- read.csv(here("QSU_Data/raw_data_dna_full.csv"), head=TRUE) 
dna <- raw
dna <- dna %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove = FALSE)
dna <- dna %>% mutate(hiddenLabel=ifelse(Condition == "OR", "OMNIgene", ifelse(Condition== "ZR", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples

#model <- read.csv(here("QSU_Data/model_data_dna_full.csv"), header=TRUE) # means and confidence intervals for all conditions

model <- model %>% mutate(Condition = Sample_Type)
#sig <- read.csv(here("QSU_Data/sig_data_dna_full.tsv"), header=TRUE, sep="\t") # p-values and percent enrichment/depletion for all tests
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values

dna <- mutate(dna, PlotOrder=ifelse(Condition == "NF", 1, 
                                    ifelse(Condition == "OF", 2, 
                                           ifelse(Condition == "OR", 3, 
                                                  ifelse(Condition == "OH", 4,
                                                         ifelse(Condition == "ZF", 5,
                                                                ifelse(Condition == "ZR", 6, 7)))))))


p <- ggplot(dna, aes(x = reorder(Condition, PlotOrder), y=DNAConcentration, fill=Condition)) + 
  geom_jitter(width=0.2, aes(color=Condition), shape=16, size=1) +
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) +
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  theme_bw() +
  ylab("DNA Concentration (ng/ul)") +
  geom_errorbar(data=model %>% filter(feature == "DNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "DNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  stat_pvalue_manual(sig %>% filter(feature == "DNAConcentration") %>% filter(group1 == "NF" | group1 == "OF" | group1 == "ZF") %>% filter(p.adj <= 0.05),
                     tip.length=0, label = "p.signif", inherit.aes = FALSE, y.position = c(160, 140, 150)) +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p
b <- ggplot(dna %>% filter(Donor == "D01" & Replicate == "R1"), aes(x=reorder(Condition, PlotOrder), y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold", size=3.5) + 
  ylim(-0.5, 0.5) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

dna_plot <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='lr')
dna_plot

# RNA

rna <- read.csv(here("QSU_Data/raw_data_rna_full.csv"), head=TRUE) #FIXME WHy only 170 - ten samples have 0 concentration
rna <- rna %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove = FALSE)
rna <- rna %>% mutate(hiddenLabel=ifelse(Condition == "OR", "OMNIgene", ifelse(Condition== "ZR", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples

model <- read.csv(here("QSU_Data/model_data_rna_full.csv"), header=TRUE) # means and confidence intervals for all conditions
model <- model %>% mutate(Condition = Sample_Type)
sig <- read.csv(here("QSU_Data/sig_data_rna_full.tsv"), header=TRUE, sep="\t") # p-values and percent enrichment/depletion for all tests
sig <- sig %>% mutate(y.position=15) # significance table requires a column called y.position for plotting - 15 is a dummy value
sig <- sig %>% mutate(p.signif=ifelse(p.signif == "", "ns", p.signif)) # add a label called "ns" for non-significant p-values

rna <- mutate(rna, PlotOrder=ifelse(Condition == "NF", 1, 
                                    ifelse(Condition == "OF", 2, 
                                           ifelse(Condition == "OR", 3, 
                                                  ifelse(Condition == "ZF", 4,
                                                         ifelse(Condition == "ZR", 5, 6))))))

p <- ggplot(rna, aes(x = reorder(Condition, PlotOrder), y=RNAConcentration)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=3.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Condition), shape=16, size=1) + 
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  theme_bw() + 
  ylab("RNA Concentration (ng/ul)") + 
  ylim(0, 350) +
  geom_errorbar(data=model %>% filter(feature == "RNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "RNAConcentration"), inherit.aes=FALSE, aes(x=Sample_Type, y=prediction), size=2.5) +
  stat_pvalue_manual(sig %>% filter(feature == "RNAConcentration") %>% filter(p.adj <= 0.05), y.position = c(325, 350, 275, 275),
                     tip.length=0, label = "p.signif") +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

p

foodf <- data.frame(xvals = c(0.3, 2, 4.6), labels=c("None", "OMNIgene", "Zymo"))

test <- ggplot(foodf, aes(x=xvals, y=0)) + 
  geom_text(aes(y=0, label=labels), fontface = "bold", size=3.5) + 
  ylim(-0.5, 0.5) +
  xlim(0,6) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme_void()
test


rna_plot <- plot_grid(p,test, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='lr')
rna_plot

concentration_plot <- plot_grid(dna_plot, rna_plot, ncol=2, nrow=1, scale=0.9, labels=c("a", "b"))
plot_grid(concentration_plot, miniplot_legend, nrow=2, ncol=1, rel_heights=c(1,0.07))

ggsave(here("outputs/figures/SupplementaryFigure1_Concentration.jpeg"), dpi=300, w=8.5, h=4)
ggsave(here("outputs/figures/SupplementaryFigure1_Concentration.pdf"), dpi=300, w=8.5, h=4)


##### REVISIONS #####
#Absolute Count across Conditions, By Donor
raw_extended <- raw %>% mutate(Donor = gsub("-.*", "", Sample))
p <- ggplot(raw_extended, aes(x = Sample_Type, y=MicrobesPerGram)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_point(aes(color=Sample_Type),  shape=16, size=1.5) + 
  scale_color_manual(values=condition_palette, labels = c("OMNI 40°C", "OMNI 23°C", "OMNI -80°C", "No Preservative -80°C", "Zymo -80°C", "Zymo 23°C", "Zymo 40°C")) +
  scale_x_discrete(labels=condition_labels) +
  scale_y_log10() +
  theme_bw() + 
  ylab("Total Microbes per Gram") + 
  xlab("Donor and Condition") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        text = element_text(size=12), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.title = element_blank(), 
        legend.position = "bottom", legend.direction="horizontal") + 
  facet_wrap(vars(Donor), nrow=2) 
p

ggsave(here("outputs/figures/SupplementaryFigure6_AbsoluteClustering.pdf"), dpi=300, w=6, h=3.5)
ggsave(here("outputs/figures/SupplementaryFigure6_AbsoluteClustering.jpeg"), dpi=300, w=6, h=3.5)

#qPCR Standard Curves
#Plate 1: Standard Curve
#Create and edit unified sample datafile
qPCR_samplewellp1 <- read.csv(here("qPCR/Plate_1/Plate1_Layout.csv"), sep=",", header=TRUE)
qPCR_datap1 <- read.csv(here("qPCR/Plate_1/HotPooqPCR_Plate1.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp1)[colnames(qPCR_samplewellp1) == "Well384"] <- "Well"
qPCRplate1 <- merge(qPCR_samplewellp1,qPCR_datap1, by="Well")

#Make standard curve for plate 1
qPCRplate1 <- qPCRplate1 %>% filter(grepl("standard", SampleName))
qPCRplate1 <- filter(qPCRplate1, PCR_Replicate!="Rep4")
qPCRplate1 <- qPCRplate1 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate1 <- mutate(qPCRplate1, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate1 <- mutate(qPCRplate1, Cq=as.numeric(Cq))
qPCRplate1 <- mutate(qPCRplate1, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate1 <- mutate(qPCRplate1, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

#Standard curve using both B vulgatus and F prausnitizii data
plate1 <- ggplot(data = qPCRplate1, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="black") + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE, label.y = "top", label.x = "right") +         
  geom_point(aes(color=Standard)) +
  theme_bw() +
  labs(title="Plate 1") +
  theme(legend.position="none")
plate1

#Plate 2: Standard Curve
#Create and edit unified sample datafile
qPCR_samplewellp2 <- read.csv(here("qPCR/Plate_2/qPCR_plate2_layout.csv"), sep=",", header=TRUE)
qPCR_datap2 <- read.csv(here("qPCR/Plate_2/HotPooqPCR_Plate2.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp2)[colnames(qPCR_samplewellp2) == "Well384"] <- "Well"
qPCRplate2 <- merge(qPCR_samplewellp2,qPCR_datap2, by="Well")
#Make standard curve for plate 2
qPCRplate2 <- qPCRplate2 %>% filter(grepl("standard", SampleName))
qPCRplate2 <- filter(qPCRplate2, PCR_Replicate!="Rep4")
qPCRplate2 <- qPCRplate2 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate2 <- mutate(qPCRplate2, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate2 <- mutate(qPCRplate2, Cq=as.numeric(Cq))
qPCRplate2 <- mutate(qPCRplate2, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate2 <- mutate(qPCRplate2, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Standard curve using both B vulgatus and F prausnitizii data
qPCRplate2 <- qPCRplate2 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

plate2 <- ggplot(data = qPCRplate2, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="black") + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE, label.y = "top", label.x = "right") +         
  geom_point(aes(color=Standard)) +
  theme_bw() +
  labs(title="Plate 2") +
  theme(legend.position="none")
plate2

# Plate 3: Standard Curve
#Create and edit unified sample datafile
qPCR_samplewellp3 <- read.csv(here("qPCR/Plate_3/qPCR_plate3_layout.csv"), sep=",", header=TRUE)
qPCR_datap3 <- read.csv(here("qPCR/Plate_3/HotPooqPCR_Plate3.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp3)[colnames(qPCR_samplewellp3) == "Well384"] <- "Well"
qPCRplate3 <- merge(qPCR_samplewellp3,qPCR_datap3, by="Well")

#Make standard curve for plate 3
qPCRplate3 <- qPCRplate3 %>% filter(grepl("standard", SampleName))
qPCRplate3 <- filter(qPCRplate3, PCR_Replicate!="Rep4")
qPCRplate3 <- qPCRplate3 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate3 <- mutate(qPCRplate3, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate3 <- mutate(qPCRplate3, Cq=as.numeric(Cq))
qPCRplate3 <- mutate(qPCRplate3, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate3 <- mutate(qPCRplate3, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Standard curve using both B vulgatus and F prausnitizii data
qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

plate3 <- ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="black") + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE, label.y = "top", label.x = "right") +         
  geom_point(aes(color=Standard)) +
  theme_bw() +
  labs(title="Plate 3")  +
  xlab("log(Copy Number)") +
  theme(legend.position="none")
plate3

#Plot legend
p <- ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=Standard)) +
  theme_bw()
l <- get_legend(p)

#Plot standard curves together
sup6 <- plot_grid(plate1, plate2, plate3, nrow=1, ncol=3, labels =c("a","b","c"))
sup6

s6 <- plot_grid(sup6, l, nrow=1, ncol=2, rel_widths=c(1, 0.1))
s6

ggsave(here("outputs/figures/SupplementaryFigure6_StandardCurves.jpeg"), dpi=300, w=13, h=5)
ggsave(here("outputs/figures/SupplementaryFigure6_StandardCurves.pdf"), dpi=300, w=13, h=5)

