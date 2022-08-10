library(ggplot2)
library(ggpubr)
library(here)
library(dplyr)
library(forcats)
library(gtable)
library(cowplot)

# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# read in the QSU data as three separate dataframes 
raw <- read.csv(here("QSU_Data/raw_data_all.csv"), header=TRUE) # raw per sample information
model <- read.csv(here("QSU_Data/raw_data_model_means.csv"), header=TRUE) # means and confidence intervals for all conditions
sig <- read.csv(here("QSU_Data/raw_data_significance.csv"), header=TRUE) # p-values and percent enrichment/depletion for all tests

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

p <- ggplot(raw, aes(x = Sample_Type, y=MicrobesPerGram)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type),  shape=16, size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "MicrobesPerGram"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "MicrobesPerGram"), inherit.aes=FALSE, aes(x=Sample_Type, y=Mean), size=2.5) +
  scale_y_log10() +
  # stat_pvalue_manual(sig %>% filter(Feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(14.5, 14.3, 14, 13.8, 13, 12.8), 
  #                    tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sig %>% filter(Feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(13.8, 13, 12.8), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))

b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

absolute <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='l')



##### Richness 
p <- ggplot(raw, aes(x = Sample_Type, y=Richness.0.01.)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type), shape=16, size=2) + 
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Richness 0.01%"), inherit.aes=FALSE, aes(x=Sample_Type, y=Mean), size=2.5) +
  #stat_pvalue_manual(sig %>% filter(Feature == "Richness 0.01%"), y.position=c(127, 124, 118, 121, 118, 121), 
  #                   tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sig %>% filter(Feature == "Richness 0.01%") %>% filter(p.adj <= 0.05), y.position=c(127, 124, 121), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Number of Genera") + 
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


plot_grid(absolute, richness, nrow=1, ncol=2, scale=0.9, labels=c("A","B"))
ggsave(here("QSU_Data/Figure2.pdf"), w=12, h=4, dpi=300)
ggsave(here("QSU_Data/Figure2.png"), w=12, h=4, dpi=300)



##### FIGURE 3 #####
