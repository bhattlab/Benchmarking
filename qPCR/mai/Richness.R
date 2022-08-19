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
colnames(model)[colnames(model) == "prediction"] <- "Mean"

#####Figure 2: Total microbes/gram#####
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
  stat_pvalue_manual(sig %>% filter(feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(13.8, 13, 12.8), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"))
p

b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

absolute <- plot_grid(p,b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='l')
absolute

##### Richness 
p <- ggplot(raw, aes(x = Sample_Type, y=Richness.100M)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  geom_jitter(width=0.2, aes(color=Sample_Type), shape=16, size=2) + 
  scale_fill_manual(values=condition_palette) +
  scale_color_manual(values=condition_palette) + 
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=model %>% filter(feature == "Richness 100M"), inherit.aes=FALSE, aes(x=Sample_Type, ymin=CI_low, ymax=CI_high), width=0.1, size=1) +
  geom_point(data=model %>% filter(feature == "Richness 100M"), inherit.aes=FALSE, aes(x=Sample_Type, y=Mean), size=2.5) +
  #stat_pvalue_manual(sig %>% filter(Feature == "Richness 0.01%"), y.position=c(127, 124, 118, 121, 118, 121), 
  #                   tip.length=0, label = "p.signif") +
  stat_pvalue_manual(sig %>% filter(feature == "Richness 100M") %>% filter(p.adj <= 0.05), y.position=c(1400, 1300), 
                     tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Number of Genera above \n100M counts/gram") + 
  ylim(0,1500) + 
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

plot_grid(absolute, richness, nrow=1, ncol=2, scale=0.9, labels=c("A","B"))
ggsave(here("QSU_Data/Figure2.pdf"), w=12, h=4, dpi=300)
ggsave(here("QSU_Data/Figure2.png"), w=12, h=4, dpi=300)

####Figure 2C: Firmcutes and Bacteroidetes Ratio
abs <- raw[c("Sample_Type","Absolute.Abundance..Firmicutes", "Absolute.Abundance..Bacteroidetes")]
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

fb <- ggplot(meltabs, aes(x=Sample_Type, y=MicrobesPerGram)) + 
  geom_vline(aes(xintercept=1.5), alpha=0.2, size=0.3) +
  geom_vline(aes(xintercept=4.5), alpha=0.2, size=0.3) + 
  #geom_jitter( aes(color=Sample_Type, shape=Phyla), position=position_jitterdodge(), size=2) + 
  scale_color_manual(values=condition_palette, guide="none") +
  scale_x_discrete(labels=condition_labels) +
  geom_errorbar(data=fbmodel, inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), ymin=CI_low, ymax=CI_high, linetype=feature), size=1, position=position_dodge(width=1), width=0) +
  geom_point(data=fbmodel, inherit.aes=FALSE, aes(x=reorder(Sample_Type, PlotOrder), y=Mean, fill=feature, color=Sample_Type), size=2.5, position=position_dodge(width=1), show.legend = FALSE) +
  scale_y_log10() +
  # stat_pvalue_manual(sig %>% filter(Feature == "MicrobesPerGram") %>% filter(p.adj <= 0.05), y.position=c(14.5, 14.3, 14, 13.8, 13, 12.8), 
  #                    tip.length=0, label = "p.signif") +
  #stat_pvalue_manual(fbsig %>% filter(p.adj <= 0.05), y.position=(c(13,13.25,13.5,13.75)),
                     #tip.length=0, label = "p.signif") +
  theme_bw() + 
  ylab("Microbes per Gram") + 
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      text = element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"), legend.position = c(0.85, 0.85), legend.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
fb


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
  geom_jitter(width=0.2, shape=16, size=2, aes(color=Sample_Type))+
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

# b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
#   geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
#   ylim(-0.05, 0.05) +
#   theme_void() + 
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# b

# rb <- plot_grid(r,b, nrow=2, ncol=1, rel_heights=c(1, 0.1), align="v", axis='lr')
# rb

ggsave(here("QSU_Data/Figure2ratio.jpeg"), w=5, h=4, dpi=300)

#Plot Figure 3C together
b <- ggplot(raw %>% filter(Patient == "D01" & Replication == "R1"), aes(x=Sample_Type, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.05, 0.05) +
  theme_void() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b
threec<-plot_grid(fb, r,b, nrow = 3, ncol = 1,rel_heights=c(1, 2,0.2), align="v", axis="l")
threec

ggsave(here("QSU_Data/Figure2C.jpeg"), w=6.5, h=6.5, dpi=300)

#####Richness Analysis#####
#For OMNI and Zymo which have lower richness, which genera are missing?
genus <- model %>% filter(!grepl("Richness", "Shannon", "Inv", "Relative", "DNAConcentration", "MicrobesPerGram"))
mean <- model %>% filter(grepl("Relative Abundance", feature))
mean <- mean %>% filter(Mean>0.0001) #filter for taxa above 0.01% Richness


