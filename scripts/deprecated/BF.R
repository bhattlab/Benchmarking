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
condition_labels <- c("40°C","23°C","-80°C","-80°C\n No Preservative","-80°C","23°C","40°C")
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

#####Figure 2 BF#####
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
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),) +
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

ggsave(here("QSU_Data/Fig2Bline.jpeg"), dpi=300, h=5, w=8)
