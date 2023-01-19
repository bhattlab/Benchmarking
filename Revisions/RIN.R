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

# build a named color palette for each condition
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

# build a named set of labels for each condition (just the temperatures)
condition_labels <- c("40°C","23°C","-80°C","-80°C","-80°C","23°C","40°C")
names(condition_labels) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

rin <- read.csv(here("Revisions/122022RIN.csv"), header=TRUE) 
rin <- rin %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)
rin <- rin %>% mutate(hiddenLabel=ifelse(Condition == "OH", "OMNIgene", ifelse(Condition== "ZH", "Zymo", ifelse(Condition == "NF", "None", "")))) # make a label just for the "R" samples

#Plot RIN by condition and shape by donor
r <- ggplot(rin, aes(x=factor(Condition, levels=c('NF','OF','OR','OH','ZF','ZR','ZH')), y=RIN)) + 
  geom_vline(aes(xintercept=1.5), color="lightgrey") + 
  geom_vline(aes(xintercept=4.5), color="lightgrey") +
  geom_point(width=0.2, aes(color=Condition, shape=Donor), size=2) + 
  scale_color_manual(values=condition_palette) +
  scale_x_discrete(labels=condition_labels) +
  ylab("RIN") +
  theme_bw() +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        text = element_text(size=12)) + 
  guides(color = FALSE)
r

#Plot x axis preservative title
b <- ggplot(rin %>% filter(Donor == "D01" & Replicate == "R1"), aes(x=Condition, y=0)) + 
  geom_text(aes(y=0, label=hiddenLabel), fontface="bold") + 
  ylim(-0.5, 0.5) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b

leg <- get_legend(r)

#Plot together
p <- plot_grid(r + theme(legend.position = "none"),
               b, nrow=2, ncol=1, rel_heights=c(1, 0.05 ), align="v", axis='l')
p

p_leg <- plot_grid(p, leg, ncol=2, nrow=1, rel_widths=c(1, 0.15))
p_leg
ggsave(here("Revisions/rin.pdf"), dpi=300, h=4, w=5)
ggsave(here("Revisions/rin.jpg"), dpi=300, h=4, w=5)
