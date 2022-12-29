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

library(paletteer) 
paletteer_d("nbapalettes::pacers_venue")

efa <- read.csv(here("Revisions/OmniZymoEfaecalisGrowth.csv"), header=TRUE) 
eco <- read.csv(here("Revisions/OmniZymoEcoliGrowth.csv"), header=TRUE) 
bth <- read.csv(here("Revisions/OmniZymoBthetaGrowth.csv"), header=TRUE) 

condition_palette <- paletteer_d("nbapalettes::pacers_venue")

#Plot B theta growth
t <- ggplot(bth, aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  geom_line(aes(color=Condition), alpha=0.5) +
  theme_bw() + 
  xlab("Time (hours)") +
  ylab("CFU/mL") +
  ggtitle("B. thetaiotaomicron") +
  scale_color_manual(values=condition_palette) +
  theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "cm"), plot.title=element_text(face="italic"))

t

#Plot E coli growth
c <- ggplot(eco, aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  geom_line(aes(color=Condition), alpha=0.4) +
  theme_bw() + 
  xlab("Time (hours)") +
  ylab("CFU/mL") +
  ggtitle("E. coli") +
  scale_color_manual(values=condition_palette) +
  theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "cm"), plot.title=element_text(face="italic"))

c

#Plot E faecalis growth
f <- ggplot(efa %>% filter(Time==0 | Time==72), aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  geom_line(aes(color=Condition), alpha=0.4) +
  theme_bw() + 
  xlab("Time (hours)") +
  ylab("CFU/mL") +
  scale_color_manual(values=condition_palette) +
  ggtitle("E. faecalis") +
  theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "cm"), plot.title=element_text(face="italic"))

f

#Plot together
plot_grid(t, c, f, nrow=1, ncol=3, scale=0.9)

ggsave(here("Revisions/growthcurve.pdf"), dpi=300, h=4, w=16)
ggsave(here("Revisions/growthcurve.jpg"), dpi=300, h=4, w=16)
