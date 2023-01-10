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

efa <- read.csv(here("CFU/OmniZymoEfaecalisGrowth.csv"), header=TRUE) 
eco <- read.csv(here("CFU/OmniZymoEcoliGrowth.csv"), header=TRUE) 
bth <- read.csv(here("CFU/OmniZymoBthetaGrowth.csv"), header=TRUE) 

condition_palette <- paletteer_d("nbapalettes::pacers_venue")

#Plot B theta growth
t <- ggplot(bth, aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  geom_line(aes(color=Condition), alpha=0.5) +
  theme_bw() + 
  xlab("Time (hours)") +
  ylab("CFU/mL") +
  scale_color_manual(values=condition_palette) +
  theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "cm"))

t

#Plot E coli growth
c <- ggplot(eco, aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  geom_line(aes(color=Condition), alpha=0.4) +
  theme_bw() + 
  xlab("Time (hours)") +
  ylab("CFU/mL") +
  scale_color_manual(values=condition_palette) +
  theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "cm"))

c

#Plot E faecalis growth
f <- ggplot(efa, aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  geom_line(aes(color=Condition), alpha=0.4) +
  theme_bw() + 
  xlab("Time (hours)") +
  ylab("CFU/mL") +
  scale_color_manual(values=condition_palette) +
  theme(text = element_text(size=16), plot.margin = unit(c(0,0,0,0), "cm"))

f

#Plot together
plot_grid(t, c, f, nrow=1, ncol=3, scale=0.9, labels=c("A","B", "C"))

ggsave(here("CFU/growthcurve.pdf"), dpi=300, h=4, w=16)
ggsave(here("CFU/growthcurve.jpg"), dpi=300, h=4, w=16)


# alternative option

# build a named color palette for each condition
condition_palette <- c("#acaaaf","#9a6faa","#4d9222") 
names(condition_palette) <- c("Media", "OMNIgene", "Zymo")

efa <- read.csv(here("CFU/OmniZymoEfaecalisGrowth.csv"), header=TRUE) 
eco <- read.csv(here("CFU/OmniZymoEcoliGrowth.csv"), header=TRUE) 
bth <- read.csv(here("CFU/OmniZymoBthetaGrowth.csv"), header=TRUE) 

efa <- mutate(efa, taxon="E. faecalis")
eco <- mutate(eco, taxon="E. coli")
bth <- mutate(bth, taxon="B. theta")
efa <- mutate(efa, Condition=gsub("E faecalis", "Media", Condition))
eco <- mutate(eco, Condition = gsub("E coli", "Media", Condition))
bth <- mutate(bth, Condition = gsub("B theta", "Media", Condition))

results <- rbind(efa, eco, bth)
results <- mutate(results, Time=as.factor(Time))
ggplot(results, aes(x=Time, y=CFU)) + 
  geom_point(aes(color = Condition, fill = Condition), position=position_dodge(width=0.35), size=2.3) + 
  geom_linerange(aes(x=Time, ymin=0, ymax=CFU, color = Condition), position = position_dodge(width=0.35)) + 
  facet_wrap(~taxon, scales = "free") + 
  theme_bw() +
  ylab("CFU per Microliter") + 
  xlab("Hours") + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  theme(panel.grid.major.x = element_blank(), strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"), 
        strip.background = element_rect(fill="white", color=NA), text = element_text(size=12))


ggsave(here("CFU/growthplot.pdf"), dpi=300, h=2.75, w=8)
ggsave(here("CFU/growthplot.jpg"), dpi=300, h=2.75, w=8)
