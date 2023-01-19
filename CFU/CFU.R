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
library(ggbeeswarm)

library(paletteer) 
paletteer_d("nbapalettes::pacers_venue")

efa <- read.csv(here("CFU/OmniZymoEfaecalisGrowthReps.csv"), header=TRUE) 
eco <- read.csv(here("CFU/OmniZymoEcoliGrowthReps.csv"), header=TRUE) 
bth <- read.csv(here("CFU/OmniZymoBthetaGrowthReps.csv"), header=TRUE) 

condition_palette <- paletteer_d("nbapalettes::pacers_venue")

#Plot B theta growth
t <- ggplot(bth, aes(x = Time, y=CFU)) + 
  geom_point(width=0.2, aes(color=Condition), shape=16, size=2.5, alpha=0.4) + 
  #geom_line(aes(color=Condition), alpha=0.5) +
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

cfu_styling <-   theme(panel.grid.major.x = element_blank(), text = element_text(size=12))

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
  geom_beeswarm(dodge.width=1, size=2, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("B. theta"))) + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  cfu_styling + 
  theme(legend.position = "bottom", legend.direction = "horizontal")

ec <- ggplot(results %>% filter(taxon == "E. coli"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=2, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("E. coli"))) + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  cfu_styling

ef <- ggplot(results %>% filter(taxon == "E. faecalis"), aes(x=Time, y=CFU, color=Condition)) + 
  geom_beeswarm(dodge.width=1, size=2, cex = 4, method="square") +
  theme_bw() +
  ylab("Million CFU per Milliliter") + 
  xlab("Hours") + 
  ggtitle(expression(~italic("E. faecalis"))) + 
  scale_color_manual(values = condition_palette) + 
  scale_fill_manual(values = condition_palette) +
  cfu_styling

cfuleg <- get_legend(bt)

temp <- plot_grid(bt + theme(legend.position = "none"), 
          ec + theme(axis.title.y = element_blank(), legend.position = "none"), 
          ef + theme(axis.title.y = element_blank(), legend.position = "none"), 
          nrow=1, ncol=3, rel_widths = c(1,1,1))

plot_grid(temp, cfuleg, nrow=2, ncol=1, rel_heights = c(1, 0.1))

ggsave(here("CFU/growthplot.pdf"), dpi=300, h=3, w=8)
ggsave(here("CFU/growthplot.jpg"), dpi=300, h=3, w=8)
