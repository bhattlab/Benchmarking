library(ggplot2)
library(tidyverse)
library(here)
library(data.table)
library(paletteer) 

# save a paletteer palette, and subset it to the first ten values
flexxt <- c("#6c7dff", "#ff8c6c")

# add "names" to the colors in the palette, corresponding to the donor names
names(flexxt) <- c("flex", "xt")


windows <- read.table(here("XTvFlex/03_GC/03_windows/concat_windowvals_100k.tsv"), header=FALSE, sep="\t")
names(windows) <- c("Sample", "Contig", "Position", "Coverage", "GC")
cov_vals <- windows %>% group_by(Sample) %>% summarise(mean = mean(Coverage))
windows <- merge(windows, cov_vals)
windows <- mutate(windows, RelCoverage = Coverage/mean)
windows <- mutate(windows, RelCoverage = ifelse(RelCoverage == 0, 0.001, RelCoverage))
windows <- mutate(windows, GCbin = round(GC, digits = 1))
windows <- mutate(windows, GCbinWhole = round(GC, digits = 0))

windows <- mutate(windows, donor=gsub(".*_", "", Sample))
windows <- mutate(windows, kit=gsub("_.*", "", Sample))
windows <- windows %>% select(Sample, donor, kit, RelCoverage, GCbin, GCbinWhole)
windows <- windows %>% mutate(relcovlog2 = log2(RelCoverage))

#windows$GCbin <- as.factor(windows$GCbin)
ggplot(windows %>% filter(GCbin > 20 & GCbin < 80), aes(x=GCbin, y=RelCoverage)) +
  geom_hline(yintercept=1, color="grey") +
  stat_summary(aes(group=kit, color=kit), fun.y=mean, geom="line") + 
  xlim(0,100) + 
  theme_bw() + 
  ylab("Relative Coverage") + 
  xlab("% GC")
  #scale_y_continuous(trans="log2")

# would be nice to add ribbon
ggplot(windows %>% filter(GCbin > 20 & GCbin < 80), aes(x=GCbinWhole, y=RelCoverage)) +
  geom_hline(yintercept=1, color="grey") +
  stat_summary(geom='ribbon', fun.data = mean_cl_normal, fun.args=list(conf.int=0.95), aes(group =kit, fill=kit), alpha = 0.3) +
  stat_summary(aes(group=kit, color=kit), fun.y=mean, geom="line") + 
  xlim(25,75) + 
  theme_bw() + 
  ylab("Relative Coverage") + 
  xlab("% GC") + 
  scale_color_manual(values=flexxt, labels=c("Illumina DNA Prep", "Nextera XT")) + 
  scale_fill_manual(values= flexxt, labels=c("Illumina DNA Prep", "Nextera XT"))

ggsave(here("outputs/figures/XTvFlex_GCrelcoverage.pdf"), width=8, h =5, dpi=300)
ggsave(here("outputs/figures/XTvFlex_GCrelcoverage.jpg"), width=8, h =5, dpi=300)

windows %>% filter(GCbinWhole==50) %>% nrow()
