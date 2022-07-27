library(ggplot2)
library(tidyverse)
library(paletteer)
library(ggsignif)
library(here)
library(RColorBrewer)
library(reshape2)
library(vegan)
library(cowplot)
library(ggpubr)
library(ggnewscale)
library(scales)


#####PLate 4 Irregular Amplification#####
qPCR_samplewellp4 <- read.csv(here("qPCR/Plate_4/Plate4_Layout.tsv"), sep="\t", header=TRUE)
qPCR_datap4 <- read.csv(here("qPCR/Plate_4/HotPooqPCR_Plate4.csv"), sep=",", header=TRUE)
colnames(qPCR_samplewellp4)[colnames(qPCR_samplewellp4) == "Well384"] <- "Well"
qPCRplate4 <- merge(qPCR_samplewellp4,qPCR_datap4, by="Well")

#Samples with Irregular mid dilution amplification: D05_ZH_R3 and D09_ZR_R1
qPCRplate4 <- qPCRplate4 %>% separate(SampleName, c("Donor", "Condition", "Replicate", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate4 <- qPCRplate4 %>% unite(Sample, Donor, Condition, Replicate, sep ="_")
qPCRplate4 <- filter(qPCRplate4, PCR_Replicate!="Rep4")
qPCRplate4 <- mutate(qPCRplate4, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate4 <- filter(qPCRplate4, Cq!="Undetermined")
qPCRplate4 <- mutate(qPCRplate4, Cq=as.numeric(Cq))

#D05_ZH_R3
donor5 <- qPCRplate4 %>% filter(Sample=="D05_ZH_R3")
ggplot(donor5, aes(x=DilutionFactor, y=Cq)) +
  geom_point(width=0.1) + 
  geom_smooth(method=lm, se=FALSE) +
  labs(
    x = "Dilution",
    y = "Cq"
  ) +
  theme_bw() +
  scale_x_log10() +
  scale_y_continuous(limits=c(10,25))
ggsave(here("qPCR/mai/D05ZHR3dilutioncq.jpg"), dpi=300, w = 5, h = 5)

#D09_ZR_R1
donor9 <- qPCRplate4 %>% filter(Sample=="D09_ZR_R1")
ggplot(donor9, aes(x=DilutionFactor, y=Cq)) +
  geom_point(width=0.1) + 
  geom_smooth(method=lm, se=FALSE) +
  labs(
    x = "Dilution",
    y = "Cq"
  ) +
  theme_bw() +
  scale_x_log10() +
  scale_y_continuous(limits=c(10,30))
ggsave(here("qPCR/mai/D09ZRR1dilutioncq.jpg"), dpi=300, w = 5, h = 5)
