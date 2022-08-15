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

#####Dry Weight Analysis#####
#Read in dataframe for buffer only weights
buffer <- read.csv(here("qPCR/mai/PreservativeDryWeight.csv"), sep=",", header=TRUE)
colnames(buffer)[colnames(buffer) == "X"] <- "SampleName"
buffer <- buffer %>% separate(SampleName, c("Preservative", "Replicate"), remove=FALSE, sep='_')

#Read in dataframe for stool weights
stool <- read.csv(here("qPCR/mai/DryStoolWeight.csv"), sep=",", header=TRUE)
library(data.table)
stoolt <- transpose(stool)
rownames(stoolt) <- colnames(stool)
colnames(stoolt) <- rownames(stool)
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
stool <-header.true(stoolt)
library(tibble)
stool <- tibble::rownames_to_column(stool, "SampleName")
stool <- stool %>% separate(SampleName, c("Preservative", "Replicate"), remove=FALSE, sep='_')

#Zymo analysis: corrected stool weight with Zymo
#Zymo: 9mL of buffer to 1mL of stool (10% stool)
####FIX - use NF wet weight not zymo wet weight
zcw <- filter(stool, Preservative=="none")
zcw$Corrected_wet <- as.numeric(zcw$Corrected_wet)
zcw <- mutate(zcw, ZWeight=(0.1*Corrected_wet)*0.23)
zcw$ZWeight <- as.numeric(zcw$ZWeight)
meanzstool <- mean(zcw$ZWeight)

#OMNI analysis: corrected stool weight with OMNIgene
#OMNI: 2mL of buffer to 520mg (or 0.52mL) stool (21% stool)
ocw <- filter(stool, Preservative=="none")
ocw$Corrected_wet <- as.numeric(ocw$Corrected_wet)
ocw <- mutate(zcw, OWeight=(0.21*Corrected_wet)*0.23)
ocw$OWeight <- as.numeric(ocw$OWeight)
meanostool <- mean(ocw$OWeight)

