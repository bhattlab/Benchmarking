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
library(ggpattern)
library(ggpmisc)
library(dplyr)

#####Plate 1: Standard Curve and Absolute Count#####
#Create and edit unified sample datafile
qPCR_samplewellp1 <- read.csv(here("qPCR/2023/Plate1/plate1_384layout.csv"), sep=",", header=TRUE)
qPCR_datap1 <- read.csv(here("qPCR/2023/Plate1/2023qpcrcq_plate1.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp1)[colnames(qPCR_samplewellp1) == "Well384"] <- "Well"
qPCRplate1 <- merge(qPCR_samplewellp1,qPCR_datap1, by="Well")

#Export file
#write.table(qPCRplate1, here("qPCR/2023/Plate1/qPCRplate1.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 1
standard1 <- qPCRplate1 %>% filter(grepl("standard", SampleName))
standard1 <- standard1 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                "DilutionFactor"), remove=FALSE, sep='_')
standard1 <- mutate(standard1, DilutionFactor=as.numeric(DilutionFactor))
standard1 <- mutate(standard1, Cq=as.numeric(Cq))
standard1 <- mutate(standard1, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.692*10^11
copyb <- 1.972*10^11
standard1 <- mutate(standard1, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard1 <- standard1 %>% mutate(logCopyNumber=log10(CopyNumber))
standard1 <- standard1 %>% filter(Cq<22.11) #filter out samples above limit of blank
standard1 <- standard1 %>% filter(DilutionFactor!=1e+02 & DilutionFactor!=1e+08)

stdp1model <- lm(Cq~logCopyNumber, data=standard1)
stdplate1 <- ggplot(data = standard1, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 22.11, linetype = "dotted", color = "black") +
  annotate("text", x = max(standard1$logCopyNumber), y = 24, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7, y = 20, label = "y = 36.3 - 3.56x \nR^2=0.99", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 1")

stdplate1

#ggsave(here("qPCR/2023/outputs/stdcurve_plate1.jpg"), dpi=300, w = 5, h = 5)

#Edit dataframe
qPCRplate1 <- read.csv(here("qPCR/2023/Plate1/qPCRplate1.csv"), sep=",", header=TRUE)
qPCRplate1 <- qPCRplate1 %>% filter(!grepl("standard", SampleName))
qPCRplate1 <- mutate(qPCRplate1, Cq=as.numeric(Cq))
lobplate1 <- qPCRplate1 %>% filter(Cq>21.77) #max Cq standard curve

#Analysis of technical qPCR replicates
qPCRplate1 <- qPCRplate1 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
techrepnf1 <- ggplot(qPCRplate1 %>% filter(Condition=="NF"), aes(x=Donor, y=Cq, color=PCR_Replicate)) +
  geom_jitter(width=0.1) +
  labs(
    x = "Sample",
    y = "Cq"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

techrepnf1

#ggsave(here("qPCR/2023/outputs/techrep_p1nf.jpg"), dpi=300, w = 6, h = 5)
techrepp1 <- spread(qPCRplate1 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)
techrepp1 <- mutate(techrepp1, Cqdiff=(abs(Rep2-Rep1)))

spreadplate1 <- ggplot(techrepp1, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Cq 1",
    y = "Cq 2"
  ) +
  ggtitle("Plate 1")
spreadplate1

#ggsave(here("qPCR/2023/outputs/techrep_p1.jpg"), dpi=300, w = 5, h = 5)

#Calculation using combined standard curve model
#Standard curve is: y=36.3-3.56x
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=(Cq-36.3)/-3.56)
qPCRplate1 <- qPCRplate1 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate1 <- qPCRplate1 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
qPCRplate1 <- mutate(qPCRplate1, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
qPCRplate1 <- mutate(qPCRplate1, label="none")
qPCRplate1 <- select(qPCRplate1,-SampleType & -Tm1 & -Tm2 & -Tm3 & -Tm4)
published <- read.csv(here("qPCR/Plate_1/qPCR_plate1.csv"), sep=",", header=TRUE)
published_colorder <- names(published)
qPCRplate1 <- qPCRplate1 %>% select(all_of(published_colorder))
#Filter out samples that are re-run on plate 4 or 5
qPCRplate1 <- filter(qPCRplate1, SampleName!="D07_ZF_R1", SampleName!="D03_NF_R1",
                     SampleName!="D10_OF_R1", SampleName!="D04_ZF_R1", SampleName!="D01_ZH_R1",
                     SampleName!="D06_NF_R1", SampleName!="D01_ZF_R1", SampleName!="D03_ZR_R1",
                     SampleName!="D10_ZR_R1", SampleName!="D07_ZH_R1", SampleName!="D04_ZH_R1",
                     SampleName!="D04_OH_R1", SampleName!="D03_ZH_R1", SampleName!="D05_ZH_R1",
                     SampleName!="D06_ZH_R1", SampleName!="D07_ZR_R1") 
#Export file
#write.table(qPCRplate1, here("qPCR/2023/Plate1/qPCRplate1_copynu.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Ordering conditions
qPCRplate1 <- mutate(qPCRplate1, PlotOrder=ifelse(Condition == "NF", 1, 
                                                        ifelse(Condition == "OF", 2, 
                                                               ifelse(Condition == "OR", 3, 
                                                                      ifelse(Condition == "OH", 4, 
                                                                             ifelse(Condition == "ZF", 5,
                                                                                    ifelse(Condition == "ZR", 6, 7)))))))
#Color palette for conditions
condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plotting copies/uL by condition
ggplot(qPCRplate1 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=CopyNumber, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/uL"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+14), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")

ggsave(here("qPCR/2023/outputs/plate1copiesuL.jpg"), dpi=300, w = 5, h = 5)

#Copies/g Absolute Count
massnf <- (50/0.062)
massomni <- (50/0.037)
masszymo <- (50/0.08)

qPCRplate1 <- mutate(qPCRplate1, Preservative=substr(Condition,1,1))
qPCRplate1 <- mutate(qPCRplate1, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                                   ifelse(Preservative=="O", CopyNumber*massomni,
                                                          ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                                          )))

ggplot(qPCRplate1 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=Copiesgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+18), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")
ggsave(here("qPCR/2023/outputs/plate1copiesg.jpg"), dpi=300, w = 5, h = 5)

#Export final plate 1 datafile
#write.table(qPCRplate1, here("qPCR/Plate_1/qPCR_plate1_final.csv"), row.names = FALSE, sep=",", quote=FALSE)  

######Plate 2: Standard Curve and Absolute Count#####
#Create and edit unified sample datafile
qPCR_samplewellp2 <- read.csv(here("qPCR/2023/Plate2/plate2_384layout.csv"), sep=",", header=TRUE)
qPCR_datap2 <- read.csv(here("qPCR/2023/Plate2/2023qpcrcq_plate2.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp2)[colnames(qPCR_samplewellp2) == "Well384"] <- "Well"
qPCRplate2 <- merge(qPCR_samplewellp2,qPCR_datap2, by="Well")

#Export file
#write.table(qPCRplate2, here("qPCR/2023/Plate2/qPCRplate2.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 2
standard2 <- qPCRplate2 %>% filter(grepl("standard", SampleName))
standard2 <- standard2 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                "DilutionFactor"), remove=FALSE, sep='_')
standard2 <- mutate(standard2, DilutionFactor=as.numeric(DilutionFactor))
standard2 <- mutate(standard2, Cq=as.numeric(Cq))
standard2 <- mutate(standard2, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.692*10^11
copyb <- 1.972*10^11
standard2 <- mutate(standard2, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard2 <- standard2 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
standard2 <- standard2 %>% filter(Cq<22.83) #filter out samples above limit of blank
standard2 <- standard2 %>% filter(DilutionFactor!=1e+02 & DilutionFactor!=1e+08)


stdp2model <- lm(Cq~logCopyNumber, data=standard2)
stdplate2 <- ggplot(data = standard2, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 22.83, linetype = "dotted", color = "black") +
  annotate("text", x = max(standard2$logCopyNumber), y = 24, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7.5, y = 20, label = "y = 36.4 - 3.59x \nR^2=0.98", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 2")

stdplate2

#ggsave(here("qPCR/2023/outputs/stdcurve_plate2.jpg"), dpi=300, w = 5, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate2 <- read.csv(here("qPCR/2023/Plate2/qPCRplate2.csv"), sep=",", header = TRUE)
qPCRplate2 <- qPCRplate2 %>% filter(!grepl("standard", SampleName))
qPCRplate2 <- mutate(qPCRplate2, Cq=as.numeric(Cq))
lobplate2 <- qPCRplate2 %>% filter(Cq>22.19) #max Cq standard curve

#Analysis of technical qPCR replicates
qPCRplate2 <- qPCRplate2 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
techrep2 <- ggplot(qPCRplate2 %>% filter(Condition=="NF"), aes(x=Donor, y=Cq, color=PCR_Replicate)) +
  geom_jitter(width=0.1) +
  labs(
    x = "Sample",
    y = "Cq"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

techrep2
ggsave(here("qPCR/2023/outputs/techrep_p2nf.jpg"), dpi=300, w = 6, h = 5)

techrepp2 <- spread(qPCRplate2 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)
techrepp2 <- mutate(techrepp2, Cqdiff=(abs(Rep2-Rep1)))

spreadplate2 <- ggplot(techrepp2, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Cq 1",
    y = "Cq 2"
  )  +
  ggtitle("Plate 2")

spreadplate2

ggsave(here("qPCR/2023/outputs/techrep_p2.jpg"), dpi=300, w = 5, h = 5)

#Calculation using combined standard curve model
#Standard curve is: y=36.4-3.59x
qPCRplate2 <- qPCRplate2 %>% mutate(logCopyNumber=(Cq-36.4)/-3.59)
qPCRplate2 <- qPCRplate2 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate2 <- qPCRplate2 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
qPCRplate2 <- mutate(qPCRplate2, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
qPCRplate2 <- select(qPCRplate2,-SampleType & -Tm1 & -Tm2 & -Tm3 & -Tm4)
#Filter out samples that are re-run on plate 4 or 5
qPCRplate2 <- filter(qPCRplate2, SampleName!="D09_ZR_R2", SampleName!="D05_ZF_R2",
                     SampleName!="D08_ZH_R2", SampleName!="D07_ZR_R2")
#Export file
#write.table(qPCRplate2, here("qPCR/2023/Plate2/qPCRplate2_copynu.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Ordering conditions
qPCRplate2 <- mutate(qPCRplate2, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
#Color palette for conditions
condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plotting copies/uL by condition
ggplot(qPCRplate2 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=CopyNumber, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(
    x = "Condition",
    y = "Copies/uL"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5,1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif")

ggsave(here("qPCR/2023/outputs/plate2copiesuL.jpg"), dpi=300, w = 5, h = 5)

#Copies/g Absolute Count
massnf <- (50/0.062)
massomni <- (50/0.037)
masszymo <- (50/0.08)

qPCRplate2 <- mutate(qPCRplate2, Preservative=substr(Condition,1,1))
qPCRplate2 <- mutate(qPCRplate2, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                                   ifelse(Preservative=="O", CopyNumber*massomni,
                                                          ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                                   )))

ggplot(qPCRplate2 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=Copiesgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+18), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")
ggsave(here("qPCR/2023/outputs/plate2copiesg.jpg"), dpi=300, w = 5, h = 5)

##### Plate 3: Standard Curve and Absolute Count #####
#Create and edit unified sample datafile
qPCR_samplewellp3 <- read.csv(here("qPCR/2023/Plate3/plate3_384layout.csv"), sep=",", header=TRUE)
qPCR_datap3 <- read.csv(here("qPCR/2023/Plate3/2024qpcrcq_plate3_r3.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp3)[colnames(qPCR_samplewellp3) == "Well384"] <- "Well"
qPCRplate3 <- merge(qPCR_samplewellp3,qPCR_datap3, by="Well")

#Export file
#write.table(qPCRplate3, here("qPCR/2023/Plate3/qPCRplate3_r3.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 3
standard3 <- qPCRplate3 %>% filter(grepl("standard", SampleName))
standard3 <- standard3 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                "DilutionFactor"), remove=FALSE, sep='_')
standard3 <- mutate(standard3, DilutionFactor=as.numeric(DilutionFactor))
standard3 <- mutate(standard3, Cq=as.numeric(Cq))
standard3 <- mutate(standard3, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.692*10^11
copyb <- 1.972*10^11
standard3 <- mutate(standard3, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard3 <- standard3 %>% mutate(logCopyNumber=log10(CopyNumber))

#Filter out standards above limit of blank
standard3 <- standard3 %>% filter(Cq<25.16)
standard3 <- filter(standard3, DilutionFactor!=1e+02 & DilutionFactor!=1e+03 & DilutionFactor!=1e+08)

stdp3model <- lm(Cq~logCopyNumber, data=standard3)
stdplate3 <- ggplot(data = standard3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 25.16, linetype = "dotted", color = "black") +
  annotate("text", x = max(standard3$logCopyNumber), y = 24, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7, y = 20, label = "y = 37.8 - 3.69x \nR^2=0.93", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 3")

stdplate3

#ggsave(here("qPCR/2023/outputs/stdcurve_plate3_lob_r2.jpg"), dpi=300, w = 6, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate3 <- read.csv(here("qPCR/2023/Plate3/qPCRplate3_r3.csv"), sep=",", header = TRUE)
qPCRplate3 <- qPCRplate3 %>% filter(!grepl("standard", SampleName))
qPCRplate3 <- mutate(qPCRplate3, Cq=as.numeric(Cq))
lobplate3 <- qPCRplate3 %>% filter(Cq>24.83) #max Cq standard curve
qPCRplate3 <- qPCRplate3 %>% filter(!grepl("Buffer", SampleName))
qPCRplate3 <- qPCRplate3 %>% filter(!grepl("Water", SampleName))
qPCRplate3 <- qPCRplate3 %>% filter(!grepl("CO", SampleName))

#Analysis of technical qPCR replicates
qPCRplate3 <- qPCRplate3 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
techrep3 <- ggplot(qPCRplate3 %>% filter(Condition=="NF"), aes(x=Donor, y=Cq, color=PCR_Replicate)) +
  geom_jitter(width=0.1) +
  labs(
    x = "Sample",
    y = "Cq"
  ) +
  theme_bw() +
  theme(legend.position = "none") 

techrep3
ggsave(here("qPCR/2023/outputs/techrep_p3nf.jpg"), dpi=300, w = 6, h = 5)

techrepp3 <- spread(qPCRplate3 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)
techrepp3 <- mutate(techrepp3, Cqdiff=(abs(Rep2-Rep1)))
rerunp3samples <- filter(techrepp3, Cqdiff>3)

spreadplate3 <- ggplot(techrepp3, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Cq 1",
    y = "Cq 2"
  )  +
  ggtitle("Plate 3")

spreadplate3

ggsave(here("qPCR/2023/outputs/techrep_p3.jpg"), dpi=300, w = 5, h = 5)

#Calculation using combined standard curve model
#Standard curve is: y=37.8-3.69x
qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=(Cq-37.8)/-3.69)
qPCRplate3 <- qPCRplate3 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate3 <- qPCRplate3 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
qPCRplate3 <- mutate(qPCRplate3, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
qPCRplate3 <- select(qPCRplate3,-SampleType & -Tm1 & -Tm2 & -Tm3 & -Tm4)
#Filter out samples that are re-run on plate 4 or 5
qPCRplate3 <- filter(qPCRplate3, SampleName!="D05_ZR_R3", SampleName!="D04_ZH_R3",
                     SampleName!="D05_ZH_R3")
#Export file
#write.table(qPCRplate3, here("qPCR/2023/Plate3/qPCRplate3_copynu.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Ordering conditions
qPCRplate3 <- mutate(qPCRplate3, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
#Color palette for conditions
condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plotting copies/uL by condition
ggplot(qPCRplate3 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=CopyNumber, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/uL"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

ggsave(here("qPCR/2023/outputs/plate3copiesuL.jpg"), dpi=300, w = 5, h = 5)

#Copies/g Absolute Count
massnf <- (50/0.062)
massomni <- (50/0.037)
masszymo <- (50/0.08)

qPCRplate3 <- mutate(qPCRplate3, Preservative=substr(Condition,1,1))
qPCRplate3 <- mutate(qPCRplate3, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                                   ifelse(Preservative=="O", CopyNumber*massomni,
                                                          ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                                   )))

ggplot(qPCRplate3 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=Copiesgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+18), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")
ggsave(here("qPCR/2023/outputs/plate3copiesg.jpg"), dpi=300, w = 5, h = 5)

######Plate 4: Standard Curve and Absolute Count of missing Cq samples#####
#Create and edit unified sample datafile
qPCR_samplewellp4 <- read.csv(here("qPCR/2023/Plate4/plate4_384layout.csv"), sep=",", header=TRUE)
qPCR_datap4 <- read.csv(here("qPCR/2023/Plate4/2024_qpcrcq_plate4.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp4)[colnames(qPCR_samplewellp4) == "Well384"] <- "Well"
qPCRplate4 <- merge(qPCR_samplewellp4,qPCR_datap4, by="Well")

#Export file
#write.table(qPCRplate4, here("qPCR/2023/Plate4/qPCRplate4.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 4
standard4 <- qPCRplate4 %>% filter(grepl("standard", SampleName))
standard4 <- standard4 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                "DilutionFactor"), remove=FALSE, sep='_')
standard4 <- mutate(standard4, DilutionFactor=as.numeric(DilutionFactor))
standard4 <- mutate(standard4, Cq=as.numeric(Cq))
standard4 <- mutate(standard4, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
#Note: using a new stock of standards
copyf <- 1.077*10^11
copyb <- 1.071*10^11
standard4 <- mutate(standard4, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard4 <- standard4 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
standard4 <- standard4 %>% filter(Cq<27.5) #filter out samples above limit of blank
standard4 <- standard4 %>% filter(DilutionFactor!=1e+02
                                  & DilutionFactor!=1e+08)

stdp4model <- lm(Cq~logCopyNumber, data=standard4)
stdplate4 <- ggplot(data = standard4, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 27.5, linetype = "dotted", color = "black") +
  annotate("text", x = max(standard4$logCopyNumber), y = 27, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 6.5, y = 20, label = "y = 35.9 - 3.63x \nR^2=0.99", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 4")

stdplate4

#ggsave(here("qPCR/2023/outputs/stdcurve_plate4.jpg"), dpi=300, w = 5, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate4 <- read.csv(here("qPCR/2023/Plate4/qPCRplate4.csv"), sep=",", header = TRUE)
qPCRplate4 <- qPCRplate4 %>% filter(!grepl("standard", SampleName))
qPCRplate4 <- mutate(qPCRplate4, Cq=as.numeric(Cq))
lobplate4 <- qPCRplate4 %>% filter(Cq>22.01) #max Cq standard curve
qPCRplate4 <- qPCRplate4 %>% filter(!grepl("Buffer", SampleName))
qPCRplate4 <- qPCRplate4 %>% filter(!grepl("Water", SampleName))
test <- qPCRplate4[!is.na(qPCRplate4$Cq), ] #determine which samples are missing Cq values
#Remove samples missing one or more Cq values
qPCRplate4 <- qPCRplate4 %>% filter(SampleName!= "D03_ZR_R1", SampleName!="D04_ZF_R1", 
  SampleName!="D04_ZH_R1", SampleName!="D04_ZH_R3", SampleName!="D05_ZH_R1", 
  SampleName!="D07_ZR_R1", SampleName!="D07_ZR_R2", SampleName!="D09_ZR_R2") 

#Analysis of technical qPCR replicates
qPCRplate4 <- qPCRplate4 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
techrepp4 <- spread(qPCRplate4 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)
techrepp4 <- mutate(techrepp4, Cqdiff=(Rep2-Rep1))

spreadplate4 <- ggplot(techrepp4, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Cq 1",
    y = "Cq 2"
  )  +
  ggtitle("Plate 4")

spreadplate4


ggsave(here("qPCR/2023/outputs/techrep_p4.jpg"), dpi=300, w = 5, h = 5)

techrows <- ggplot(techrepp4, aes(x=SampleName, y=Cqdiff)) +
  geom_point(size=0.9) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "grey") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y= expression(paste("ΔCq")))
techrows
ggsave(here("qPCR/2023/outputs/techreprows_p4.jpg"), dpi=300, w = 8, h = 4)

#Calculation using combined standard curve model
#Standard curve is: y = 35.9 - 3.63x
qPCRplate4 <- qPCRplate4 %>% mutate(logCopyNumber=(Cq-35.9)/-3.63)
qPCRplate4 <- qPCRplate4 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate4 <- qPCRplate4 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
qPCRplate4 <- mutate(qPCRplate4, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
qPCRplate4 <- select(qPCRplate4,-SampleType & -Tm1 & -Tm2 & -Tm3 & -Tm4)

#Export file
#write.table(qPCRplate4, here("qPCR/2023/Plate4/qPCRplate4_copynu.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Ordering conditions
qPCRplate4 <- mutate(qPCRplate4, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
#Color palette for conditions
condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plotting copies/uL by condition
ggplot(qPCRplate4 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=CopyNumber, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(
    x = "Condition",
    y = "Copies/uL"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5,1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif")

#ggsave(here("qPCR/2023/outputs/plate4copiesuL.jpg"), dpi=300, w = 5, h = 5)

#Copies/g Absolute Count
massnf <- (50/0.062)
massomni <- (50/0.037)
masszymo <- (50/0.08)

qPCRplate4 <- mutate(qPCRplate4, Preservative=substr(Condition,1,1))
qPCRplate4 <- mutate(qPCRplate4, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                                   ifelse(Preservative=="O", CopyNumber*massomni,
                                                          ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                                   )))

ggplot(qPCRplate4 %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=Copiesgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+18), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")
#ggsave(here("qPCR/2023/outputs/plate4copiesg.jpg"), dpi=300, w = 5, h = 5)

#####Plate 5: Standard Curve#####
#Create and edit unified sample datafile
qPCR_samplewellp5 <- read.csv(here("qPCR/2023/Plate5/plate5_384layout.csv"), sep=",", header=TRUE)
qPCR_datap5 <- read.csv(here("qPCR/2023/Plate5/2024_qpcrcq_plate5.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp5)[colnames(qPCR_samplewellp5) == "Well384"] <- "Well"
qPCRplate5 <- merge(qPCR_samplewellp5,qPCR_datap5, by="Well")

#Export file
#write.table(qPCRplate5, here("qPCR/2023/Plate5/qPCRplate5.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 5
standard5 <- qPCRplate5 %>% filter(grepl("standard", SampleName))
standard5 <- standard5 %>% separate(SampleName, c("Standard", "StandardReplicate", 
                                                "DilutionFactor"), remove=FALSE, sep='_')
standard5 <- mutate(standard5, DilutionFactor=as.numeric(DilutionFactor))
standard5 <- mutate(standard5, Cq=as.numeric(Cq))
standard5 <- mutate(standard5, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
#NOTE: new stock of standards used
copyf <- 1.077*10^11
copyb <- 1.071*10^11
standard5 <- mutate(standard5, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Color palette for standards
standard_palette <- c("#e28743","#76b5c5")
names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#Standard curve using both B vulgatus and F prausnitizii data
standard5 <- standard5 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
standard5 <- standard5 %>% filter(Cq<28.98) #filter out samples above limit of blank
standard5 <- standard5 %>% filter(DilutionFactor!=1e+08 & DilutionFactor!=1e+02)

stdp5model <- lm(Cq~logCopyNumber, data=standard5)
stdplate5 <- ggplot(data = standard5, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x, color="Grey") + 
  geom_jitter(aes(color=Standard, shape=StandardReplicate), width=0.1) +
  theme_bw() +
  geom_hline(yintercept = 28.98, linetype = "dotted", color = "black") +
  annotate("text", x = max(standard5$logCopyNumber), y = 29, label = "LoB", color = "black", hjust=2, vjust=1.5) +
  scale_color_manual(values=standard_palette) +
  annotate("text", x = 7.5, y = 20, label = "y = 35.6 - 3.59x \nR^2=0.99", size=3) +
  theme(legend.position="none") +
  ggtitle("Plate 5")

stdplate5

#ggsave(here("qPCR/2023/outputs/stdcurve_plate5.jpg"), dpi=300, w = 5, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate5 <- read.csv(here("qPCR/2023/Plate5/qPCRplate5.csv"), sep=",", header = TRUE)
qPCRplate5 <- qPCRplate5 %>% filter(!grepl("standard", SampleName))
qPCRplate5 <- mutate(qPCRplate5, Cq=as.numeric(Cq))
lobplate5 <- qPCRplate5 %>% filter(Cq>22.33) #max Cq standard curve
qPCRplate5 <- qPCRplate5 %>% filter(!grepl("Buffer", SampleName))
qPCRplate5 <- qPCRplate5 %>% filter(!grepl("Water", SampleName))
#Remove overlapping samples
qPCRplate5 <- qPCRplate5 %>% filter(SampleName!="D07_ZH_R1", SampleName!="D10_ZR_R1")

#Analysis of technical qPCR replicates
qPCRplate5 <- qPCRplate5 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
techrepp5 <- spread(qPCRplate5 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)
techrepp5 <- mutate(techrepp5, Cqdiff=(Rep2-Rep1))

spreadplate5 <- ggplot(techrepp5, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Cq 1",
    y = "Cq 2"
  )  +
  ggtitle("Plate 5")

spreadplate5


#ggsave(here("qPCR/2023/outputs/techrep_p5.jpg"), dpi=300, w = 5, h = 5)

techrows <- ggplot(techrepp5, aes(x=SampleName, y=Cqdiff)) +
  geom_point(size=0.9) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "grey") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y= expression(paste("ΔCq")))
techrows
#ggsave(here("qPCR/2023/outputs/techreprows_p5.jpg"), dpi=300, w = 8, h = 4)
#Plot together the tech rep Cq comparison across the five plates
#Plot Cq differential for all plates
techrepp1 <- mutate(techrepp1, PlateNu=1)
techrepp2 <- mutate(techrepp2, PlateNu=2)
techrepp3 <- mutate(techrepp3, PlateNu=3)
techrepp4 <- mutate(techrepp4, PlateNu=4)
techrepp5 <- mutate(techrepp5, PlateNu=5)
allplatetr <- rbind(techrepp1, techrepp2, techrepp3, techrepp4, techrepp5)

allplatetr <- allplatetr %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
allplatetr <- unite(allplatetr, DonorCondition, c(Donor, Condition))
allplatetr <- allplatetr %>% filter(!grepl("Water", SampleName) 
                                    & !grepl("Buffer", SampleName) 
                                    & !grepl("PCO", SampleName)
                                    & !grepl("NCO", SampleName))
allplatetr$PlateNu <- as.factor(allplatetr$PlateNu)
allplatetr$Cqdiff <- abs(allplatetr$Cqdiff)

allplatecqdiff <- ggplot(allplatetr, aes(x=DonorCondition, y=Cqdiff, color=PlateNu)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 3, linetype = "dashed", color = "black") +
  scale_color_paletteer_d("rcartocolor::Temps") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() 
  ) +
  labs(
    x= "Sample",
    y = expression(Delta~"Cq"),
    color = "Plate")
allplatecqdiff

ggsave(here("qPCR/2023/outputs/cqdiff_allplate.jpg"), dpi=300, w = 9, h = 3.5)

#Calculation using combined standard curve model
#Standard curve is: y = 35.6 - 3.59x
qPCRplate5 <- qPCRplate5 %>% mutate(logCopyNumber=(Cq-35.6)/-3.59)
qPCRplate5 <- qPCRplate5 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate5 <- qPCRplate5 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')
qPCRplate5 <- mutate(qPCRplate5, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
qPCRplate5 <- select(qPCRplate5,-SampleType & -Tm1 & -Tm2 & -Tm3 & -Tm4)
#Filter out samples missing 1 or more Cq values
qPCRplate5 <- filter(qPCRplate5, SampleName!="D07_ZR_R2", SampleName!="D05_ZH_R1",
                     SampleName!="D04_ZH_R3", SampleName!="D04_ZH_R1", SampleName!="D04_ZF_R1")
#Export file
#write.table(qPCRplate5, here("qPCR/2023/Plate5/qPCRplate5_copynu.csv"), row.names = FALSE, sep=",", quote=FALSE)  

##### All plate analysis: Standards and absolute count#####
#Cowplot all standards
#Make legend - to include into cowplot without warping plots
stdlegend <- ggplot(standard5, aes(x=logCopyNumber, y=Cq)) +
  geom_jitter(width=0.1, aes(color=Standard, shape=StandardReplicate)) +
  labs(
    x = "Sample",
    y = "Cq"
  ) +
  theme_bw() +
  scale_color_manual(values=standard_palette)
stdlegend
legend <- get_legend(stdlegend)
legend

allstd <- plot_grid(stdplate1, stdplate2, stdplate3, 
                    stdplate4, stdplate5,legend, nrow=2, ncol=3, 
                    rel_widths=c(1, 1, 1, 1, 1, 0.4))
allstd
ggsave(here("qPCR/2023/outputs/allstdcurves.jpg"), dpi=300, w=10, h=7)
#Compare standards from 2023 run to 2022 publication run
#Create dataframe with linear regressions
x <- seq(0, 15, by = 1)
data <- data.frame(x = x, 
                   y23plate1 = 36.3 - 3.56 * x, 
                   y23plate2 = 36.4 - 3.59 * x, 
                   y23plate3 = 37.8 - 3.69 * x,
                   y23plate4 = 35.9 - 3.63 * x,
                   y23plate5 = 35.6 - 3.59 * x,
                   y22plate1 = 35 - 2.93 * x,
                   y22plate2 = 35.1 - 3.01 * x,
                   y22plate3 = 35.5 - 3.13 * x)


# Create plot with 22 and 23 qPCR standard curves
ggplot(data, aes(x = x)) +
  geom_line(aes(y = y23plate1, color = "'23 Plate 1")) +
  geom_line(aes(y = y23plate2, color = "'23 Plate 2")) +
  geom_line(aes(y = y23plate3, color = "'23 Plate 3")) +
  geom_line(aes(y = y23plate4, color = "'23 Plate 4")) +
  geom_line(aes(y = y23plate5, color = "'23 Plate 5")) +
  labs(x = "logCopyNumber", y = "Cq", color="") +
  xlim(0,15) +
  ylim(-5, 30) +
  theme_bw() 
  #scale_color_paletteer_d("colorBlindness::Brown2Blue10Steps")

#ggsave(here("qPCR/2023/outputs/stdregression.jpg"), dpi=300, w = 6, h = 5)

#Comparison of standard Cq values across plates
allstandard <- rbind(standard1, standard2, standard3, standard4, standard5)
allstandard <- allstandard %>% mutate(Plate=as.character(Plate))
allstandard <- allstandard %>% mutate(DilutionFactor=as.character(DilutionFactor))
standardcq <- ggplot(allstandard, aes(x=DilutionFactor, y=Cq, color=Plate)) +
  geom_point(size=0.9) +
  theme_bw() 
standardcq
#ggsave(here("qPCR/2023/outputs/standardcq.jpg"), dpi=300, w = 6, h = 5)
#Calculate Cq difference between each dilution across plates
s1 <- standard1 %>%
  group_by(DilutionFactor) %>%
  summarise_at(vars(Cq), list(Cq=mean))
s1 <- s1 %>%
  mutate(deltaCq = Cq - lag(Cq), Plate=1)
s2 <- standard2 %>%
  group_by(DilutionFactor) %>%
  summarise_at(vars(Cq), list(Cq=mean))
s2 <- s2 %>%
  mutate(deltaCq = Cq - lag(Cq), Plate=2)
s3 <- standard3 %>%
  group_by(DilutionFactor) %>%
  summarise_at(vars(Cq), list(Cq=mean))
s3 <- s3 %>%
  mutate(deltaCq = Cq - lag(Cq), Plate=3)
s4 <- standard4 %>%
  group_by(DilutionFactor) %>%
  summarise_at(vars(Cq), list(Cq=mean))
s4 <- s4 %>%
  mutate(deltaCq = Cq - lag(Cq), Plate=4)
s5 <- standard5 %>%
  group_by(DilutionFactor) %>%
  summarise_at(vars(Cq), list(Cq=mean))
s5 <- s5 %>%
  mutate(deltaCq = Cq - lag(Cq), Plate=5)
allstandarddcq <- rbind(s1, s2, s3, s4, s5)
allstandarddcq <- allstandarddcq %>% mutate(DilutionFactor=as.character(DilutionFactor),
                                            Plate=as.character(Plate))
dstandardcq <- ggplot(allstandarddcq, aes(x=DilutionFactor, y=deltaCq, color=Plate)) +
  geom_point(size=0.9) +
  theme_bw() +
  ylim(2, 5)
dstandardcq
#ggsave(here("qPCR/2023/outputs/standarddeltacq.jpg"), dpi=300, w = 6, h = 5)

#Color palette for conditions
regression_palette <- c("#0D47A1FF", "#2166ACFF","#4393C3FF","#92C5DEFF","#D1E5F0FF","#F4A582FF", "#D6604DFF")
names(regression_palette) <- c("'23 Plate 1", "'23 Plate 2", "'23 Plate 3", 
                               "'23 Plate 4", "'23 Plate 5", "FOS Plate 1", "FOS Plate 2")

#Linear regression standard comparison with FOS and 2023 benchmarking
x <- seq(0, 15, by = 1)
data <- data.frame(x = x, 
                   yfosplate1 = 34.6 - 3.34 * x,
                   yfosplate2 = 35.1 - 3.48 * x,
                   y23plate1 = 36.3 - 3.56 * x, 
                   y23plate2 = 36.4 - 3.59 * x, 
                   y23plate3 = 37.8 - 3.69 * x,
                   y23plate4 = 35.9 - 3.63 * x,
                   y23plate5 = 35.6 - 3.59 * x)

ggplot(data, aes(x = x)) +
  geom_line(aes(y = y23plate1, color = "'23 Plate 1")) +
  geom_line(aes(y = y23plate2, color = "'23 Plate 2")) +
  geom_line(aes(y = y23plate3, color = "'23 Plate 3")) +
  geom_line(aes(y = y23plate4, color = "'23 Plate 4")) +
  geom_line(aes(y = y23plate5, color = "'23 Plate 5")) +
  geom_line(aes(y = yfosplate1, color = "FOS Plate 1")) +
  geom_line(aes(y = yfosplate2, color = "FOS Plate 2")) +
  labs(x = "logCopyNumber", y = "Cq", color="") +
  xlim(0,15) +
  ylim(0, 30) +
  theme_bw() +
  scale_color_manual(values=regression_palette)

#ggsave(here("qPCR/2023/outputs/stdcurves23fos.jpg"), dpi=300, w = 6, h = 5)

#Copies/uL of all plates - combined data
qPCRplate1 <- read.csv(here("qPCR/2023/Plate1/qPCRplate1_copynu.csv"), sep=",", header=TRUE)
qPCRplate1 <- select(qPCRplate1, -label)
qPCRplate2 <- read.csv(here("qPCR/2023/Plate2/qPCRplate2_copynu.csv"), sep=",", header=TRUE)
qPCRplate3 <- read.csv(here("qPCR/2023/Plate3/qPCRplate3_copynu.csv"), sep=",", header=TRUE)
qPCRplate4 <- read.csv(here("qPCR/2023/Plate4/qPCRplate4_copynu.csv"), sep=",", header=TRUE)
qPCRplate5 <- read.csv(here("qPCR/2023/Plate5/qPCRplate5_copynu.csv"), sep=",", header=TRUE)
allcopies <- rbind(qPCRplate1, qPCRplate2, qPCRplate3, qPCRplate4, qPCRplate5)
qPCRplate1 <- qPCRplate1 %>% filter(!grepl("standard", SampleName))
allcopies <- filter(allcopies, !grepl("PCO", SampleName) & !grepl("Buffer", SampleName)
                    & !grepl("NCO", SampleName) & !grepl("Water", SampleName))
#210 samples run by qPCR, 5 samples unsuccessful qPCR-> 410 reactions is seen in this dataframe
#write.table(allcopies, here("qPCR/2023/allsamplecopiesul.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Color palette for conditions
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")
#Ordering conditions
allcopies <- mutate(allcopies, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))

allcopiesul <- ggplot(allcopies %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=CopyNumber, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/uL"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 
allcopiesul

#ggsave(here("qPCR/2023/outputs/allplatecopiesuL.jpg"), dpi=300, w = 5, h = 5)
#Calculate copies/gram for each sample
massnf <- (50/0.062)
massomni <- (50/0.013)
masszymo <- (50/0.006)

allcopies <- mutate(allcopies, Preservative=substr(Condition,1,1))
allcopies <- mutate(allcopies, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                               ifelse(Preservative=="O", CopyNumber*massomni,
                                                      ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                               )))


ggplot(allcopies %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=Copiesgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Copies/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+9, 1e+19), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")
#ggsave(here("qPCR/2023/outputs/allplatecopiesg.jpg"), dpi=300, w = 5, h = 5)

##### Absolute Counts at Phylum Level#####
#Create a unified dataframe containing all Plate data
allcopies <- read.csv(here("qPCR/2023/allsamplecopiesul.csv"), sep=",", header=TRUE)

allplate <- allcopies %>% select(SampleName,CopyNumber, SampleNameFull, Condition, Cq)
#Calculate copies/gram for each sample
massnf <- (50/0.062)
massomni <- (50/0.013)
masszymo <- (50/0.006)

allplate <- mutate(allplate, Preservative=substr(Condition,1,1))
allplate <- mutate(allplate, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                                   ifelse(Preservative=="O", CopyNumber*massomni,
                                                          ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                                   )))

#Read in kraken outputs and rrnDB copy number data
krakenphyla <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_phylum_percentage.txt"), sep="\t", header=TRUE)
copyndb <- read.csv(here("qPCR/phylum_copynumber_table.tsv"), sep="\t", header=TRUE)

krakenphyla$Phylum <- row.names(krakenphyla)
row.names(krakenphyla) <- NULL
kraken_long <- melt(krakenphyla, id.vars = "Phylum", variable.name = "Sample", value.name = "rel_abundance")
kraken_long <- mutate(kraken_long, Sample=gsub("\\.", "_", Sample))
colnames(copyndb)[colnames(copyndb) == "name"] <- "Phylum"
kraken_long <- merge(kraken_long, copyndb, by="Phylum", all.x=TRUE)

#Calculating total bacteria per sample (median and corrected rrnDB)
#kraken_long <- replace(kraken_long,is.na(kraken_long),1.94)
kraken_long <- mutate(kraken_long, tempabundance=rel_abundance*mean16S/100)
donor_total <- kraken_long %>%
  group_by(Sample) %>%
  summarise_at(vars(tempabundance), list(totaltempabundance=sum))
colnames(allplate)[colnames(allplate) == "SampleName"] <- "Sample"
donor_total <- merge(donor_total, allplate %>% select(Sample, Copiesgram, Condition), by="Sample")
donor_total <- mutate(donor_total, TotalBacteria=Copiesgram/totaltempabundance)
donor_total <- donor_total %>%
  group_by(Sample) %>%
  summarise_at(vars(TotalBacteria), list(TotalBacteria=mean))
donor_total <- donor_total %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)

#Ordering conditions
donor_total <- mutate(donor_total, PlotOrder=ifelse(Condition == "NF", 1, 
                                                  ifelse(Condition == "OF", 2, 
                                                         ifelse(Condition == "OR", 3, 
                                                                ifelse(Condition == "OH", 4, 
                                                                       ifelse(Condition == "ZF", 5,
                                                                              ifelse(Condition == "ZR", 6, 7)))))))
#Color palette for conditions
condition_palette <- c("#762983","#9a6faa","#c3a5d0","#acaaaf","#7ebd42","#4d9222","#26641a") 
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Export total bacteria/g of sample file
#write.table(donor_total, here("qPCR/totalbacteria_pergram.tsv"), row.names = FALSE, sep="\t", quote=FALSE)  

#Plotting total bacteria/g of sample (not per taxon)
ggplot(donor_total %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=TotalBacteria, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Total Bacteria/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+8, 1e+19), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

#ggsave(here("qPCR/2023/outputs/totalbacteriamean.jpg"), dpi=300, w = 5, h = 5)

#Calculate bacteria/taxon
mega_total <- merge(donor_total,kraken_long, by="Sample")
mega_total <- mutate(mega_total, Countgram=((rel_abundance*TotalBacteria)/100))
mega_total <- mega_total %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)

f <- ggplot(mega_total %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")) %>%filter(Phylum=="Firmicutes"| Phylum=="Bacteroidetes"), aes(x=reorder(Condition, PlotOrder), y=Countgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  #geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_boxplot_pattern(aes(pattern=Phylum)) +
  #stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     #tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Total Bacteria/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+20), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 
f
#ggsave(here("qPCR/Plots/qPCRFirmicutesBacteroidetes.jpg"), dpi=300, w = 5, h = 5)

#Calculating Firmicutes/Bacteroidetes Ratio by Condition
fb <- mega_total %>% filter(Phylum %in% c("Firmicutes","Bacteroidetes"))
fb <- fb %>% select(Sample, Phylum, Countgram)
fb <- spread(fb, key="Phylum", value="Countgram") #change df format to have F and B as columns
fbratio <- mutate(fb, FBRatio=Firmicutes/Bacteroidetes)
fbratio <- fbratio %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)

fbratio <- mutate(fbratio, PlotOrder=ifelse(Condition == "NF", 1, 
                                        ifelse(Condition == "OF", 2, 
                                               ifelse(Condition == "OR", 3, 
                                                      ifelse(Condition == "OH", 4, 
                                                             ifelse(Condition == "ZF", 5,
                                                                    ifelse(Condition == "ZR", 6, 7)))))))

fb <- ggplot(fbratio %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=FBRatio, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "# of Firmicutes/# of Bacteroidetes"
  ) +
  theme_bw() +
  ylim(0,20) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 
fb
ggsave(here("qPCR/Plots/qPCRFBRatio.jpg"), dpi=300, w=5, h=5)

fandb <- plot_grid(f,b, nrow=1, ncol=2, align="v", axis='lr')
fandb
fbratio<-plot_grid(fb, fandb, nrow=2, ncol=1)
fbratio

##### Absolute Counts at Genus Level#####
#Create a unified dataframe containing all Plate data
allcopies <- read.csv(here("qPCR/2023/allsamplecopiesul.csv"), sep=",", header=TRUE)
allplate <- allcopies %>% select(SampleName,CopyNumber, SampleNameFull, Condition, Cq)

#Calculate copies/gram for each sample
massnf <- (50/0.062)
massomni <- (50/0.013)
masszymo <- (50/0.006)

allplate <- mutate(allplate, Preservative=substr(Condition,1,1))
allplate <- mutate(allplate, Copiesgram=ifelse(Preservative =="Z",CopyNumber*masszymo,
                                               ifelse(Preservative=="O", CopyNumber*massomni,
                                                      ifelse(Preservative=="N", CopyNumber*massnf, NA)
                                               )))


#Read in kraken outputs and rrnDB copy number data
krakengenus <- read.table(here("DNA/2.kraken/kraken2_classification/processed_results_krakenonly/taxonomy_matrices_classified_only/kraken_genus_percentage.txt"), sep="\t", header=TRUE)
copyndb <- read.csv(here("qPCR/genus_copynumber_table.tsv"), sep="\t", header=TRUE)

krakengenus$Genera <- row.names(krakengenus)
row.names(krakengenus) <- NULL
kraken_long <- melt(krakengenus, id.vars = "Genera", variable.name = "Sample", value.name = "rel_abundance")
kraken_long <- mutate(kraken_long, Sample=gsub("\\.", "_", Sample))
colnames(kraken_long)[colnames(kraken_long) == "Genera"] <- "Genus"
colnames(copyndb)[colnames(copyndb) == "name"] <- "Genus"
kraken_long <- merge(kraken_long, copyndb, by="Genus", all.x=TRUE)

#Calculating total bacteria per sample (median and corrected rrnDB)
#FIXME
#kraken_long <- replace(kraken_long,is.na(kraken_long),1.94)
kraken_long <- mutate(kraken_long, tempabundance=rel_abundance*mean16S/100)
donor_total <- kraken_long %>%
  group_by(Sample) %>%
  summarise_at(vars(tempabundance), list(totaltempabundance=sum))
colnames(allplate)[colnames(allplate) == "SampleName"] <- "Sample"
donor_total <- merge(donor_total, allplate %>% select(Sample, Copiesgram, Condition), by="Sample")
donor_total <- mutate(donor_total, TotalBacteria=Copiesgram/totaltempabundance)
donor_total <- donor_total %>%
  group_by(Sample) %>%
  summarise_at(vars(TotalBacteria), list(TotalBacteria=median))
donor_total <- donor_total %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)

#Ordering conditions
donor_total <- mutate(donor_total, PlotOrder=ifelse(Condition == "NF", 1, 
                                                    ifelse(Condition == "OF", 2, 
                                                           ifelse(Condition == "OR", 3, 
                                                                  ifelse(Condition == "OH", 4, 
                                                                         ifelse(Condition == "ZF", 5,
                                                                                ifelse(Condition == "ZR", 6, 7)))))))
#Color palette for conditions
condition_palette <- c("#7c1836","#a22a5e","#c36599","#acaaaf","#bfcd6e","#9dad34","#85950f")
names(condition_palette) <- c("OH", "OR", "OF", "NF", "ZF", "ZR", "ZH")

#Plotting total bacteria/g of sample (not per taxon)
ggplot(donor_total %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=TotalBacteria, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Total Bacteria/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+19), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

#ggsave(here("qPCR/Plots/qPCRtotalmediangenus.jpg"), dpi=300, w = 5, h = 5)

#Calculate bacteria/taxon
mega_total <- merge(donor_total,kraken_long, by="Sample")
mega_total <- mutate(mega_total, Countgram=((rel_abundance*TotalBacteria)/100))
mega_total <- mega_total %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)
mega_total <- mutate(mega_total, Countgram=ifelse(Countgram==0, 1, Countgram))

ggplot(mega_total %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")) %>%filter(Genus=="Treponema"), aes(x=reorder(Condition, PlotOrder), y=Countgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Total Bacteria/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1, 1e+14), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

ggsave(here("qPCR/Plots/qPCRTreponema.jpg"), dpi=300, w = 5, h = 5)
