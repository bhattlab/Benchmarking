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

#####Plate 1: Standard Curve#####
#Create and edit unified sample datafile
qPCR_samplewellp1 <- read.csv(here("qPCR/Plate_1/Plate1_Layout.csv"), sep=",", header=TRUE)
qPCR_datap1 <- read.csv(here("qPCR/Plate_1/HotPooqPCR_Plate1.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp1)[colnames(qPCR_samplewellp1) == "Well384"] <- "Well"
qPCRplate1 <- merge(qPCR_samplewellp1,qPCR_datap1, by="Well")

#Export file
#write.table(qPCRplate1, here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/qPCRplate1.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 1
qPCRplate1 <- qPCRplate1 %>% filter(grepl("standard", SampleName))
qPCRplate1 <- filter(qPCRplate1, PCR_Replicate!="Rep4")
qPCRplate1 <- qPCRplate1 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate1 <- mutate(qPCRplate1, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate1 <- mutate(qPCRplate1, Cq=as.numeric(Cq))
qPCRplate1 <- mutate(qPCRplate1, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate1 <- mutate(qPCRplate1, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Filter for B vulgatus and calculate log10
qPCRplate1 <- qPCRplate1 %>% filter(Standard=="B vulgatus")
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate1, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point() +
  theme_bw()

#ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/Plots/qPCRplate1_stdcurve_bv1.jpg"), dpi=300, w = 5, h = 5)

#Filter for F prausnitizii and calculate log10
qPCRplate1 <- qPCRplate1 %>% filter(Standard=="F prausnitzii")
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number


library(ggpmisc)
ggplot(data = qPCRplate1, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point() +
  theme_bw()

#ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/Plots/qPCRplate1_stdcurve_fp1.jpg"), dpi=300, w = 5, h = 5)
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

#Standard curve using both B vulgatus and F prausnitizii data
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate1, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=Standard)) +
  theme_bw()

#ggsave(here("qPCR/Plate_1/Plots/qPCRplate1_stdcurve.jpg"), dpi=300, w = 6, h = 5)
#To introduce custom colors in the plot

#To introduce custom colors in the plot
# library(ggpmisc)
# standard_palette <- c("#BCA88E", "#9DA993")
# names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#####Plate 1: Absolute Abundance Calculation#####
#Edit dataframe
qPCRplate1 <- read.csv(here("qPCR/Plate_1/qPCR_plate1.csv"), sep=",", header=TRUE)
qPCRplate1 <- qPCRplate1 %>% filter(!grepl("standard", SampleName))
qPCRplate1 <- filter(qPCRplate1, PCR_Replicate!="Rep4")
qPCRplate1 <- mutate(qPCRplate1, Cq=as.numeric(Cq))

#Calculation using combined standard curve model
#Standard curve is: y=35-2.93x
qPCRplate1 <- qPCRplate1 %>% mutate(logCopyNumber=(Cq-35)/-2.93)
qPCRplate1 <- qPCRplate1 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate1 <- qPCRplate1 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')

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

#Filter out irregular signal samples
irregularqPCR <- qPCR_weird$SampleNameFull
qPCRplate1 <- qPCRplate1 %>% filter(!(SampleNameFull %in% irregularqPCR))

#Export datafile removing irregular samples
#write.table(qPCRplate1, here("qPCR/Plate_1/qPCR_plate1.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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
  scale_y_log10(limits=c(1e+5, 1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none")

#ggsave(here("qPCR/Plate_1/Plots/qPCRplate1_copiespercondition.jpg"), dpi=300, w = 5, h = 5)

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
ggsave(here("qPCR/Plate_1/Plots/qPCRplate1_copiesgramcondition.jpg"), dpi=300, w = 5, h = 5)

#Export final plate 1 datafile
#write.table(qPCRplate1, here("qPCR/Plate_1/qPCR_plate1_final.csv"), row.names = FALSE, sep=",", quote=FALSE)  


######Plate 2: Standard Curve and Absolute Count#####
#Create and edit unified sample datafile
qPCR_samplewellp2 <- read.csv(here("qPCR/Plate_2/qPCR_plate2_layout.csv"), sep=",", header=TRUE)
qPCR_datap2 <- read.csv(here("qPCR/Plate_2/HotPooqPCR_Plate2.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp2)[colnames(qPCR_samplewellp2) == "Well384"] <- "Well"
qPCRplate2 <- merge(qPCR_samplewellp2,qPCR_datap2, by="Well")

#Export file
#write.table(qPCRplate2, here("Documents/Bhatt/Benchmarking/qPCR/Plate_2/qPCRplate2.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 2
qPCRplate2 <- qPCRplate2 %>% filter(grepl("standard", SampleName))
qPCRplate2 <- filter(qPCRplate2, PCR_Replicate!="Rep4")
qPCRplate2 <- qPCRplate2 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate2 <- mutate(qPCRplate2, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate2 <- mutate(qPCRplate2, Cq=as.numeric(Cq))
qPCRplate2 <- mutate(qPCRplate2, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate2 <- mutate(qPCRplate2, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Filter for B vulgatus and calculate log10
qPCRplate2 <- qPCRplate2 %>% filter(Standard=="B vulgatus")
qPCRplate2 <- qPCRplate2 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate2, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point() +
  theme_bw()

#ggsave(here("qPCR/Plate_2/Plots/qPCRplate2_stdcurve_bv2.jpg"), dpi=300, w = 5, h = 5)

#Filter for F prausnitizii and calculate log10
qPCRplate2 <- qPCRplate2 %>% filter(Standard=="F prausnitzii")
qPCRplate2 <- qPCRplate2 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number


library(ggpmisc)
ggplot(data = qPCRplate2, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point() +
  theme_bw()

#ggsave(here("qPCR/Plate_2/Plots/qPCRplate2_stdcurve_fp2.jpg"), dpi=300, w = 5, h = 5)

#Standard curve using both B vulgatus and F prausnitizii data
qPCRplate2 <- qPCRplate2 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate2, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=Standard)) +
  theme_bw()

ggsave(here("qPCR/Plate_2/Plots/qPCRplate2_stdcurve_all.jpg"), dpi=300, w = 6, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate2 <- read.csv(here("qPCR/Plate_2/qPCRplate2.csv"), sep=",", header = TRUE)
qPCRplate2 <- qPCRplate2 %>% filter(!grepl("standard", SampleName))
qPCRplate2 <- filter(qPCRplate2, PCR_Replicate!="Rep4")
qPCRplate2 <- mutate(qPCRplate2, Cq=as.numeric(Cq))

#Calculation using combined standard curve model
#Standard curve is: y=35.1-3.01x
qPCRplate2 <- qPCRplate2 %>% mutate(logCopyNumber=(Cq-35.1)/-3.01)
qPCRplate2 <- qPCRplate2 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate2 <- qPCRplate2 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')

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

#Filter out irregular signal samples
irregularqPCR <- qPCR_weird$SampleNameFull
qPCRplate2 <- qPCRplate2 %>% filter(!(SampleNameFull %in% irregularqPCR))

#Export datafile removing irregular samples
#write.table(qPCRplate2, here("qPCR/Plate_2/qPCR_plate2.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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

ggsave(here("qPCR/Plate_2/Plots/qPCRplate2_copiespercondition.jpg"), dpi=300, w = 5, h = 5)

##### Plate 3: Standard Curve and Absolute Count #####
#Create and edit unified sample datafile
qPCR_samplewellp3 <- read.csv(here("qPCR/Plate_3/qPCR_plate3_layout.csv"), sep=",", header=TRUE)
qPCR_datap3 <- read.csv(here("qPCR/Plate_3/HotPooqPCR_Plate3.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp3)[colnames(qPCR_samplewellp3) == "Well384"] <- "Well"
qPCRplate3 <- merge(qPCR_samplewellp3,qPCR_datap3, by="Well")

#Export file
#write.table(qPCRplate3, here("qPCR/Plate_3/qPCRplate3.csv"), row.names = FALSE, sep=",", quote=FALSE)  

#Make standard curve for plate 3
qPCRplate3 <- qPCRplate3 %>% filter(grepl("standard", SampleName))
qPCRplate3 <- filter(qPCRplate3, PCR_Replicate!="Rep4")
qPCRplate3 <- qPCRplate3 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate3 <- mutate(qPCRplate3, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate3 <- mutate(qPCRplate3, Cq=as.numeric(Cq))
qPCRplate3 <- mutate(qPCRplate3, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate3 <- mutate(qPCRplate3, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

#Filter for B vulgatus and calculate log10
qPCRplate3 <- qPCRplate3 %>% filter(Standard=="B vulgatus")
qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point() +
  theme_bw()

ggsave(here("qPCR/Plate_3/Plots/qPCRplate3_stdcurve_bv.jpg"), dpi=300, w = 5, h = 5)

#Filter for F prausnitizii and calculate log10
qPCRplate3 <- qPCRplate3 %>% filter(Standard=="F prausnitzii")
qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number


library(ggpmisc)
ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point() +
  theme_bw()

ggsave(here("qPCR/Plate_3/Plots/qPCRplate3_stdcurve_fp3.jpg"), dpi=300, w = 5, h = 5)

#Standard curve using both B vulgatus and F prausnitizii data
qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=Standard)) +
  theme_bw()

ggsave(here("qPCR/Plate_3/Plots/qPCRplate3_stdcurve_all.jpg"), dpi=300, w = 6, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate3 <- read.csv(here("qPCR/Plate_3/qPCRplate3.csv"), sep=",", header = TRUE)
qPCRplate3 <- qPCRplate3 %>% filter(!grepl("standard", SampleName))
qPCRplate3 <- filter(qPCRplate3, PCR_Replicate!="Rep4")
qPCRplate3 <- mutate(qPCRplate3, Cq=as.numeric(Cq))

#Calculation using combined standard curve model
#Standard curve is: y=35.5-3.13x
qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=(Cq-35.5)/-3.13)
qPCRplate3 <- qPCRplate3 %>% mutate(CopyNumber=((10^(logCopyNumber))*1000)/6) #copies/uL
qPCRplate3 <- qPCRplate3 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE, sep='_')

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

#Filter out irregular signal samples
irregularqPCR <- qPCR_weird$SampleNameFull
qPCRplate3 <- qPCRplate3 %>% filter(!(SampleNameFull %in% irregularqPCR))

#Export datafile removing irregular samples
write.table(qPCRplate3, here("qPCR/Plate_3/qPCR_plate3.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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

ggsave(here("qPCR/Plate_3/Plots/qPCRplate3_copiespercondition.jpg"), dpi=300, w = 5, h = 5)

##### Plate 4: Standard Curve and Absolute Count#####
#Create and edit unified sample datafile
qPCR_samplewellp4 <- read.csv(here("qPCR/Plate_4/Plate4_Layout.tsv"), sep="\t", header=TRUE)
qPCR_datap4 <- read.csv(here("qPCR/Plate_4/HotPooqPCR_Plate4.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp4)[colnames(qPCR_samplewellp4) == "Well384"] <- "Well"
qPCRplate4 <- merge(qPCR_samplewellp4,qPCR_datap4, by="Well")

#Export file
#write.table(qPCRplate4, here("qPCR/Plate_4/qPCRplate4.tsv"), row.names = FALSE, sep="\t", quote=FALSE)  

#Make standard curve for plate 4
qPCRplate4 <- qPCRplate4 %>% filter(grepl("standard", SampleName))
qPCRplate4 <- filter(qPCRplate4, PCR_Replicate!="Rep4")
qPCRplate4 <- qPCRplate4 %>% separate(SampleName, c("Standard", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate4 <- mutate(qPCRplate4, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate4 <- mutate(qPCRplate4, Cq=as.numeric(Cq))
qPCRplate4 <- mutate(qPCRplate4, Standard=gsub("standardB","B vulgatus",gsub("standardF","F prausnitzii", Standard)))

#Include copy number for the standards
copyf <- 2.456*10^11
copyb <- 2.171*10^11
qPCRplate4 <- mutate(qPCRplate4, CopyNumber=ifelse(Standard=="B vulgatus", copyb/DilutionFactor, copyf/DilutionFactor))

# #Filter for B vulgatus and calculate log10
# qPCRplate3 <- qPCRplate3 %>% filter(Standard=="B vulgatus")
# qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
# 
# library(ggpmisc)
# ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
#   geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
#   stat_poly_eq(formula  = y~x,
#                eq.with.lhs = "italic(y)~`=`~",
#                eq.x.rhs    = "~italic(x)",
#                aes(label   = paste(..eq.label..)), 
#                parse = TRUE) +         
#   geom_point() +
#   theme_bw()
# 
# ggsave(here("qPCR/Plate_3/Plots/qPCRplate3_stdcurve_bv.jpg"), dpi=300, w = 5, h = 5)
# 
# #Filter for F prausnitizii and calculate log10
# qPCRplate3 <- qPCRplate3 %>% filter(Standard=="F prausnitzii")
# qPCRplate3 <- qPCRplate3 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number
# 
# 
# library(ggpmisc)
# ggplot(data = qPCRplate3, aes(x = logCopyNumber, y = Cq)) +
#   geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
#   stat_poly_eq(formula  = y~x,
#                eq.with.lhs = "italic(y)~`=`~",
#                eq.x.rhs    = "~italic(x)",
#                aes(label   = paste(..eq.label..)), 
#                parse = TRUE) +         
#   geom_point() +
#   theme_bw()
# 
# ggsave(here("qPCR/Plate_3/Plots/qPCRplate3_stdcurve_fp3.jpg"), dpi=300, w = 5, h = 5)

#Standard curve using both B vulgatus and F prausnitizii data
qPCRplate4 <- qPCRplate4 %>% mutate(logCopyNumber=log10(CopyNumber)) #log10 of the copy number

library(ggpmisc)
ggplot(data = qPCRplate4, aes(x = logCopyNumber, y = Cq)) +
  geom_smooth(method = "lm", se=FALSE, formula = y~x) + 
  stat_poly_eq(formula  = y~x,
               eq.with.lhs = "italic(y)~`=`~",
               eq.x.rhs    = "~italic(x)",
               aes(label   = paste(..eq.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=Standard)) +
  theme_bw()

ggsave(here("qPCR/Plate_4/Plots/qPCRplate4_stdcurve_all.jpg"), dpi=300, w = 6, h = 5)

# #Check dilutions by logcopynumber
# #Edit dataframe for absolute abundance calculations
# qPCRplate4 <- read.csv(here("qPCR/Plate_4/qPCRplate4.tsv"), sep="\t", header = TRUE)
# qPCRplate4 <- qPCRplate4 %>% filter(!grepl("standard", SampleName))
# qPCRplate4 <- filter(qPCRplate4, PCR_Replicate!="Rep4")
# qPCRplate4 <- mutate(qPCRplate4, Cq=as.numeric(Cq))
# 
# qPCRplate4 <- qPCRplate4 %>% mutate(logCopyNumber=(Cq-31.9)/-3.02)
# qPCRplate4 <- qPCRplate4 %>% separate(SampleName, c("Donor", "Condition", "Replicate", "DilutionFactor"), remove=FALSE, sep='_')
# #qPCRplate4 <- mutate(qPCRplate4, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
# #qPCRplate4 <- qPCRplate4 %>% mutate(CopyNumber=((10^(logCopyNumber))*DilutionFactor)/6) #copies/uL
# qPCRplate4 <- qPCRplate4 %>% filter(!grepl("water", SampleName))
# qPCRplate4 <- qPCRplate4 %>% filter(!grepl("buffer", SampleName))
# qPCRplate4 <- qPCRplate4 %>% unite(Sample, Donor, Condition, Replicate, sep ="_")
# qPCRplate4 <- qPCRplate4[!is.na(qPCRplate4$logCopyNumber),]
# #qPCRplate4 <- mutate(qPCRplate4, DilutionFactor=as.factor(DilutionFactor))
# 
# #Plotting copies/uL by Sample
# ggplot(qPCRplate4, aes(x=Sample, y=logCopyNumber, color=DilutionFactor)) +
#   geom_point() + 
#   labs(
#     x = "Sample",
#     y = "logCopyNumber"
#   ) +
#   theme_bw() +
#   scale_color_discrete(breaks=c("10","100","1,000", "10,000", "100,000", "1,000,000")) #reorder legend
#   #scale_y(limits=c(1, 10), n.breaks=1, minor_breaks=NULL)
# 
# 
# ggsave(here("qPCR/Plate_4/Plots/qPCRplate4_checkdilution.jpg"), dpi=300, w = 7, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate4 <- read.csv(here("qPCR/Plate_4/qPCRplate4.tsv"), sep="\t", header = TRUE)
qPCRplate4 <- qPCRplate4 %>% filter(!grepl("standard", SampleName))
qPCRplate4 <- filter(qPCRplate4, PCR_Replicate!="Rep4")
qPCRplate4 <- mutate(qPCRplate4, Cq=as.numeric(Cq))
#Calculation using combined standard curve model
#Standard curve is: y=31.9-3.02x
qPCRplate4 <- qPCRplate4 %>% mutate(logCopyNumber=(Cq-31.9)/-3.02)
qPCRplate4 <- qPCRplate4 %>% separate(SampleName, c("Donor", "Condition", "Replicate", "DilutionFactor"), remove=FALSE, sep='_')
qPCRplate4 <- mutate(qPCRplate4, DilutionFactor=as.numeric(gsub(",","",DilutionFactor)))
qPCRplate4 <- qPCRplate4 %>% mutate(CopyNumber=((10^(logCopyNumber))*DilutionFactor)/6) #copies/uL
qPCRplate4 <- qPCRplate4 %>% filter(!grepl("water", SampleName))
qPCRplate4 <- qPCRplate4 %>% filter(!grepl("buffer", SampleName))
qPCRplate4 <- qPCRplate4 %>% unite(Sample, Donor, Condition, Replicate, sep ="_")
qPCRplate4 <- qPCRplate4[!is.na(qPCRplate4$CopyNumber),]
qPCRplate4 <- mutate(qPCRplate4, DilutionFactor=as.factor(DilutionFactor))

#Plotting copies/uL by Sample
ggplot(qPCRplate4, aes(x=Sample, y=CopyNumber, color=DilutionFactor)) +
  geom_point() + 
  labs(
    x = "Sample",
    y = "Copies/uL"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+11), n.breaks=1e1, minor_breaks=NULL)

ggsave(here("qPCR/Plate_4/Plots/qPCRplate4_copiespersample.jpg"), dpi=300, w = 7, h = 5)

#Remove unnecessary dilutions - keep 1:1000 or 1:100 for D05_ZH_R3 and D09_ZR_R1
plate4 <- qPCRplate4 %>% filter(SampleName %in% c("D01_ZH_R1_1,000", "D03_ZF_R3_1,000", "D05_ZH_R3_100", "D07_NF_R3_1,000", "D09_ZH_R1_1,000", "D09_ZR_R1_100"))

#Cleanup plate4 dataframe: change column names and create variables needed for merging all data together
colnames(plate4)[colnames(plate4) == "SampleNameFull"] <- "DiluteSampleName"
colnames(plate4)[colnames(plate4) == "SampleName"] <- "OGSampleName"
colnames(plate4)[colnames(plate4) == "Sample"] <- "SampleName"
plate4$SampleNameFull <- paste(plate4$SampleName, plate4$PCR_Replicate, sep="_") #create samplename without dilution factor in name
plate4 <- plate4 %>% separate(SampleName, c("Donor", "Condition", "Replicate"), remove=FALSE)

#Export datafile with copy numbers
write.table(plate4, here("qPCR/Plate_4/qPCR_plate4.tsv"), row.names = FALSE, sep="\t", quote=FALSE)  

#Plate 4: Scatterplots for consistency
qPCR4 <- spread(qPCRplate4 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)

ggplot(qPCR4, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Replicate 1 Cq",
    y = "Replicate 2 Cq"
  ) 

ggplot(qPCR4, aes(x=Rep2, y=Rep3)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Replicate 2 Cq",
    y = "Replicate 3 Cq"
  ) 

ggplot(qPCR4, aes(x=Rep1, y=Rep3)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Replicate 1 Cq",
    y = "Replicate 3 Cq"
  ) 

#####Scatterplots for qPCR Consistency#####
qPCR2 <- spread(qPCRplate2 %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)

ggplot(qPCR2, aes(x=Rep1, y=Rep2)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Replicate 1 Cq",
    y = "Replicate 2 Cq"
  ) 

ggsave(here("qPCR/Plots/RepCompare_12.jpg"), dpi=300, w = 5, h = 5)

ggplot(qPCR2, aes(x=Rep1, y=Rep3)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Replicate 1 Cq",
    y = "Replicate 3 Cq"
  ) 
ggsave(here("qPCR/Plots/RepCompare_13.jpg"), dpi=300, w = 5, h = 5)

ggplot(qPCR2, aes(x=Rep2, y=Rep3)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept=0) +
  labs(
    x = "Replicate 2 Cq",
    y = "Replicate 3 Cq"
  ) 
ggsave(here("qPCR/Plots/RepCompare_23.jpg"), dpi=300, w = 5, h = 5)

#Irregular curves: samples with 1/3 failed replicates
qPCR_w <- spread(qPCR_weird %>% select(SampleName, Cq, PCR_Replicate), key=PCR_Replicate, value=Cq)

##### Absolute Counts at Phylum Level#####
#Create a unified dataframe containing all Plate data
plate1 <- read.csv(here("qPCR/Plate_1/qPCR_plate1.csv"), sep=",", header=TRUE)
plate2 <- read.csv(here("qPCR/Plate_2/qPCR_plate2.csv"), sep=",", header=TRUE)
plate3 <- read.csv(here("qPCR/Plate_3/qPCR_plate3.csv"), sep=",", header=TRUE)
plate4 <- read.csv(here("qPCR/Plate_4/qPCR_plate4.tsv"), sep="\t", header=TRUE)

allplate <- rbind(plate1 %>% select(SampleName,CopyNumber, SampleNameFull, Condition), plate2 %>% select(SampleName,CopyNumber, SampleNameFull, Condition), plate3 %>% select(SampleName,CopyNumber, SampleNameFull, Condition), plate4 %>% select(SampleName,CopyNumber, SampleNameFull, Condition))
#Calculate copies/gram for each sample
massnf <- (50/0.062)
massomni <- (50/0.037)
masszymo <- (50/0.08)

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
  scale_y_log10(limits=c(1e+5, 1e+19), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

ggsave(here("qPCR/Plots/qPCRtotalbacteriamedian.jpg"), dpi=300, w = 5, h = 5)

#Calculate bacteria/taxon
mega_total <- merge(donor_total,kraken_long, by="Sample")
mega_total <- mutate(mega_total, Countgram=((rel_abundance*TotalBacteria)/100))
mega_total <- mega_total %>% separate(Sample, c("Donor", "Condition", "Replicate"), remove=FALSE)

ggplot(mega_total %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")) %>%filter(Phylum=="Bacteroidetes"), aes(x=reorder(Condition, PlotOrder), y=Countgram, color=Condition)) +
  geom_jitter(width=0.1) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(comparisons=list(c("OF","OR"), c("OR", "OH"), c("OF","OH"), c("ZF", "ZR"), c("ZR", "ZH"), c("ZF", "ZH"), c("ZF", "NF"), c("OF", "NF")),
                     tip.length = 0, label="p.signif") +
  labs(
    x = "Condition",
    y = "Total Bacteria/g"
  ) +
  theme_bw() +
  scale_y_log10(limits=c(1e+5, 1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

ggsave(here("qPCR/Plots/qPCRBacteroidetes.jpg"), dpi=300, w = 5, h = 5)

#Calculating Firmicute/Bacteroidetes Ratio by Condition
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

ggplot(fbratio %>% filter(Condition %in% c("NF", "OF", "ZF", "OR", "ZR", "OH", "ZH")), aes(x=reorder(Condition, PlotOrder), y=FBRatio, color=Condition)) +
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

ggsave(here("qPCR/Plots/qPCRFBRatio.jpg"), dpi=300, w=5, h=5)

##### Absolute Counts at Genus Level#####
#Create a unified dataframe containing all Plate data
plate1 <- read.csv(here("qPCR/Plate_1/qPCR_plate1.csv"), sep=",", header=TRUE)
plate2 <- read.csv(here("qPCR/Plate_2/qPCR_plate2.csv"), sep=",", header=TRUE)
plate3 <- read.csv(here("qPCR/Plate_3/qPCR_plate3.csv"), sep=",", header=TRUE)
plate4 <- read.csv(here("qPCR/Plate_4/qPCR_plate4.tsv"), sep="\t", header=TRUE)

allplate <- rbind(plate1 %>% select(SampleName,CopyNumber, SampleNameFull, Condition), plate2 %>% select(SampleName,CopyNumber, SampleNameFull, Condition), plate3 %>% select(SampleName,CopyNumber, SampleNameFull, Condition), plate4 %>% select(SampleName,CopyNumber, SampleNameFull, Condition))
#Calculate copies/gram for each sample
massnf <- (50/0.062)
massomni <- (50/0.037)
masszymo <- (50/0.08)

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
  scale_y_log10(limits=c(1e+5, 1e+16), n.breaks=1e1, minor_breaks=NULL) +
  scale_color_manual(values=condition_palette) +
  theme(legend.position = "none") 

ggsave(here("qPCR/Plots/qPCRtotalmediangenus.jpg"), dpi=300, w = 5, h = 5)

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
#####Irregular curve analysis#####
#Merge dataframes with irregular curves, plate layout, and DNA concentration
qPCR_weird <- read.csv(here("qPCR/Irregular qPCR Signal.csv"), sep=",", header=TRUE)
qPCR_layout <-read.csv(here("qPCR/qPCR384layout.csv"), sep=",", header=TRUE)
qPCR_weird <- merge(qPCR_weird, qPCR_layout, c("Well384","Plate"))
DNA_conc <- read.csv(here("qPCR/BenchmarkingDNAConcentration.csv"), sep=",", header=TRUE)
qPCR_weird <- merge(qPCR_weird, DNA_conc, by="SampleName")

#Count the number of times a sample fails
qPCR_weirdd <- qPCR_weird %>% count(SampleName)
qPCR_weird <- merge(qPCR_weird, qPCR_weirdd, by="SampleName")
qPCR_weirddd <- qPCR_weirdd %>% count(n, name="FailCount") #across all plates, how many times do 1/3, 2/3, or 3/3 replicates fail

#Formatting
qPCR_weird <- mutate(qPCR_weird, n=as.character(n))

#Plot the DNA concentration by sample replicate
ggplot(qPCR_weird, aes(x=Plate, y=DNAConcentration, color=n)) +
  #geom_point() + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  labs(
    x = "Plate",
    y = "DNA Concentration (ng/uL)",
    color = "# of Failed Replicates"
  ) +
  theme_bw()

ggsave(here("Plots/qPCRweird_byDNAconc2.jpg"), dpi=300, w = 5, h = 5)


#####Outlier Analaysis#####
plate3 <- plate3 %>% filter(SampleName %in% c("D10_NF_R3","D03_OF_R3"))
plate2 <- plate2 %>% filter(SampleName %in% c("D02_OH_R2","D03_OR_R2","D09_OF_R2"))

outliers <- rbind(plate3 %>% select(SampleName, SampleNameFull, Cq, CopyNumber, PCR_Replicate), plate2 %>% select(SampleName, SampleNameFull, Cq, CopyNumber, PCR_Replicate))

ggplot(outliers, aes(x=SampleName, y=Cq, color=PCR_Replicate)) +
  geom_point(width=0.1) + 
  labs(
    x = "Sample",
    y = "Cq"
  ) +
  theme_bw() 

ggsave(here("qPCR/Plots/qPCRoutliers.jpg"), dpi=300, w = 7, h = 5)

#Merge with DNA concentration table
outliers <- merge(DNA_conc, outliers, by="SampleName")
ggplot(outliers, aes(x=SampleName, y=DNAConcentration)) +
  geom_point(width=0.1) + 
  labs(
    x = "Sample",
    y = "DNA Concentration"
  ) +
  theme_bw() 

ggsave(here("qPCR/Plots/qPCRoutliersconc.jpg"), dpi=300, w = 5, h = 5)

#####DNA Concentration v Absolute Count#####
#Absolute count data
abscount <- read.csv(here("qPCR/totalbacteria_pergram.tsv"), sep="\t", header=TRUE)
colnames(abscount)[colnames(abscount) == "Sample"] <- "SampleName"


#DNA Concentration data
dnacon <- read.csv(here("qPCR/BenchmarkingDNAConcentration.csv"), sep=",", header=TRUE)

absdna <- merge(abscount, dnacon, by="SampleName")

ggplot(absdna, aes(x=DNAConcentration, y=TotalBacteria)) +
  geom_point() +
  theme_bw() +
  #geom_abline(slope=1, intercept=0) +
  labs(
    x = "DNA Concentration ng/uL",
    y = "Total Bacteria/g Dry Stool"
  ) +
  geom_smooth(method=lm, se=FALSE) +
  scale_y_log10()
  

