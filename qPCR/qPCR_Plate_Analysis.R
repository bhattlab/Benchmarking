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
qPCR_samplewellp1 <- read.csv(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/Plate1_Layout.csv"), sep=",", header=TRUE)
qPCR_datap1 <- read.csv(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/HotPooqPCR_Plate1.csv"), sep=",", header=TRUE)

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

ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/Plots/qPCRplate1_stdcurve_bv1.jpg"), dpi=300, w = 5, h = 5)

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

ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/Plots/qPCRplate1_stdcurve_fp1.jpg"), dpi=300, w = 5, h = 5)
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
  geom_point() +
  theme_bw()

ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_1/Plots/qPCRplate1_stdcurve.jpg"), dpi=300, w = 5, h = 5)
#To introduce custom colors in the plot

#To introduce custom colors in the plot
# library(ggpmisc)
# standard_palette <- c("#BCA88E", "#9DA993")
# names(standard_palette) <- c("B vulgatus", "F prausnitzii")

#####Plate 1: Absolute Abundance Calculation#####
#Edit dataframe
qPCRplate1 <- read.csv(here("Plate_1/qPCRplate1.csv"), sep=",", header=TRUE)
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
#write.table(qPCRplate1, here("Plate_1/qPCR_plate1.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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

ggsave(here("Plate_1/Plots/qPCRplate1_copiespercondition.jpg"), dpi=300, w = 5, h = 5)

######Plate 2: Standard Curve and Absolute Count#####
#Create and edit unified sample datafile
qPCR_samplewellp2 <- read.csv(here("Documents/Bhatt/Benchmarking/qPCR/Plate_2/qPCR_plate2_layout.csv"), sep=",", header=TRUE)
qPCR_datap2 <- read.csv(here("Documents/Bhatt/Benchmarking/qPCR/Plate_2/HotPooqPCR_Plate2.csv"), sep=",", header=TRUE)

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

ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_2/Plots/qPCRplate2_stdcurve_bv2.jpg"), dpi=300, w = 5, h = 5)

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

ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_2/Plots/qPCRplate2_stdcurve_fp2.jpg"), dpi=300, w = 5, h = 5)

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
  geom_point() +
  theme_bw()
ggsave(here("Documents/Bhatt/Benchmarking/qPCR/Plate_2/Plots/qPCRplate2_stdcurve_all.jpg"), dpi=300, w = 5, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate2 <- read.csv(here("Plate_2/qPCRplate2.csv"), sep=",", header = TRUE)
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
#write.table(qPCRplate2, here("Plate_2/qPCR_plate2.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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

ggsave(here("Plate_2/Plots/qPCRplate2_copiespercondition.jpg"), dpi=300, w = 5, h = 5)

##### Plate 3: Standard Curve and Absolute Count #####
#Create and edit unified sample datafile
qPCR_samplewellp3 <- read.csv(here("Plate_3/qPCR_plate3_layout.csv"), sep=",", header=TRUE)
qPCR_datap3 <- read.csv(here("Plate_3/HotPooqPCR_Plate3.csv"), sep=",", header=TRUE)

colnames(qPCR_samplewellp3)[colnames(qPCR_samplewellp3) == "Well384"] <- "Well"
qPCRplate3 <- merge(qPCR_samplewellp3,qPCR_datap3, by="Well")

#Export file
#write.table(qPCRplate3, here("Plate_3/qPCRplate3.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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

ggsave(here("Plate_3/Plots/qPCRplate3_stdcurve_bv.jpg"), dpi=300, w = 5, h = 5)

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

ggsave(here("Plate_3/Plots/qPCRplate3_stdcurve_fp3.jpg"), dpi=300, w = 5, h = 5)

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
  geom_point() +
  theme_bw()
ggsave(here("Plate_3/Plots/qPCRplate3_stdcurve_all.jpg"), dpi=300, w = 5, h = 5)

#Edit dataframe for absolute abundance calculations
qPCRplate3 <- read.csv(here("Plate_3/qPCRplate3.csv"), sep=",", header = TRUE)
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
#write.table(qPCRplate3, here("Plate_3/qPCR_plate3.csv"), row.names = FALSE, sep=",", quote=FALSE)  

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

ggsave(here("Plate_3/Plots/qPCRplate3_copiespercondition.jpg"), dpi=300, w = 5, h = 5)

#####Irregular curve analysis#####
#Merge dataframes with irregular curves, plate layout, and DNA concentration
qPCR_weird <- read.csv(here("Irregular qPCR Signal.csv"), sep=",", header=TRUE)
qPCR_layout <-read.csv(here("qPCR384layout.csv"), sep=",", header=TRUE)
qPCR_weird <- merge(qPCR_weird, qPCR_layout, c("Well384","Plate"))
DNA_conc <- read.csv(here("BenchmarkingDNAConcentration.csv"), sep=",", header=TRUE)
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

