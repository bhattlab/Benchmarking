library(tidyverse)
library(here)
# # Pilot Conversion
# layout96 <- read.csv("/Users/dylanmaghini/scg4/projects/benchmarking/pilot_96.txt", sep="\t", header = FALSE)
# layout_conversion <- read.csv("/Users/dylanmaghini/scg4/projects/benchmarking/format_conversion.tsv", sep="\t", header = FALSE)
# names(layout96) <- c("Well96", "SampleName")
# names(layout_conversion) <- c("Well96", "Replicate", "Well384")
# 
# full_layout <- layout_conversion %>% full_join(layout96)
# full_layout <- full_layout %>% mutate(SampleNameFull = paste(SampleName, Replicate, sep="_"))
# 
# write.table(full_layout, "/Users/dylanmaghini/scg4/projects/benchmarking/pilot_384.tsv", sep="\t", row.names=FALSE)

# Full Conversion
layout96 <- read.csv(here("qPCR/PlateLayout/qPCR_plate4_96welllayout.csv"), sep=",", header = FALSE)
layout_conversion <- read.csv(here("qPCR/PlateLayout/format_conversion.tsv"), sep="\t", header = FALSE)

names(layout96) <- c("Plate", "Well96", "SampleName")
names(layout_conversion) <- c("Well96", "Replicate", "Well384")

full_layout <- layout_conversion %>% full_join(layout96)
full_layout <- na.omit(full_layout)
full_layout <- full_layout %>% mutate(SampleNameFull = paste(SampleName, Replicate, sep="_"))
names(full_layout) <- c("Well96", "PCR_Replicate", "Well384", "Plate", "SampleName", "SampleNameFull")

write.table(full_layout, here("qPCR/Plate_4/Plate4_Layout.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
