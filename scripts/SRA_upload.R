library(here)
library(dplyr)
library(tidyverse)

df <- read.table(here("outputs/tables/DNA_metadata.tsv"), header=TRUE, sep="\t")
df <- df %>% filter(Replicate == "R1")
df <- df %>% mutate(SampleName = paste(Donor, Condition, sep="_"))
df <- df %>% mutate(SampleTitle = gsub("F", "-80C", 
                                       gsub("R", "23C", 
                                            gsub("H", "40C", 
                                                 gsub("N", "No Preservative ", 
                                                      gsub("O", "omnigene ", 
                                                           gsub("Z", "zymo ", 
                                                                gsub("D", "Donor", gsub("_", " ", SampleName)))))))))
df <- df %>% mutate(Organism = "human feces metagenome")
df <- df %>% mutate(Host = "Homo sapiens")
df <- df %>% mutate(IsolationSource = "not applicable")
df <- df %>% mutate(CollectionDate = "2021")
df <- df %>% mutate(GeoLocName = "USA:California")
df <- df %>% mutate(LatLong = "37.43 N 122.17 W")
df <- df %>% mutate(SampCollectDevice = ifelse(Kit == "Omnigene", "DNA Genotek OMR-200",
                                               ifelse(Kit == "Zymo", "Zymo DNA/RNA Shield", "NA")))
df <- df %>% mutate(SampMatProcess = ifelse(Temperature == "Hot", "40C for 7 days", 
                                            ifelse(Temperature == "Room Temp", "23C for 7 days", "Immediate storage at -80C")))

df <- df %>% select(SampleName, SampleTitle, Organism, Host, IsolationSource, CollectionDate, GeoLocName, LatLong, SampCollectDevice, SampMatProcess, Donor)
write.table(df, here("outputs/tables/SRA_Biosample.tsv"), sep="\t", quote=FALSE, row.names=FALSE)




########
df <- read.table(here("outputs/tables/DNA_metadata.tsv"), header=TRUE, sep="\t")
df <- df %>% select(SampleID)
df <- df %>% mutate(SampleName=gsub("_R.", "", SampleID))
df <- df %>% mutate(LibraryID=paste(SampleID, "_DNA", sep=""))
df <- df %>% mutate(LibraryLayout="paired")
df <- df %>% mutate(Title=paste("Paired end metagenomic sequencing of ", SampleID, sep=""))

df <- df %>% mutate(LibraryStrategy="WGS")
df <- df %>% mutate(LibrarySource="METAGENOMIC")
df <- df %>% mutate(LibrarySelection="RANDOM")
df <- df %>% mutate(LibrarySource="METAGENOMIC")
df <- df %>% mutate(Platform="ILLUMINA")
df <- df %>% mutate(InstrumentModel="Illumina NovaSeq 6000")
df <- df %>% mutate(DesignDescription="2 x 150 bp metagenomic sequencing on a NovaSeq 6000 following library preparation with Illumina DNA Prep")
df <- df %>% mutate(filetype="fastq")
df <- df %>% mutate(filename=paste(gsub("_", "-", SampleID), "_1.fq.gz", sep=""))
df <- df %>% mutate(filename2=paste(gsub("_", "-", SampleID), "_2.fq.gz", sep=""))
df <- df %>% mutate(filename3=paste(gsub("_", "-", SampleID), "_orphans.fq.gz", sep=""))


# metatranscriptomic
df_t <- read.table(here("outputs/tables/DNA_metadata.tsv"), header=TRUE, sep="\t")
df_t <- df_t %>% select(SampleID)
df_t <- df_t %>% mutate(SampleName=gsub("_R.", "", SampleID))
df_t <- df_t %>% mutate(LibraryID=paste(SampleID, "_RNA", sep=""))
df_t <- df_t %>% mutate(LibraryLayout="paired")
df_t <- df_t %>% mutate(Title=paste("Paired end metatranscriptomic sequencing of ", SampleID, sep=""))

df_t <- df_t %>% mutate(LibraryStrategy="RNA-Seq")
df_t <- df_t %>% mutate(LibrarySource="METATRANSCRIPTOMIC")
df_t <- df_t %>% mutate(LibrarySelection="RANDOM")
df_t <- df_t %>% mutate(LibrarySource="METAGENOMIC")
df_t <- df_t %>% mutate(Platform="ILLUMINA")
df_t <- df_t %>% mutate(InstrumentModel="Illumina NovaSeq 6000")
df_t <- df_t %>% mutate(DesignDescription="2 x 150 bp metatranscriptomic sequencing on a NovaSeq 6000 following library preparation with Illumina Stranded Total RNA Prep, Ligation with Ribo-Zero Plus Microbiome")
df_t <- df_t %>% mutate(filetype="fastq")
df_t <- df_t %>% mutate(filename=paste("RNA_", gsub("_", "-", SampleID), "_1.fq.gz", sep=""))
df_t <- df_t %>% mutate(filename2=paste("RNA_",gsub("_", "-", SampleID), "_2.fq.gz", sep=""))
df_t <- df_t %>% mutate(filename3="")
df_t <- df_t %>% filter(!grepl("OH", SampleID))
df_t <- df_t %>% filter(!(SampleID %in% c("D02_OR_R3", "D04_OR_R1", "D05_OR_R1", "D05_OR_R2", "D05_OR_R3", "D10_OR_R1", "D10_OR_R2", "D10_OR_R3")))
df_t %>% filter(grepl("OR", SampleID)) %>% select(SampleID) %>% arrange(SampleID)

# nanopore
df_n <- read.table(here("outputs/tables/DNA_metadata.tsv"), header=TRUE, sep="\t")
df_n <- df_n %>% select(SampleID)
df_n <- df_n %>% mutate(SampleName=gsub("_R.", "", SampleID))
df_n <- df_n %>% mutate(LibraryID=paste(SampleID, "_nanopore", sep=""))
df_n <- df_n %>% mutate(LibraryLayout="single")
df_n <- df_n %>% mutate(Title=paste("Nanopore metagenomic sequencing of ", SampleID, sep=""))

df_n <- df_n %>% mutate(LibraryStrategy="WGS")
df_n <- df_n %>% mutate(LibrarySource="METAGENOMIC")
df_n <- df_n %>% mutate(LibrarySelection="RANDOM")
df_n <- df_n %>% mutate(LibrarySource="METAGENOMIC")
df_n <- df_n %>% mutate(Platform="OXFORD_NANOPORE")
df_n <- df_n %>% mutate(InstrumentModel="PromethION")
df_n <- df_n %>% mutate(DesignDescription="Nanopore metagenomic sequencing on a PromethION following library preparation with Q20+ SQK-NBD112.24 kit.")
df_n <- df_n %>% mutate(filetype="fastq")
df_n <- df_n %>% mutate(filename2="")
df_n <- df_n %>% mutate(filename3="")

df_n <- df_n %>% mutate(filename=paste(SampleID, "_nanopore.fastq.gz", sep=""))
df_n <- df_n %>% filter(SampleID %in% c("D01_NF_R1", "D02_NF_R2", "D03_NF_R1", "D04_NF_R1", "D05_NF_R1", "D06_NF_R1", "D07_NF_R1", "D08_NF_R1", "D09_NF_R1", "D10_NF_R1"))
df
# swap out D02_NF_R1 with R2

df_full <- rbind(df, df_t, df_n)
df_full <- df_full %>% select(SampleName, LibraryID, Title, LibraryStrategy, LibrarySource, LibrarySelection, LibraryLayout, Platform, InstrumentModel, DesignDescription, filetype, filename, filename2, filename3)

write.table(df_full, here("outputs/tables/SRA_Metadata.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
