library(ggplot2)
library(here)
library(dplyr)

readcounts <- read.table(here("DNA/1.preprocess/preprocessing/01_processing/readcounts.tsv"), header=TRUE, sep="\t")

readcounts <- readcounts %>% filter(!grepl("CO", Sample))
median(readcounts$raw_reads)
min(readcounts$raw_reads)
max(readcounts$raw_reads)

median(readcounts$host_removed_reads)
min(readcounts$host_removed_reads)
max(readcounts$host_removed_reads)

### RNA
readcounts <- read.table(here("RNA/01_rrna/readcounts.tsv"), header=TRUE, sep="\t")

median(readcounts$raw_reads)
min(readcounts$raw_reads)
max(readcounts$raw_reads)

median(readcounts$After_RRNA)
min(readcounts$After_RRNA)
max(readcounts$After_RRNA)
