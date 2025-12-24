library(tidyverse)
library(data.table)

datadir <- "/research/groups/northcgrp/projects/northcgrp_hartwell/common/illumina/northcgrp_864520_Amplicon_premade-1/"
#datadir <- "/home/vmartin/projects/pipeline/smartseq3/input"
allfiles <- list.files(datadir, pattern = "\\.gz$", recursive=TRUE, full.names=TRUE)

df <- data.frame(fastq_file = allfiles) %>%
  mutate(
    sample =  sapply(strsplit(basename(fastq_file),"_S"), `[`, 1),
    type = case_when(
      grepl("_R1_", fastq_file) ~ "fastq1", 
      grepl("_R2_", fastq_file) ~ "fastq2", 
      TRUE ~ "other"
    )
  ) %>% 
  dplyr::filter(type != "other") %>%
  pivot_wider(names_from = type, values_from = fastq_file,values_fn = list) %>%
  unnest(dplyr::starts_with("fastq")) %>%  
  dplyr::select(sample, everything())

fwrite(df,"/home/vmartin/projects/pipeline/smartseq3/samplesheet.csv")
