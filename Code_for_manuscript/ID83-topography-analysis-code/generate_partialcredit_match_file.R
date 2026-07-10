#!/usr/bin/env Rscript
rm(list = ls()); gc()
library(data.table)
library(dplyr)
library(tibble)
ID83.vcf.file = "/public/data/Topography_analysis_2026/results/ID83/01_createBed/02.unique.single.indel.aggregate.strand.added.matrix.txt"
## generate match reference for ID83
ID83.vcfs <- fread(ID83.vcf.file)
ID83.vcfs$match <- paste(ID83.vcfs$Sample, ID83.vcfs$ID.class, sep="_")
## keep columns: Sample, type4, sample_type, signatures probability
ID83.formatch <- ID83.vcfs[!duplicated(ID83.vcfs$match), c(5, 6, 4,8:41)]
out.ID83 = "/public/data/Topography_analysis_2026/simulation_data/ID83_partialcredit_for_match.txt"
fwrite(ID83.formatch,out.ID83,quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")