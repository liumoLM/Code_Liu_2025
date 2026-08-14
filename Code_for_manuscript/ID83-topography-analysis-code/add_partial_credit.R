#!/usr/bin/env Rscript
rm(list = ls()); gc()
##add partial credit for SimulatedData
library(data.table)
library(dplyr)
library(tibble)
simulated.vcfs.File = "/public/data/Topography_analysis_2026/simulation_data/all_simulated_with_annotation.txt"
simulated.vcfs <- fread(simulated.vcfs.File) 
simulated.vcfs$match <- paste(simulated.vcfs$Sample, simulated.vcfs$ID.class, sep="_")

ID83.formatch.File = "/public/data/Topography_analysis_2026/simulation_data/ID83_partialcredit_for_match.txt"
ID83.formatch <- fread(ID83.formatch.File)

ID83.annotated.simulated.vcfs <- left_join(simulated.vcfs, ID83.formatch, by = "match")

## problem 1: samples with low reconstruction were removed 
## problem 2: simulated mutations with different mutation types not in the original sample
## solution: remove all the NAs in the prbability matrix??
ID83.annotated.simulated.vcfs <- ID83.annotated.simulated.vcfs[!is.na(ID83.annotated.simulated.vcfs$C_ID1), ]
dput(colnames(ID83.annotated.simulated.vcfs))
## keep the useful info
ID83.annotated.simulated <- ID83.annotated.simulated.vcfs[,c("sample.simu", "position", "Sample.x","REF", "ALT", "matGenClass","Strand",
                                                              "trans.start.pos", "trans.end.pos","trans.strand", "trans.gene.symbol", 
                                                              "ID.class.x","dna.region", "chr", "change","indel.length", "type_4",
                                                              "C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7",
                                                              "C_ID8", "C_ID9", "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", 
                                                              "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
                                                              "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
                                                              "ID_L", "ID_M", "ID_N", "ID_O", "ID_P")]
names(ID83.annotated.simulated)[1:17] <- names(simulated.vcfs)[c(1:17)]
out.ID83.annotated.File = "/public/data/Topography_analysis_2026/simulation_data/new.ID83.simualted.partialcredit.annotated.txt"
fwrite(ID83.annotated.simulated, out.ID83.annotated.File, quote=F, col.names=TRUE, row.names=FALSE, sep="\t")

df_simulated <- fread("/public/data/Topography_analysis_2026/simulation_data/new.ID83.simualted.partialcredit.annotated.txt")
df_simulated <- df_simulated %>% select(chr, position, everything())
df_simulated <- df_simulated %>%
  add_column(end = df_simulated$position, .after = 2)
fwrite(df_simulated,"/public/data/Topography_analysis_2026/results/ID83/01_createBed/04.SimulatedData.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
fwrite(df_simulated,"/public/data/Topography_analysis_2026/results/ID83/01_createBed/04.SimulatedData.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")