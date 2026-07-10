rm(list = ls()); gc()
setwd("/Users/yaoruntian/HuangLab/Topography_analysis/result/ID83_0322/")
if (!dir.exists("06_computeOddsRatio")) {dir.create("06_computeOddsRatio")}
setwd("06_computeOddsRatio")
# The colour intensity reflects the odds ratio between the ratio of real mutations and the ratio of simulated
# mutations, where each ratio is calculated using the number of mutations in the genic regions and the
# number of mutations in the intergenic regions.
library(openxlsx)
library(dplyr)
library(tidyr)
library(reshape2)
source("/Users/yaoruntian/HuangLab/Topography_analysis/scripts/ID83/scripts/ID83/r.computeOddsRatioFunctions.R")
#01.repli.strand--------------------------------------------------------------------------------------
#load the data
simulate_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Simulated_50cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Mutation_50cutoff.xlsx")
#set type list
type_list <- c(-1,1)
new_type_list <- c("Lagging","Leading")
#set cut count value
cut_count_value <- 1000
createOddsRatioMatrix(mutation_dat,simulate_dat,type_list,new_type_list,cut_count_value,dname="repli.strand")
#02.trans.strand-------------------------------------------------------------------------------------
#load the data
simulate_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Simulated_50cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Mutation_50cutoff.xlsx")
#set type list
type_list <- c(-1,1)
new_type_list <- c("Transcribed","UnTranscribed")
#set cut count value
cut_count_value <- 1000
createOddsRatioMatrix(mutation_dat,simulate_dat,type_list,new_type_list,cut_count_value,dname="trans.strand")
#03.DNA.region---------------------------------------------------------------------------------------
#load the data
simulate_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Simulated_50cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Mutation_50cutoff.xlsx")
#set type list
type_list <- c("I","G")
new_type_list <- c("Intergenic","Genic")
#set cut count value
cut_count_value <- 1000
createOddsRatioMatrix(mutation_dat,simulate_dat,type_list,new_type_list,cut_count_value,dname="DNA.region")

# #75_cutoff--------------------------------------------------------------------
# #load the data
# simulate_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Simulated_75cutoff.xlsx")
# mutation_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Mutation_75cutoff.xlsx")
# #set type list
# type_list <- c(-1,1)
# new_type_list <- c("Lagging","Leading")
# #set cut count value
# cut_count_value <- 1000
# createOddsRatioMatrix(mutation_dat,simulate_dat,type_list,new_type_list,cut_count_value,dname="repli.strand")
# #02.trans.strand-------------------------------------------------------------------------------------
# #load the data
# simulate_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Simulated_75cutoff.xlsx")
# mutation_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Mutation_75cutoff.xlsx")
# #set type list
# type_list <- c(-1,1)
# new_type_list <- c("Transcribed","UnTranscribed")
# #set cut count value
# cut_count_value <- 1000
# createOddsRatioMatrix(mutation_dat,simulate_dat,type_list,new_type_list,cut_count_value,dname="trans.strand")
# #03.DNA.region---------------------------------------------------------------------------------------
# #load the data
# simulate_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Simulated_75cutoff.xlsx")
# mutation_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Mutation_75cutoff.xlsx")
# #set type list
# type_list <- c("I","G")
# new_type_list <- c("Intergenic","Genic")
# #set cut count value
# cut_count_value <- 1000
# createOddsRatioMatrix(mutation_dat,simulate_dat,type_list,new_type_list,cut_count_value,dname="DNA.region")
