rm(list = ls()); gc()
setwd("/public/data/Topography_analysis_2026/results0322/ID83")
if (!dir.exists("04_create4InfMatrix")) {dir.create("04_create4InfMatrix")}
setwd("04_create4InfMatrix")

#This step aims to organize the four information(rtGroup,trans.strand,repli.strand,DNA.region)
#to differennt martix,which contain the information,Tissues and Signatures proportion.

library(openxlsx)
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tibble)
options(scipen = 200)

# function
createClassifyMatrix <- function(class_type,dat1,id.map,cut_value = 0.75,dname,fname){
  print(paste(class_type,"Inf"))
  data1 <- dat1[,c("Tissue",class_type,id.map)]
  data1 <- data1 %>% filter(data1[[class_type]] != 0)
  data1 <- na.omit(data1)
  data1[,class_type] <- factor(data1[,class_type])
  data1 <- data1[order(data1[,class_type]),]
  print(paste("Raw data:",nrow(data1)))
  if(cut_value == 0){
    print("cut value = 0,no need to filter.")
    sum_result <- data1 %>%
      group_by(Tissue,data1[[class_type]]) %>%
      summarise(
        across(where(is.numeric), sum),
        .groups = "drop"
      )
    print(paste("sum_result col:",ncol(sum_result)))
    print(paste("sum_result row:",nrow(sum_result)))
  }else{      #cut value != 0 
    print(paste("cut value = ",cut_value,"Filter Start"))
    binary_mat <- (data1[-(1:2)] > cut_value) * 1L
    data1 <- cbind(data1[,1:2],as.data.frame(binary_mat))
    sum_result <- data1 %>%
      group_by(Tissue,data1[[class_type]]) %>%
      summarise(
        across(where(is.numeric), sum),
        .groups = "drop"
      )
    print(paste("sum_result col:",ncol(sum_result)))
    print(paste("sum_result row:",nrow(sum_result)))
  }
  #change class_type colname
  colnames(sum_result)[2] <- dname
  #save files
  if (!dir.exists(dname)) {dir.create(dname)}
  savefilepath <- paste0("./",dname,"/")
  xlsx_filename <- paste0(savefilepath,"tissue.",fname,".xlsx")
  Rds_filename <- paste0(savefilepath,"tissue.",fname,".Rds")
  saveRDS(sum_result, file = Rds_filename)
  write.xlsx(sum_result,xlsx_filename)
  print("--------------")
  return(sum_result)
}

#main---------------------------------------------------
#01.mutation data---------------------------------------
print("Mutation Data Start!")
#load the data
dat1 <- fread("../03_mergeInf/01.allInf.single.indel.txt",sep = "\t") %>% as.data.frame()

#get tissue count
tissue_df <- as.data.frame(table(dat1$Tissue))
write.xlsx(tissue_df,"01.Tissue_count.xlsx")

id.map <- c("C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
            "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
            "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
            "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
            "ID_L", "ID_M", "ID_N")

DNA.region.result <- createClassifyMatrix("DNA.region",dat1,id.map,0,dname = "DNA.region",fname = "DNA.region.Mutation_nocutoff")
rtGroup.result <- createClassifyMatrix("group",dat1,id.map,0,dname = "rtGroup",fname = "rtGroup.Mutation_nocutoff")
repli.strand.result <- createClassifyMatrix("new.repli.strand",dat1,id.map,0,dname = "repli.strand",fname = "repli.strand.Mutation_nocutoff")
trans.strand.result <- createClassifyMatrix("new.trans.strand",dat1,id.map,0,dname = "trans.strand",fname = "trans.strand.Mutation_nocutoff")


DNA.region.result <- createClassifyMatrix("DNA.region",dat1,id.map,0.75,dname = "DNA.region",fname = "DNA.region.Mutation_75cutoff")
rtGroup.result <- createClassifyMatrix("group",dat1,id.map,0.75,dname = "rtGroup",fname = "rtGroup.Mutation_75cutoff")
repli.strand.result <- createClassifyMatrix("new.repli.strand",dat1,id.map,0.75,dname = "repli.strand",fname = "repli.strand.Mutation_75cutoff")
trans.strand.result <- createClassifyMatrix("new.trans.strand",dat1,id.map,0.75,dname = "trans.strand",fname = "trans.strand.Mutation_75cutoff")


DNA.region.result <- createClassifyMatrix("DNA.region",dat1,id.map,0.5,dname = "DNA.region",fname = "DNA.region.Mutation_50cutoff")
rtGroup.result <- createClassifyMatrix("group",dat1,id.map,0.5,dname = "rtGroup",fname = "rtGroup.Mutation_50cutoff")
repli.strand.result <- createClassifyMatrix("new.repli.strand",dat1,id.map,0.5,dname = "repli.strand",fname = "repli.strand.Mutation_50cutoff")
trans.strand.result <- createClassifyMatrix("new.trans.strand",dat1,id.map,0.5,dname = "trans.strand",fname = "trans.strand.Mutation_50cutoff")

print("Mutation Data Done!")

#02.simulated data--------------------------------------
print("Simulated Data Start!")
#load the data
dat2 <- fread("../03_mergeInf/02.allInf.simulated.txt",sep = "\t") %>% as.data.frame()
#get tissue count
tissue_df <- as.data.frame(table(dat2$Tissue))
write.xlsx(tissue_df,"02.Saimulated_Tissue_count.xlsx")
id.map <- c("C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
            "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
            "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
            "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
            "ID_L", "ID_M", "ID_N")
DNA.region.result <- createClassifyMatrix("DNA.region",dat2,id.map,0,dname = "DNA.region",fname = "DNA.region.Simulated_nocutoff")
rtGroup.result <- createClassifyMatrix("group",dat2,id.map,0,dname = "rtGroup",fname = "rtGroup.Simulated_nocutoff")
repli.strand.result <- createClassifyMatrix("new.repli.strand",dat2,id.map,0,dname = "repli.strand",fname = "repli.strand.Simulated_nocutoff")
trans.strand.result <- createClassifyMatrix("new.trans.strand",dat2,id.map,0,dname = "trans.strand",fname = "trans.strand.Simulated_nocutoff")

DNA.region.result <- createClassifyMatrix("DNA.region",dat2,id.map,0.75,dname = "DNA.region",fname = "DNA.region.Simulated_75cutoff")
rtGroup.result <- createClassifyMatrix("group",dat2,id.map,0.75,dname = "rtGroup",fname = "rtGroup.Simulated_75cutoff")
repli.strand.result <- createClassifyMatrix("new.repli.strand",dat2,id.map,0.75,dname = "repli.strand",fname = "repli.strand.Simulated_75cutoff")
trans.strand.result <- createClassifyMatrix("new.trans.strand",dat2,id.map,0.75,dname = "trans.strand",fname = "trans.strand.Simulated_75cutoff")

DNA.region.result <- createClassifyMatrix("DNA.region",dat2,id.map,0.5,dname = "DNA.region",fname = "DNA.region.Simulated_50cutoff")
rtGroup.result <- createClassifyMatrix("group",dat2,id.map,0.5,dname = "rtGroup",fname = "rtGroup.Simulated_50cutoff")
repli.strand.result <- createClassifyMatrix("new.repli.strand",dat2,id.map,0.5,dname = "repli.strand",fname = "repli.strand.Simulated_50cutoff")
trans.strand.result <- createClassifyMatrix("new.trans.strand",dat2,id.map,0.5,dname = "trans.strand",fname = "trans.strand.Simulated_50cutoff")