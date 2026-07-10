rm(list = ls()); gc()
setwd("/public/data/Topography_analysis_2026/results0322/ID83")
if (!dir.exists("01_createBed")) {dir.create("01_createBed")}
setwd("01_createBed")

#This step aims to create bed files to the next intersectBed step from mutation data

library(openxlsx)
# library("R.matlab")
library(data.table)
library(dplyr)
library(tibble)
library(R.utils)
options(scipen = 200)

# #01.create replication time Group data.frame bed file-------------------------------------------------------
# data <- readMat("/public/data/Topography_analysis/mutation_data/per_base_territories_20kb.mat")
# df1 <- as.data.frame(data$W[,,1])
# ## remove NA
# df1 <- df1[-grep("NaN",df1$rt),]
# ## order rt time
# df1 <- df1[order(df1$rt),]
# # check row numbers
# n <- nrow(df1)
# k <- 10
# group_size <- floor(n / k)  # 基本组大小
# remainder <- n %% k         # 余数行数
# ## add group labels
# group_labels <- rep(1:k, each = group_size)
# if (remainder > 0) {
#   group_labels <- c(group_labels, 1:remainder)
# }
# df1$group <- paste0("group", group_labels)
# table(df1$group)
# ## create bed file
# df2 <- df1[,c(1,2,3,8,18)]
# colnames(df2) <- c('chrom', 'chromStart', 'chromEnd', 'rt','group')
# df2$chrom <- paste0("chr",df2$chrom)
# write.table(df2, "01.rtGroup.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

#02.create mutaion data bed file--------------------------------------------
# df3 <- fread("/public/data/Topography_analysis/reference/single.indel.aggregate.strand.added.matrix.txt")
df3 <- fread("/public/data/Topography_analysis/cap9.for.topography.20260322/all.indel.ID83.normalized.txt")
df3 <- df3[,-1]
df3 <- df3 %>%
  add_column(new_col = df3$POS, .after = 2)
df3 <- as.data.frame(df3)
## create match column for merge
df3 <- df3 %>%
  add_column(match = paste(sub(".*::", "",df3$Sample),df3$CHROM,df3$POS,sep = "_"), .before = 1)
df3 <- df3[duplicated(df3$match) == FALSE,]
length(unique(df3$match))
dput(colnames(df3))
# c("match", "Sample", "CHROM", "new_col", "POS", "REF", "ALT", 
# "region", "pos_shift", "trans.start.pos", "trans.end.pos", "trans.strand", 
# "trans.Ensembl.gene.ID", "trans.gene.symbol", "POS2", "bothstrand", 
# "count", "ins_or_del", "pre", "ins_or_del_seq", "post", "L", 
# "U_seq", "U", "U_seq_count_in_indel_seq", "indel_str_count_in_ref", 
# "R", "R_outside_ins_or_del_seq", "mh", "koh_mh", "unit", "unit_length", 
# "internal_rep", "internal_reps", "spacer", "spacer_length", "prime3_rep", 
# "prime3_reps", "original_reps", "COSMIC_83", "Koh_89", "Koh_476", 
# "dna.region", "ID166.class", "renamed_Koh_476", "C_ID1", "C_ID2", 
# "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
# "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
# "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
# "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
# "ID_L", "ID_M", "ID_N")


## get columns
df3 <- df3[,c("CHROM","new_col","POS","match","Sample","COSMIC_83","dna.region","trans.strand",
"C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
"C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
"C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
"ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
"ID_L", "ID_M", "ID_N")]

fwrite(df3,"02.unique.single.indel.aggregate.strand.added.matrix.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
fwrite(df3,"02.unique.single.indel.aggregate.strand.added.matrix.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

#03.get trans.strand information from vcfs.RData---------------------------------
# load("/public/data/Topography_analysis/mutation_data/all.annotate.ID.vcfs.RData")
# ## remove missing info lines
# all.annotate.ID.vcfs <- all.annotate.ID.vcfs[["annotated.vcfs"]][-which(names(all.annotate.ID.vcfs[["annotated.vcfs"]])=="CPCT02160009T")]
# #apply
# trans.strand_list <- lapply(seq_along(all.annotate.ID.vcfs), function(i) {
#   dat1 <- all.annotate.ID.vcfs[[i]]
#   if(nrow(dat1) < 17) return(NULL)
#   #dat1 <- dat1[, c(1,2,12,15,16,17,9,13)]
#   dat1$CHROM <- paste0("chr", dat1$CHROM)
#   dat1$Sample <- names(all.annotate.ID.vcfs)[i]
#   return(dat1)
# })
# trans.strand_df <- do.call(rbind, Filter(Negate(is.null), trans.strand_list))
# 
# trans.strand_df <- trans.strand_df[-grep(TRUE,trans.strand_df$bothstrand),]
# table(trans.strand_df$bothstrand)
# trans.strand_df$match <- paste(trans.strand_df$Sample,trans.strand_df$CHROM,trans.strand_df$POS,sep = "_" )
# length(unique(trans.strand_df$match))
# dim(trans.strand_df)
# ## check and remove dup
# # DUP <- trans.strand_df[duplicated(trans.strand_df$match) == TRUE,]
# # write.xlsx(DUP,"duplicate.xlsx")
# trans.strand_df <- trans.strand_df[duplicated(trans.strand_df$match) == FALSE,]
# trans.strand_df2 <- trans.strand_df[,c("match","trans.strand")]

#04.create signature.bed to intersectBed-----------------------------------
#merge step02 and step03 according to match column 
# merge_res <- right_join(
#   trans.strand_df2,        
#   df3,                   
#   by = "match",        
#   multiple = "all"   
# )
# length(unique(merge_res$match))
# ## change trans.strand col after sample col
# transINF <- merge_res$trans.strand
# merge_res <- merge_res[,-c(1,2)]
# merge_res <- merge_res %>%
#   add_column(trans.strand = transINF, .after = "dna.region")
# ## change Strand col after trans.strand col
# Strand <- merge_res$Strand
# merge_res <- merge_res[,-ncol(merge_res)]
# merge_res <- merge_res %>%
#   add_column(Strand = Strand, .after = "trans.strand") %>% as.data.frame()
merge_res <- df3
head(merge_res)
fwrite(merge_res,"03.signature.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
fwrite(merge_res,"03.signature.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

# save.image("01.createBed.RData")

