rm(list = ls()); gc()
setwd("/public/data/Topography_analysis_2026/results0322/ID83")
if (!dir.exists("03_mergeInf")) {dir.create("03_mergeInf")}
setwd("03_mergeInf")

#This step aims to create file,including rtGroup,trans.strand,DNA.region,repli.strand,Strand;
#Meanwhile calculate the new.trans.strand and repli.strand according to Strand.

library(openxlsx)
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tibble)
options(scipen = 200)

#01.merge rtGroup,trans.strand,DNA.region,repli.strand,Strand info-------------------------
df1 <- fread("../02_intersectBed/01.RelicateTime.TransStrand.DNA_Region.intersect_res.txt",sep = "\t") %>% as.data.frame()
df3 <- fread("../02_intersectBed/02.RepliStrand.intersect_res.txt",sep = "\t") %>% as.data.frame()
id.map <- c("C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
            "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
            "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
            "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
            "ID_L", "ID_M", "ID_N")

colnames(df1) <- c("chr_rt","start_rt","end_rt","rt","group","chr_mutation","start_mutation","end_mutation","match1","Sample","ID.class","DNA.region","trans.strand",id.map)
colnames(df3) <- c("chr_repli","start_repli","end_repli", "strand.direction","length","chr_mutation","start_mutation","end_mutation","match1","Sample","ID.class","DNA.region","trans.strand",id.map)
#add match column for merge
df1 <- df1 %>%
  add_column(match = paste(df1$Sample,df1$chr_mutation,df1$start_mutation,sep = "_"), .before = 1)

dim(df1)[1] == length(unique(df1$match))
# [1] TRUE


df3 <- df3 %>%
  add_column(match = paste(df3$Sample,df3$chr_mutation,df3$start_mutation,sep = "_"), .before = 1)
df3 <- df3[,c(1:6)]
dim(df3)[1] == length(unique(df3$match))
# [1] TRUE

length(intersect(df3$match,df1$match))
#[1] 21720151 2025
# [1] 21874375  2026
# [1] 17453515 03222026

merge_res <- full_join(
  df3,        
  df1,                   
  by = "match",           
  multiple = "all"    
)
dim(merge_res)
# [1] 23195573       52

#02.add Tissue---------------------
all.meta.0528 <- readRDS("/public/data/Topography_analysis/mutation_data/all.meta.0528.rds") %>% as.data.frame()
sample_dat <- data.frame("old" = all.meta.0528$new_sample_ID,
                         "new"  = all.meta.0528$final_ID,
                         "Tissue" = all.meta.0528$cancertype)
sample_dat$old_sample <- sub(".*::", "",sample_dat$old)
sample_dat$new_sample <- sub(".*::", "",sample_dat$new)

sample_dat <- sample_dat[,c(-1,-2)]
## match new sample
sample_dat2 <- sample_dat[sample_dat$old_sample %in% unique(merge_res$Sample),]
dat1 <- sample_dat2 %>%
  left_join(merge_res,
            by = c("old_sample" = "Sample"))
head(dat1)

#03.add Strand--------------------
strand_df <- fread("/public/data/Topography_analysis/reference/single.indel.aggregate.strand.added.matrix.txt")
colnames(strand_df)

strand_df_formatch2 <- strand_df[,c("Sample","Strand","CHROM","POS")]
strand_df_formatch2$match <- paste(sub(".*::", "",strand_df_formatch2$Sample),paste0("chr",strand_df_formatch2$CHROM),strand_df_formatch2$POS,sep = "_")

length(intersect(dat1$match,strand_df_formatch2$match))

dat2 <- right_join(
  strand_df_formatch2,        
  dat1,                   
  by = "match",            
  multiple = "all"    
)
head(dat2)

#04.calculate new.trans.strand and new.repli.strand-------------------------
#The Strand information is to calculate the new +/- of trans.strand and repli.strand
#When Strand == Q ,give NA; When Strand == +, no change; When Strand == - ,get the opposite direction(+ --> -;- --> +)
#So,We gsub the Strand information(Q to 0,+ to 1,- to -1),and gsub the trans.strand.

#change the 12th colname
colnames(dat2)[12] <- "repli.strand"

#gsub the Strand inf
table(dat2$Strand)
#       -       +       Q 
# 9838346 9880429 2330980 

table(dat2$trans.strand)
#                 -        + 
# 12530386  4734635  4942331 

dat2$Strand <- case_when(
  dat2$Strand == "Q" ~ 0,
  dat2$Strand == "+" ~ 1,
  dat2$Strand == "-" ~ -1,
  TRUE ~ NA_real_  
)
dat2$trans.strand <- case_when(
  dat2$trans.strand == names(table(dat2$trans.strand))[1]  ~ 0,
  dat2$trans.strand == "+" ~ 1,
  dat2$trans.strand == "-" ~ -1,
  TRUE ~ NA_real_  
)
dat2 <- dat2 %>%
  add_column(new.trans.strand = dat2$trans.strand*dat2$Strand, .after = "Strand")
dat2 <- dat2 %>%
  add_column(new.repli.strand = dat2$repli.strand*dat2$Strand, .after = "Strand")
fwrite(dat2,"01.allInf.single.indel.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
saveRDS(dat2, file = "01.allInf.single.indel.Rds")

