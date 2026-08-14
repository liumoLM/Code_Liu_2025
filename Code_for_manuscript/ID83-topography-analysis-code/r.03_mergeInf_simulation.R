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
df1 <- fread("../02_intersectBed/03.RelicateTime.simulated.intersect_res.txt",sep = "\t") %>% as.data.frame()
df3 <- fread("../02_intersectBed/04.RepliStrand.simulated.intersect_res.txt",sep = "\t") %>% as.data.frame()
id.map <- c("C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
            "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
            "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
            "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
            "ID_L", "ID_M", "ID_N")


colnames(df1) <- c("chr_rt","start_rt","end_rt","rt","group","chr_mutation","start_mutation","end_mutation",
                    "sample.simu","Sample","REF","ALT","matGenClass","Strand","trans.start.pos","trans.end.pos",
                    "trans.strand","trans.gene.symbol","ID.class","DNA.region","change","indel.length","type_4",id.map)

colnames(df3) <- c("chr_repli","start_repli","end_repli", "strand.direction","length","chr_mutation","start_mutation","end_mutation",
                    "sample.simu","Sample","REF","ALT","matGenClass","Strand","trans.start.pos","trans.end.pos",
                    "trans.strand","trans.gene.symbol","ID.class","DNA.region","change","indel.length","type_4",id.map)


#add match column for merge
df1 <- df1 %>%
  add_column(match = paste(df1$Sample,df1$chr_mutation,df1$start_mutation,sep = "_"), .before = 1)
df1 <- df1[duplicated(df1$match) == FALSE,]
dim(df1)[1] == length(unique(df1$match))



df3 <- df3 %>%
  add_column(match = paste(df3$Sample,df3$chr_mutation,df3$start_mutation,sep = "_"), .before = 1)
df3 <- df3[duplicated(df3$match) == FALSE,]
df3 <- df3[,c(1:6)]
dim(df3)[1] == length(unique(df3$match))


length(intersect(df3$match,df1$match))
# [1] 194012872
# 193930220
merge_res <- full_join(
  df3,        
  df1,                   
  by = "match",            
  multiple = "all"    
)
dim(merge_res)
# [1] 240582765        63


#02.add Tissue---------------------
all.meta.0528 <- readRDS("/public/data/Topography_analysis/mutation_data/all.meta.0528.rds") %>% as.data.frame()
sample_dat <- data.frame("old" = all.meta.0528$new_sample_ID,
                         "new"  = all.meta.0528$final_ID,
                         "Tissue" = all.meta.0528$cancertype)
sample_dat$old_sample <- sub(".*::", "",sample_dat$old)
sample_dat$new_sample <- sub(".*::", "",sample_dat$new)

duplicate_check <- sample_dat %>%
  group_by(old_sample) %>%
  summarise(
    n_tissues = n_distinct(Tissue),
    tissues = ifelse(n_distinct(Tissue) > 1, 
                    paste(unique(Tissue), collapse = ", "), 
                    first(Tissue))
  ) %>%
  filter(n_tissues > 1)

if (nrow(duplicate_check) > 0) {
  cat("There are", nrow(duplicate_check), "new samples with multi tissues:\n")
  print(duplicate_check)
} else {
  cat("✓ all new samples have only one tissue\n")
}

sample_dat <- sample_dat[,c(-1,-2)]
## match new sample
sample_dat2 <- sample_dat[sample_dat$old_sample %in% unique(merge_res$Sample),]
dat1 <- sample_dat2 %>%
  left_join(merge_res,
            by = c("old_sample" = "Sample"))
head(dat1)

#04.calculate new.trans.strand and new.repli.strand-------------------------
#The Strand information is to calculate the new +/- of trans.strand and repli.strand
#When Strand == Q ,give NA; When Strand == +, no change; When Strand == - ,get the opposite direction(+ --> -;- --> +)
#So,We gsub the Strand information(Q to 0,+ to 1,- to -1),and gsub the trans.strand.

#change the 8th colname
colnames(dat1)[8] <- "repli.strand"

#gsub the Strand inf
table(dat1$Strand)
 #        -         +         Q 
 # 82458229 91370077 20655718 


table(dat1$trans.strand)
#                   -         + 
# 106330565  43015937  45137522


dat1$Strand <- case_when(
  dat1$Strand == "Q" ~ 0,
  dat1$Strand == "+" ~ 1,
  dat1$Strand == "-" ~ -1,
  TRUE ~ NA_real_  
)
dat1$trans.strand <- case_when(
  dat1$trans.strand == names(table(dat1$trans.strand))[1]  ~ 0,
  dat1$trans.strand == "+" ~ 1,
  dat1$trans.strand == "-" ~ -1,
  TRUE ~ NA_real_  
)
dat1 <- dat1 %>%
  add_column(new.trans.strand = dat1$trans.strand*dat1$Strand, .after = "Strand")
dat1 <- dat1 %>%
  add_column(new.repli.strand = dat1$repli.strand*dat1$Strand, .after = "Strand")
fwrite(dat1,"02.allInf.simulated.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
saveRDS(dat1, file = "02.allInf.simulated.Rds")

cat /public/data/Topography_analysis_2026/results/ID83/03_mergeInf/02.allInf.simulated.txt | awk '{print $5}' | sort | uniq -c | sort -k2 -V