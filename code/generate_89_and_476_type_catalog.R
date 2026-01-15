##Input
library(gridExtra)
library(ggplot2)
library(indelsig.tools.lib)
Mo_template_Koh89 <- as.data.frame(read.table(text = "
IndelType                 Indel
[Del(C):R1]A              Del(C)
[Del(C):R1]T              Del(C)
[Del(C):R2]A              Del(C)
[Del(C):R2]T              Del(C)
[Del(C):R3]A              Del(C)
[Del(C):R3]T              Del(C)
[Del(C):R(4,5)]A          Del(C)
[Del(C):R(4,5)]T          Del(C)
[Del(C):R(1,5)]G          Del(C)
Del(C):R(6,9)             Del(C)
A[Del(T):R(1,4)]A         Del(T)
A[Del(T):R(1,4)]C         Del(T)
A[Del(T):R(1,4)]G         Del(T)
C[Del(T):R(1,4)]A         Del(T)
C[Del(T):R(1,4)]C         Del(T)
C[Del(T):R(1,4)]G         Del(T)
G[Del(T):R(1,4)]A         Del(T)
G[Del(T):R(1,4)]C         Del(T)
G[Del(T):R(1,4)]G         Del(T)
A[Del(T):R(5,7)]A         Del(T)
A[Del(T):R(5,7)]C         Del(T)
A[Del(T):R(5,7)]G         Del(T)
C[Del(T):R(5,7)]A         Del(T)
C[Del(T):R(5,7)]C         Del(T)
C[Del(T):R(5,7)]G         Del(T)
G[Del(T):R(5,7)]A         Del(T)
G[Del(T):R(5,7)]C         Del(T)
G[Del(T):R(5,7)]G         Del(T)
A[Del(T):R(8,)]A          Del(T)
A[Del(T):R(8,)]C          Del(T)
A[Del(T):R(8,)]G          Del(T)
C[Del(T):R(8,)]A          Del(T)
C[Del(T):R(8,)]C          Del(T)
C[Del(T):R(8,)]G          Del(T)
G[Del(T):R(8,)]A          Del(T)
G[Del(T):R(8,)]C          Del(T)
G[Del(T):R(8,)]G          Del(T)
A[Ins(C):R0]A             Ins(C)
A[Ins(C):R0]T             Ins(C)
Ins(C):R(0,3)             Ins(C)
Ins(C):R(4,6)             Ins(C)
Ins(C):R(7,)              Ins(C)
A[Ins(T):R(0,4)]A         Ins(T)
A[Ins(T):R(0,4)]C         Ins(T)
A[Ins(T):R(0,4)]G         Ins(T)
C[Ins(T):R(0,4)]A         Ins(T)
C[Ins(T):R(0,4)]C         Ins(T)
C[Ins(T):R(0,4)]G         Ins(T)
G[Ins(T):R(0,4)]A         Ins(T)
G[Ins(T):R(0,4)]C         Ins(T)
G[Ins(T):R(0,4)]G         Ins(T)
A[Ins(T):R(5,7)]A         Ins(T)
A[Ins(T):R(5,7)]C         Ins(T)
A[Ins(T):R(5,7)]G         Ins(T)
C[Ins(T):R(5,7)]A         Ins(T)
C[Ins(T):R(5,7)]C         Ins(T)
C[Ins(T):R(5,7)]G         Ins(T)
G[Ins(T):R(5,7)]A         Ins(T)
G[Ins(T):R(5,7)]C         Ins(T)
G[Ins(T):R(5,7)]G         Ins(T)
A[Ins(T):R(8,)]A          Ins(T)
A[Ins(T):R(8,)]C          Ins(T)
A[Ins(T):R(8,)]G          Ins(T)
C[Ins(T):R(8,)]A          Ins(T)
C[Ins(T):R(8,)]C          Ins(T)
C[Ins(T):R(8,)]G          Ins(T)
G[Ins(T):R(8,)]A          Ins(T)
G[Ins(T):R(8,)]C          Ins(T)
G[Ins(T):R(8,)]G          Ins(T)
Del(2,4):R1               Del(2,):R(0,9)
Del(5,):R1                Del(2,):R(0,9)
Del(2,8):U(1,2):R(2,4)    Del(2,):R(0,9)
Del(2,):U(1,2):R(5,)      Del(2,):R(0,9)
Del(3,):U(3,):R2          Del(2,):R(0,9)
Del(3,):U(3,):R(3,)       Del(2,):R(0,9)
Ins(2,4):R0               Ins(2,)
Ins(5,):R0                Ins(2,)
Ins(2,4):R1               Ins(2,)
Ins(5,):R1                Ins(2,)
Ins(2,):R(2,4)            Ins(2,)
Ins(2,):R(5,)             Ins(2,)
Del(2,5):M1               Del(2,):M(1,)
Del(3,5):M2               Del(2,):M(1,)
Del(4,5):M(3,4)           Del(2,):M(1,)
Del(6,):M1                Del(2,):M(1,)
Del(6,):M2                Del(2,):M(1,)
Del(6,):M3                Del(2,):M(1,)
Del(6,):M(4,)             Del(2,):M(1,)
Complex                    Complex
", header = TRUE, sep = "", fill = TRUE, stringsAsFactors = FALSE))




## this function directly takes output from ICAMS::AnnovateIDVcf. 
#' @param muts_list ICAMS annotated ID vcf
#' @param sample_col  the column containing the sample name
GenerateKoh89CatalogfromAnnotateVcf <- function (muts_list, sample_col){
  muts_list <- as.data.frame(muts_list)
  indel_catalogue <- data.frame(table(muts_list[, sample_col], 
                                      muts_list$Koh_89))
  names(indel_catalogue) <- c("Sample", "IndelType", "freq")
  indel_catalogue <- reshape2::dcast(indel_catalogue, IndelType ~ 
                                       Sample, value.var = "freq")
  indel_catalogue <- merge(Mo_template_Koh89, indel_catalogue, 
                           by = "IndelType", all.x = T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  rownames(indel_catalogue) <- indel_catalogue[, "IndelType"]
  return(indel_catalogue[Mo_template_Koh89$IndelType, -c(1:2), drop = FALSE])
}



Mo_template_Koh476 <- indelsig.tools.lib::indel_template_type_4_full
Mo_template_Koh476$IndelType <- gsub("R9","R(9,)",Mo_template_Koh476$IndelType)
Mo_template_Koh476$IndelType[403] <- "Del5:U1:R(5,9)"


## this function directly takes output from ICAMS::AnnovateIDVcf. 
#' @param muts_list ICAMS annotated ID vcf
#' @param sample_col  the column containing the sample name

GenerateKoh476CatalogfromAnnotateVcf <- function (muts_list, sample_col){
  muts_list <- as.data.frame(muts_list)
  indel_catalogue <- data.frame(table(muts_list[, sample_col], 
                                      muts_list$Koh_476))
  names(indel_catalogue) <- c("Sample", "IndelType", "freq")
  indel_catalogue <- reshape2::dcast(indel_catalogue, IndelType ~ 
                                       Sample, value.var = "freq")
  indel_catalogue <- merge(Mo_template_Koh476, indel_catalogue, 
                           by = "IndelType", all.x = T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  rownames(indel_catalogue) <- indel_catalogue[, "IndelType"]
  return(indel_catalogue[Mo_template_Koh476$IndelType, -c(1:2), drop = FALSE])
}
