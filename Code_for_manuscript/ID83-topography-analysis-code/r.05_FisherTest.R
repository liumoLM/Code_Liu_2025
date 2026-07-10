rm(list = ls()); gc()
setwd("/Users/yaoruntian/HuangLab/Topography_analysis/result/ID83_0322/")
# setwd("/public/data/Topography_analysis/result/ID83")
if (!dir.exists("05_FisherTest")) {dir.create("05_FisherTest")}
setwd("05_FisherTest")
library(openxlsx)
library(dplyr)
library(data.table)
library(tibble)
library(reshape2)
library(tidyr)
run_fisher_tests <- function(data, tissue_col = "Tissue", strand_col = "repli.strand", 
                             type_col = "type", id_cols) {
  all_results <- list()
  # all cancer types
  tissues <- unique(data[[tissue_col]])
  
  for(tissue in tissues) {
    tissue_data <- data %>% filter(.data[[tissue_col]] == tissue)
    for(id_col in id_cols) {
      # create contingency_table
      contingency_table <- tissue_data %>%
        mutate(value = round(.data[[id_col]])) %>%  # 四舍五入取整
        select(all_of(c(strand_col, type_col)), value) %>%
        pivot_wider(names_from = all_of(strand_col), 
                    values_from = value,
                    values_fn = sum) %>%
        column_to_rownames(var = type_col) %>%
        as.matrix()

      if(nrow(contingency_table) == 2 && ncol(contingency_table) == 2) {
        fisher_res <- fisher.test(contingency_table, alternative = "two.sided")
        result <- data.frame(
          Tissue = tissue,
          ID = id_col,
          p.value = fisher_res$p.value,
          odds_ratio = fisher_res$estimate,
          conf_low = fisher_res$conf.int[1],
          conf_high = fisher_res$conf.int[2],
          method = fisher_res$method,
          stringsAsFactors = FALSE
        )
        all_results[[length(all_results) + 1]] <- result
      } else {
        warning(paste("can't create contingency_table for", tissue, id_col))
      }
    }
  }

  bind_rows(all_results)
}

id.map <- c("C_ID1", "C_ID2", "C_ID3", "C_ID4", "C_ID5", "C_ID6", "C_ID7", "C_ID8", "C_ID9", 
            "C_ID10", "C_ID11", "C_ID12", "C_ID13", "C_ID14", "C_ID15", "C_ID16", 
            "C_ID17", "C_ID18", "C_ID19", "C_ID23", "ID_A", "ID_B", "ID_C", 
            "ID_D", "ID_E", "ID_F", "ID_G", "ID_H", "ID_I", "ID_J", "ID_K", 
            "ID_L", "ID_M", "ID_N")

#01.trans.strand----------------------------------------------------
simulate_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Simulated_50cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Mutation_50cutoff.xlsx")
simulate_dat$type <- "simulated"
mutation_dat$type <- "mutation"
dat1 <- rbind(simulate_dat,mutation_dat)
fisher_test_res <- run_fisher_tests(
  data = dat1,
  tissue_col = "Tissue",
  strand_col = "trans.strand",
  type_col = "type",
  id_cols = id.map
)
write.xlsx(fisher_test_res,"trans.strand.fisher_test_50cutoff.xlsx")
##02.DNA.region----------------------------------------------------
#load data
simulate_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Simulated_50cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Mutation_50cutoff.xlsx")
simulate_dat$type <- "simulated"
mutation_dat$type <- "mutation"
dat1 <- rbind(simulate_dat,mutation_dat)
fisher_test_res <- run_fisher_tests(
  data = dat1,
  tissue_col = "Tissue",
  strand_col = "DNA.region",
  type_col = "type",
  id_cols = id.map
)
write.xlsx(fisher_test_res,"DNA.region.fisher_test_50cutoff.xlsx")
#03.repli.strand--------------------------------------------------
#load data
simulate_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Simulated_50cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Mutation_50cutoff.xlsx")
simulate_dat$type <- "simulated"
mutation_dat$type <- "mutation"
dat1 <- rbind(simulate_dat,mutation_dat)
fisher_test_res <- run_fisher_tests(
  data = dat1,
  tissue_col = "Tissue",
  strand_col = "repli.strand",
  type_col = "type",
  id_cols = id.map
)
write.xlsx(fisher_test_res,"repli.strand.fisher_test_50cutoff.xlsx")


#75cutoff----------------------------------------------
#01.trans.strand----------------------------------------------------
simulate_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Simulated_75cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/trans.strand/tissue.trans.strand.Mutation_75cutoff.xlsx")
simulate_dat$type <- "simulated"
mutation_dat$type <- "mutation"
dat1 <- rbind(simulate_dat,mutation_dat)
fisher_test_res <- run_fisher_tests(
  data = dat1,
  tissue_col = "Tissue",
  strand_col = "trans.strand",
  type_col = "type",
  id_cols = id.map
)
write.xlsx(fisher_test_res,"trans.strand.fisher_test_75cutoff.xlsx")
##02.DNA.region----------------------------------------------------
#load data
simulate_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Simulated_75cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/DNA.region/tissue.DNA.region.Mutation_75cutoff.xlsx")
simulate_dat$type <- "simulated"
mutation_dat$type <- "mutation"
dat1 <- rbind(simulate_dat,mutation_dat)
fisher_test_res <- run_fisher_tests(
  data = dat1,
  tissue_col = "Tissue",
  strand_col = "DNA.region",
  type_col = "type",
  id_cols = id.map
)
write.xlsx(fisher_test_res,"DNA.region.fisher_test_75cutoff.xlsx")
#03.repli.strand--------------------------------------------------
#load data
simulate_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Simulated_75cutoff.xlsx")
mutation_dat <- read.xlsx("../04_create4InfMatrix/repli.strand/tissue.repli.strand.Mutation_75cutoff.xlsx")
simulate_dat$type <- "simulated"
mutation_dat$type <- "mutation"
dat1 <- rbind(simulate_dat,mutation_dat)
fisher_test_res <- run_fisher_tests(
  data = dat1,
  tissue_col = "Tissue",
  strand_col = "repli.strand",
  type_col = "type",
  id_cols = id.map
)
write.xlsx(fisher_test_res,"repli.strand.fisher_test_75cutoff.xlsx")