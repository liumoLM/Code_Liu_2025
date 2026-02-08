library(stringr)
library(ICAMS)
library(glue)

indir1 = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs"

sigFfile = "DRUP01030028T"
sig4file = "CPCT02100089T"

file1 = dir(indir1, pattern = sigFfile, full.names = TRUE) # This is the InsdelF linking Tumor

vcf <- ICAMS::ReadVCFs(file1, filter.status = "PASS")
library(ICAMS)
F <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )
sigF = F[[1]]
sigF %>% fwrite(glue("F5_study/{sigFfile}_annotvcf.csv"))

file2 = dir(indir1, pattern = sig4file, full.names = TRUE) # This is the Insdel4 linking tumor
vcf <- ICAMS::ReadVCFs(file2, filter.status = "PASS")
library(ICAMS)
F <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )
sig4 = F[[1]]
sig4 %>% fwrite(glue("F5_study/{sig4file}_annotvcf.csv"))


deltype_counts = function(annotvcf) {
  annotvcf %>%
    dplyr::filter(nchar(ins_or_del_seq) == 2) %>%
    dplyr::filter(ins_or_del == "d") %>%
    dplyr::count(short_visual) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(vis2 = gsub("\\}", "", gsub("\\{", "", short_visual))) %>%
    dplyr::mutate(dinuc = str_extract(vis2, "(?<=<).{2}(?=>)")) %>%
    dplyr::mutate(hasmh = grepl("\\{", short_visual)) -> foo
  return(foo)
}

sigFsummary = deltype_counts(sigF)
sig4summary = deltype_counts(sig4)

write.csv("foo", "del2_in_DRUP01030028T.csv")

F[[1]] %>%
  dplyr::filter(nchar(ins_or_del_seq) == 2) %>%
  dplyr::filter(ins_or_del == "d") -> zz

foo %>%
  group_by(dinuc, hasmh) %>%
  dplyr::summarize(total_n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dinuc) -> bar

View(dplyr::filter(foo, dinuc %in% c("AG", "TC")))

vcf <- ICAMS::ReadVCFs(file2, filter.status = "PASS")
library(ICAMS)
list <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )

list[[1]] %>% data.table::fwrite("")
dplyr::filter(nchar(ins_or_del_seq) > 1) %>%
  dplyr::filter(ins_or_del == "d") %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) -> foo4

write.csv("foo4", "del_in_CPCT02100089T.csv")
