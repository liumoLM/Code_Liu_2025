library(stringr)
library(ICAMS)
library(glue)
library(dplyr)
library(data.table)

indir1 = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs"

sigFfile = "DRUP01030028T"
sig4file = "CPCT02100089T"

sfpath = glue("F5_study/{sigFfile}_annotvcf.csv")
s4path = glue("F5_study/{sig4file}_annotvcf.csv")
mousepath = "Manuscript_data/Reijn.mouse.vcfs.with.annotation.txt"
wt_cell_path = "Manuscript_data/Reijn_wt_cell_annot.vcf"
mut_cell_path = "Manuscript_data/Reijn_mut_cell_annot.vcf"

file1 = dir(indir1, pattern = sigFfile, full.names = TRUE) # This is the InsdelF linking Tumor


if (!file.exists(sfpath)) {
  vcf <- ICAMS::ReadVCFs(file1, filter.status = "PASS")
  F <-
    AnnotateIDVCF(
      vcf[[1]],
      ref.genome = "hg19",
      explain_indels = 1
    )
  sigF = F[[1]]
  sigF %>% fwrite(sfpath)
}

if (!file.exists(s4path)) {
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
  sig4 %>% fwrite(s4path)
}

fread(sfpath) -> sigF
fread(s4path) -> sig4

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


deltype_counts2 = function(annotvcf) {
  annotvcf %>%
    dplyr::filter(nchar(ins_or_del_seq) == 2) %>%
    dplyr::filter(ins_or_del == "d") %>%
    dplyr::mutate(
      xvis = regmatches(long_visual, regexpr("..\\s.*\\s..", long_visual))
    ) %>%
    dplyr::count(xvis) %>%
    arrange(desc(n))
}

sigFsummary = deltype_counts(sigF)

normit = function(xx) {
  xx %>%
    mutate(post = gsub(".*>", "", vis2)) %>%
    mutate(post = gsub("\\[", "", post)) %>%
    mutate(clean = gsub("]", "", gsub("[><\\[]", "", vis2))) %>%
    mutate(
      norm = ifelse(substr(clean, 1, 1) %in% c("C", "T"), revc(clean), clean)
    ) %>%
    relocate(norm)
}
normit(sigFsummary) -> normF

sig4summary = deltype_counts(sig4)
normit(sig4summary) -> norm4

normF %>%
  filter(grepl("AG", norm)) %>%
  select(n) %>%
  sum() %>%
  '/'(sum(normF$n))
norm4 %>%
  filter(grepl("AG", norm)) %>%
  select(n) %>%
  sum() %>%
  '/'(sum(norm4$n))

# View(deltype_counts2(sigF))

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
