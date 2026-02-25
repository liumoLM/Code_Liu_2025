quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel2b_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel2b",
    types_of_interest = c("\\QDel2:U1:R(5,9)\\E", "\\QIns2:U2:R(5,9)\\E")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel2c_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel2c",
    cosine_cutoff = 0.95,
    types_of_interest = c("\\QDel2:U1:R(5,9)\\E", "\\QIns2:U2:R(5,9)\\E")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel7_msi_study_sim_0.85.html",
  execute_params = list(
    sig_to_report = "InsDel7",
    cosine_cutoff = 0.85, # annotated VCF not avail, graylisted
    num_exemplars = 50,
    types_of_interest = c("\\QDel2:U1:R(5,9)\\E", "\\QIns2:U2:R(5,9)\\E")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel7_msi_study_sim_0.95.html",
  execute_params = list(
    sig_to_report = "InsDel7",
    cosine_cutoff = 0.95, # annotated VCF not avail, graylisted
    num_exemplars = 50,
    types_of_interest = c("\\QDel2:U1:R(5,9)\\E", "\\QIns2:U2:R(5,9)\\E")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel_J_msi_study_sim_0.95.html",
  execute_params = list(
    sig_to_report = "InsDel_J",
    num_exemplars = 50,
    cosine_cutoff = 0.95, # 0.989,
    types_of_interest = c("\\QDel2:U1:R(5,9)\\E", "\\QDe3:U1:R(5,9)\\E")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel_J_msi_study_sim_0.85_max_n_100.html",
  execute_params = list(
    sig_to_report = "InsDel_J",
    num_exemplars = 100,
    cosine_cutoff = 0.85, # 0.989,
    types_of_interest = c()
  )
)


quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel_N_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_N",
    cosine_cutoff = 0.987,
    types_of_interest = c("\\QDel2:U1:R(5,9)\\E", "\\QDel3:U1:R(5,9)\\E")
  )
)


quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel_O_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_O",
    cosine_cutoff = 0.93,
    types_of_interest = c("\\QIns2:U1:R(5,9)\\E", "\\QIns2:U2:R(5,9)\\E")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel_P_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_P",
    cosine_cutoff = 0.95,
    types_of_interest = c(
      "\\QIns2:U2:R(5,9)\\E",
      "\\QDel2:U2:R2\\E",
      "\\QDel3:U3:R2\\E"
    )
  )
)


quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel_K_beta_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_K_beta",
    types_of_interest = c(),
    cosine_cutoff = 0.85
  )
)


quarto::quarto_render(
  input = here::here("msi_study/clip_study.qmd"),
  output_file = "InsDel1d_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel1d",
    types_of_interest = c(),
    cosine_cutoff = 0.8
  )
)

# Render all signatures from the 89-to-83 mapping table
conn <- read.delim(
  here::here("Manuscript_data/89type_to_83type_connection.tsv")
)
sigids <- sub("_476$", "", conn$type476)
sigids <- sigids[sigids != ""]


for (indeltype in c("89", "476")) {
  for (sigid in sigids) {
    quarto::quarto_render(
      input = here::here("msi_study/clip_study.qmd"),
      output_file = paste0(sigid, "_", indeltype, "_clip_study.html"),
      execute_params = list(
        sig_to_report = sigid,
        cosine_cutoff = 0.0,
        num_exemplars = 30,
        indeltype = indeltype,
        types_of_interest = c(
          "\\QDel2:U1:R(5,9)\\E",
          "\\QIns2:U2:R(5,9)\\E",
          "\\QDel2:U2:R2\\E",
          "\\QDel3:U3:R2\\E",
          "Ins[^2].*:U.:R.",
          "Del.*M\\d",
          "Del[^2].*:U.:R."
        )
      )
    )
  }
}
