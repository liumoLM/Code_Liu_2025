quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel2b_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel2b",
    types_of_interest = c("Del2:U1:R(5,9)", "Ins2:U2:R(5,9)")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel2c_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel2c",
    types_of_interest = c("Del2:U1:R(5,9)", "Ins2:U2:R(5,9)")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel7_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel7",
    types_of_interest = c("Del2:U1:R(5,9)", "Ins2:U2:R(5,9)")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel_J_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_J",
    types_of_interest = c("Del2:U1:R(5,9)", "De3:U1:R(5,9)")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel_N_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_N",
    types_of_interest = c("Del2:U1:R(5,9)", "Del3:U1:R(5,9)")
  )
)


quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel_O_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_O",
    types_of_interest = c("Ins2:U1:R(5,9)", "Ins2:U2:R(5,9)")
  )
)

quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel_P_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_P",
    types_of_interest = c("Ins2:U2:R(5,9)", "Del2:U2:R2", "Del3:U3:R2")
  )
)


quarto::quarto_render(
  input = here::here("msi_study/msi_study.qmd"),
  output_file = "InsDel_K_beta_msi_study.html",
  execute_params = list(
    sig_to_report = "InsDel_K_beta",
    types_of_interest = c()
  )
)
