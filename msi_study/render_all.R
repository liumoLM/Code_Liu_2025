library(quarto)

sigs <- list(
  list(
    sig = "InsDel2b",
    types = list("Del2:U1:R(5,9)")
  ),
  list(
    sig = "InsDel2c",
    types = list("Del2:U1:R(5,9)")
  ),
  list(
    sig = "InsDel7",
    types = list("Del2:U1:R(5,9)")
  ),
  list(
    sig = "InsDel_J",
    types = list("Del2:U1:R(5,9)")
  ),
  list(
    sig = "InsDel_Kbeta",
    types = list("Del2:U1:R(5,9)")
  ),
  list(
    sig = "InsDel_O",
    types = list("Ins2:U2:R(5,9)", "Ins2:U1:R(5,9)")
  ),
  list(
    sig = "InsDel_P",
    types = list("Del2:U2:R2", "Del3:U3:R2", "Ins2:U2:R(5,9)")
  )
)

for (s in sigs) {
  message("=== Rendering ", s$sig, " ===")
  quarto_render(
    input = "msi_study.qmd",
    output_file = paste0(s$sig, "_msi_study.html"),
    execute_params = list(
      sig_to_report = s$sig,
      types_of_interest = s$types
    )
  )
}
