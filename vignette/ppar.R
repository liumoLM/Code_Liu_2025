# Plotting parameters for vignette
# Changes to this file will trigger plot regeneration

plot476_base_size <- 25
plot476_label_size <- 3
plto476_simplify_labels <- FALSE

ppar <- list(
  w476 = 19,
  h476 = 6,
  w89 = 24,
  h89 = 7,
  basesize89 = 30,
  textsize89 = 5,
  w83 = 19,
  h83 = 6,
  basesize83 = 20,
  textsize83 = 6,

  # Plot scaling for the html output
  hw = "100%", # Default
  hw89 = "80%",
  hw83 = "80%",

  extra89y = 1.1,
  cosine_digits = 4
)

getp <- function(parname) {
  ret = ppar[[parname]]
  if (is.null(ret)) {
    stop("unknown parameter: ", parname)
  }
  return(ret)
}
