table_1_col_name_mapping <- function(colnames) {
  mm = c(
    signature_id = "89-Type Classification Signature ID",
    sig89_from_476_v_89_cos = "Cosine vs Extracted 476-Type Sig Collapsed to 89-Type",
    exemplar_id = "Linking Tumor ID",
    sig89_v_exemplar_cos = "Cosine 89-Type 89 Sig vs Linking Tumor Spectrum",
    best_match_koh = "Best Match to Koh et al., 2025",
    cos_v_koh = "Cosine 89-Type Sig vs Koh et al.",
    sig476_v_exemplar_cos = "Cosine of Extracted 476-Type Sig vs Linking Tumor Spectrum",
    type83_sig_id = "Extracted 83-Type Sig ID",
    sig83_v_exemplar_cos = "Cosine of Extracted 83-Type Sig vs Linking Tumor Spectrum",
    best_match_cosmic = "Closest COSMIC Sig ID",
    cosine_v_cosmic = "Cosine of COSMIC Sig",
    best_match_jin = "Closest Sig from Jin et al., 2024",
    cosine_v_jin = "Cosine vs Jin et al.",
    sig83_from_476_v_83_cos = "Cosine 83-Type vs 476-Type Sig Approx Collapsed to 83-Type",
    is_polyT_removed = "is_polyT_removed"
  )

  return(mm[colnames])
}
