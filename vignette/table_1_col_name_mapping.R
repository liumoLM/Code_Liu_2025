table_1_col_name_mapping <- function(colnames) {
  mm = c(
    signature_id = "Type 89 Signature ID",
    sig89_from_476_v_89_cos = "Cosine vs Extracted Type 476 Sig Collapsed to Type 89",
    exemplar_id = "Linking Tumor ID",
    sig89_v_exemplar_cos = "Cosine Type-89 Sig vs Linking Tumor Spectrum",
    best_match_koh = "Best Match to Koh et al., 2025",
    cos_v_koh = "Cosine Type-89 Sig vs Koh et al.",
    sig476_v_exemplar_cos = "Cosine of Extracted 476-Type Sig vs Linking Tumor Spectrum",
    type83_sig_id = "Extracted Type-83 Sig ID",
    sig83_v_exemplar_cos = "Cosine of Extracted 83-Type Sig vs Linking Tumor Spectrum",
    best_match_cosmic = "Closest COSMIC Sig ID",
    cosine_v_cosmic = "Cosine of COSMIC Sig",
    best_match_jin = "Closest Sig from Jin et al., 2024",
    cosine_v_jin = "Cosine vs Jin et al.",
    sig83_from_476_v_83_cos = "Cosine vs Type 476 Sig Approx Collapsed to Type 83",
    is_polyT_removed = "is_polyT_removed"
  )

  return(mm[colnames])
}
