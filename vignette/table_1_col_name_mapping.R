table_1_col_name_mapping <- function(colnames) {
  mm = c(
    signature_id = "Type 89 Signature ID",
    sig89_from_476_v_89_cos = "Cosine vs<br>Extracted Type 476 Sig<br>Collapsed<br>to Type 89",
    exemplar_id = "Linking<br>Tumor ID",
    sig89_v_exemplar_cos = "Cosine vs<br>Linking Tumor<br>Spectrum",
    best_match_koh = "Best Match to<br>Koh et al., 2025",
    cos_v_koh = "Cosine vs<br>Koh et al.",
    sig476_v_exemplar_cos = "Cosine of<br>Extracted 476-Type Sig vs<br>Linking Tumor Spectrum",
    type83_sig_id = "Extracted<br>Type-83 Sig ID",
    sig83_v_exemplar_cos = "Cosine of<br>Extracted 83-Type Sig vs<br>Linking Tumor Spectrum",
    best_match_cosmic = "Closest COSMIC Sig ID",
    cosine_v_cosmic = "Cosine of<br>COSMIC Sig",
    best_match_jin = "Closest Sig from<br>Jin et al., 2024",
    cosine_v_jin = "Cosine of<br>Jin et al.",
    sig83_from_476_v_83_cos = "Cosine vs<br>Type 476 Sig Approx<br>Collapsed to Type 83",
    is_polyT_removed = "is_polyT_removed"
  )

  return(mm[colnames])
}
