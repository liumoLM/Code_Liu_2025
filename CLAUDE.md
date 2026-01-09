# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains code and data for Liu et al., a scientific paper on mutational signatures, specifically focusing on indel (insertion/deletion) signatures. The project works with three classification schemes for indel mutations:
- **ID83**: 83-type indel classification (COSMIC standard)
- **ID89**: 89-type indel classification (Koh et al.)
- **ID476**: 476-type indel classification (fine-grained classification)

## Key R Packages

This project relies on several specialized R packages:
- `ICAMS` - for mutational catalog operations and plotting
- `mSigHdp` / `hdpx` - for hierarchical Dirichlet process signature extraction
- `mSigPlot` - for plotting mutational signatures (especially `plot_89`, `plot_476`)
- `indelsig.tools.lib` - custom library for indel signature tools
- `lsa` - for cosine similarity calculations

## Common Operations

### Running Signature Extraction

**mSigHdp (R)**:
```r
# From script/mSigHdp_signature_extraction.Rmd
mSigHdp::RunHdpxParallel(
  input.catalog = catalog[,-1],  # Remove mutation type column
  seedNumber = 1234,
  K.guess = 10,
  out.dir = "output_dir",
  multi.types = TRUE,  # FALSE for ID89
  ...
)
```

**SigProfilerExtractor (Python 3.10)**:
```python
# Requires SigProfilerExtractor v1.2.1
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix", output_path, input_catalog, ...)
```

### Rendering Vignettes

```bash
# From vignette/ directory
Rscript -e "rmarkdown::render('ID89_ID83_vignette_github.Rmd')"
```

## Data Files

Key data in `Manuscript_data/`:
- `Liu_et_al_final_*_type_signatures.tsv` - Final extracted signatures (83, 89, 476 types)
- `Liu_et_al_*_type_spectra.tsv` - Mutational spectra/catalogs
- `Liu_et_al_*_type_signature_assignments.tsv` - Signature attribution to samples
- `89type_to_83type_connection.tsv` - Mapping between 89-type and 83-type signatures
- `correspondence_476_to_89_and_note.xlsx` - 476 to 89 type correspondence

Reference signatures:
- `COSMIC_v3.5_ID_GRCh37_signatures.tsv` - COSMIC indel signatures
- `Koh_signatures.tsv` - Koh et al. signatures
- `jin_2024_sup_tab_1_signatures.tsv` - Jin 2024 signatures

## Code Organization

- `code/` - Core R scripts for signature mapping and clustering
- `script/` - Signature extraction pipelines (mSigHdp, SigProfilerExtractor)
- `signature_comparisons/` - Scripts for comparing signatures across sources (cosine similarity, best matches)
- `vignette/` - R Markdown documents demonstrating signature analysis
- `test_data/` - Test catalogs for ID83, ID89, ID476
- `HEK293T/` - VCF files from HEK293T cell line experiments

## Signature Comparison Workflow

The `signature_comparisons/` directory contains tools for comparing signatures:
1. `best_matches.R` - Find best matching reference signatures via cosine similarity
2. `all_pairwise_cos.R` - Compute pairwise cosine similarities
3. `83_cluster.R` / `89_cluster.R` - Hierarchical clustering of signatures

## ID Type Conversions

To convert between classification types (e.g., 476 to 89), see:
- `code/code.to.map.476.89.R` - Uses `ID476_ID89_mapping.txt` mapping file
