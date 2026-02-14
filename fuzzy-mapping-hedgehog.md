# Plan: Create `filter_vcf_by_mutation_type()` function

## Context
We need a function that retrieves an annotated indel VCF for a given sample and filters it to rows matching a specific 476-type indel mutation type. This supports downstream analysis of individual mutation types within a sample.

## New file: `code/filter_vcf_by_mutation_type.R`

Create a function `filter_vcf_by_mutation_type(sample_id, mutation_type)` that:

1. **Validates `mutation_type`** against `ICAMS::catalog.row.order$ID476` — stop with an informative error if not found
2. **Calls `read_annotated_vcf(sample_id)`** (sourced from `code/read_annotated_vcf.R`) to get the VCF as a `data.table`
3. **Filters** the data.table to rows where `Koh_476 == mutation_type`
4. **Returns** the filtered `data.table`

Include roxygen2 documentation following existing conventions in the codebase.

## Verification
- Source the file and call with a known sample ID and a valid 476-type mutation type
- Confirm the returned data.table only contains rows matching the specified mutation type
