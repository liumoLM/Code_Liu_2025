# Plan: Remove PDF link from HTML output

## Context
Quarto auto-adds a link to alternative formats (PDF) in the HTML output when multiple formats are defined in the YAML. Add `format-links: false` under the `html` format options.

## Change
`vignette/vignette.qmd` — add `format-links: false` under the `html:` section in the YAML header.

## Verification
Render and confirm no PDF link appears.
