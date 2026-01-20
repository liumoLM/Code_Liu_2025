# Script to export prot_table_1.csv to Excel with formatting
# - Bold header row
# - Numeric columns formatted to 4 decimal places
# - Merged cells in columns type83_sig_id through cosine_v_jin for rows with identical type83_sig_id

library(openxlsx)

# Read CSV
df <- read.csv("prot_table_1.csv", stringsAsFactors = FALSE)

# Replace underscores with spaces in column names
names(df) <- gsub("_", " ", names(df))

# Append asterisk to type83 sig id where is polyT removed is TRUE
polyT_rows <- which(df$`is polyT removed` == TRUE)
df$`type83 sig id`[polyT_rows] <- paste0(df$`type83 sig id`[polyT_rows], "*")

# Create workbook and worksheet
wb <- createWorkbook()
addWorksheet(wb, "Table 1")

# Write data starting at row 1
writeData(wb, 1, df, startRow = 1, startCol = 1)

# Style header row (bold, wrap text, centered)
headerStyle <- createStyle(
  textDecoration = "bold",
  wrapText = TRUE,
  halign = "center"
)
addStyle(wb, 1, headerStyle, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

# Identify numeric columns
numeric_cols <- which(sapply(df, is.numeric))

# Style for numeric columns (4 decimal places, centered)
numericStyle <- createStyle(numFmt = "0.0000", halign = "center")

# Style for non-numeric columns (centered)
textStyle <- createStyle(halign = "center")

# Apply styles to data rows
for (col in 1:ncol(df)) {
  if (col %in% numeric_cols) {
    addStyle(
      wb,
      1,
      numericStyle,
      rows = 2:(nrow(df) + 1),
      cols = col,
      gridExpand = TRUE
    )
  } else {
    addStyle(
      wb,
      1,
      textStyle,
      rows = 2:(nrow(df) + 1),
      cols = col,
      gridExpand = TRUE
    )
  }
}

# Calculate column widths based on data only (rows 2-48), not headers
col_widths <- sapply(1:ncol(df), function(col) {
  if (col %in% numeric_cols) {
    # For numeric columns, format to 4 decimal places and find max width
    vals <- df[[col]]
    formatted <- ifelse(is.na(vals), "NA", sprintf("%.4f", vals))
  } else {
    # For text columns, use as-is
    formatted <- as.character(df[[col]])
    formatted[is.na(formatted)] <- "NA"
  }
  max(nchar(formatted), na.rm = TRUE) + 2 # Add padding
})

setColWidths(wb, 1, cols = 1:ncol(df), widths = col_widths)

# Find column indices for merge range
merge_start_col <- which(names(df) == "type83 sig id")
merge_end_col <- which(names(df) == "cosine v jin")

# Find groups of consecutive rows with same type83 sig id
type83_col <- df$`type83 sig id`

# Get run-length encoding of type83_sig_id
rle_result <- rle(type83_col)
lengths <- rle_result$lengths
values <- rle_result$values

# Compute starting positions for each group
start_pos <- cumsum(c(1, lengths[-length(lengths)]))

# Merge cells for groups with more than 1 row
for (i in seq_along(lengths)) {
  if (lengths[i] > 1) {
    # Row indices in Excel (add 1 for header row)
    row_start <- start_pos[i] + 1
    row_end <- row_start + lengths[i] - 1

    # Merge each column in the range
    for (col in merge_start_col:merge_end_col) {
      mergeCells(wb, 1, cols = col, rows = row_start:row_end)
    }
  }
}

# Add vertical and horizontal alignment to merged cells (center both)
mergeStyle <- createStyle(valign = "center", halign = "center")
for (i in seq_along(lengths)) {
  if (lengths[i] > 1) {
    row_start <- start_pos[i] + 1
    row_end <- row_start + lengths[i] - 1
    for (col in merge_start_col:merge_end_col) {
      addStyle(
        wb,
        1,
        mergeStyle,
        rows = row_start:row_end,
        cols = col,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }

  # Freeze header row and first column
  freezePane(wb, 1, firstRow = TRUE, firstCol = TRUE)
}

# Add footnote at row 50
footnote <- "* indicates that in the supporting tumor's spectrum, insertions of a single T in long poly-T contexts were set to 0 before calculating cosine similarity to the 83-type signature."
writeData(wb, 1, footnote, startRow = 50, startCol = 1)

# Save workbook
saveWorkbook(wb, "prot_table_1.xlsx", overwrite = TRUE)

message("Saved prot_table_1.xlsx")
