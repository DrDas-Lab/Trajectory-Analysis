############################################################
## 01_preprocess_counts.R
############################################################
# library(data.table)
# library(org.Hs.eg.db)

output_file <- "output/Expression_BAML_707_log2RPM.txt"

if (file.exists(output_file)) {
  
  message(sprintf("File '%s' already exists. Skipping preprocessing.", output_file))
  
} else {
  
  message("Running count preprocessing...")
  

  counts <- fread("input/beataml_waves1to4_counts_dbgap.txt", data.table = FALSE)
  rownames(counts) <- counts[, 1]
  counts_matrix <- as.matrix(counts[, sapply(counts, is.numeric)])
  
  lib_size <- colSums(counts_matrix)
  rpm <- t(t(counts_matrix) / lib_size) * 1e6
  log2rpm <- log2(rpm + 1)
  

  message("Mapping Ensembl IDs to Gene Symbols...")
  clean_ensembl <- gsub("\\..*", "", rownames(log2rpm))
  
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = clean_ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  gene_symbols[is.na(gene_symbols)] <- clean_ensembl[is.na(gene_symbols)]
  rownames(log2rpm) <- make.unique(gene_symbols)
  
  write.table(
    log2rpm,
    output_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  message(sprintf("Preprocessing complete. Processed %d genes.", nrow(log2rpm)))
}