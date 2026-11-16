# DE for the medulloblastoma molecular subgroups
# GSE85217 (Cavalli et al. 2017), microarray, 763 samples
# WNT / SHH / Group3 / Group4, one-vs-rest, using limma (the standard
# tool for microarray DE - same idea as the DESeq2 script but for
# array data instead of RNA-seq counts).
#
# I'm re-downloading GSE85217 fresh here via GEOquery rather than reusing
# any cached copy, because an earlier partial download I had was missing
# AXIN2 and PTCH1 - both of which I need for the marker gene sanity check
# below (WNT subgroup should show AXIN2/LEF1 up, SHH should show
# GLI1/GLI3/PTCH1 up - this is a quick check my supervisor asked me to do
# to confirm the subgroup labels and DE pipeline are behaving as expected
# before trusting the GSEA/GSVA results downstream).
#
# Run from the repository root:  Rscript scripts/04_limma_medulloblastoma_subgroups.R

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma", update = FALSE, ask = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

library(GEOquery)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)

cat("Downloading GSE85217 from GEO - this takes a few minutes...\n")
dir.create("data/raw/GSE85217", recursive = TRUE, showWarnings = FALSE)
gset <- getGEO("GSE85217", GSEMatrix = TRUE, getGPL = FALSE, destdir = "data/raw/GSE85217")
gset <- gset[[1]]

expr <- exprs(gset)
pdata <- pData(gset)

cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

# probe IDs come through as ENSG IDs with a "_at" suffix tacked on - strip it
ensg <- sub("_at$", "", rownames(expr))
rownames(expr) <- ensg

# map ENSG -> gene symbol so the results are actually readable
symbols <- mapIds(org.Hs.eg.db, keys = ensg, keytype = "ENSEMBL",
                   column = "SYMBOL", multiVals = "first")
cat("Probes mapped to a gene symbol:", sum(!is.na(symbols)), "of", length(symbols), "\n")

# GEOquery names the subgroup column something like "subgroup:ch1" -
# grep for it rather than hardcoding in case the exact name shifts
subgroup_col <- grep("subgroup", colnames(pdata), ignore.case = TRUE, value = TRUE)[1]
if (is.na(subgroup_col)) {
  cat("Couldn't find the subgroup column automatically. Columns available:\n")
  print(colnames(pdata))
  stop("Set subgroup_col manually and re-run.")
}
subgroup <- pdata[[subgroup_col]]
subgroup <- trimws(sub("^.*:", "", subgroup))
subgroup <- factor(subgroup)
cat("Subgroup counts:\n")
print(table(subgroup))
stopifnot(ncol(expr) == length(subgroup))

# keep the expression matrix, symbol lookup and subgroup labels for the
# ssGSEA / ROC script, so GSE85217 doesn't have to be downloaded twice
saveRDS(list(expr = expr, symbols = symbols, subgroup = subgroup),
        "data/raw/GSE85217/GSE85217_expr_symbols_subgroup.rds")

# small sample -> subgroup table, used by the classifier and figure scripts
write.csv(data.frame(sample_id = colnames(expr), subgroup = subgroup),
          "data/GSE85217_subgroups.csv", row.names = FALSE)

dir.create("results/MB_DE_results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/MB_rnk_files", recursive = TRUE, showWarnings = FALSE)

marker_check <- list(
  WNT = c("AXIN2", "LEF1"),
  SHH = c("GLI1", "GLI3", "PTCH1")
)

for (tt in levels(subgroup)) {

  cat(sprintf("\n==== %s vs rest ====\n", tt))

  grp <- ifelse(subgroup == tt, tt, "Other")
  grp <- factor(grp, levels = c("Other", tt))
  design <- model.matrix(~ grp)

  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "t")
  res$ensembl_id <- rownames(res)
  res$gene_symbol <- symbols[res$ensembl_id]

  # multiple probes can map to the same gene symbol - keep the one with
  # the biggest |t| per symbol so GSEA gets one ranked value per gene
  res_sym <- res[!is.na(res$gene_symbol) & res$gene_symbol != "", ]
  res_sym <- res_sym[order(-abs(res_sym$t)), ]
  res_sym <- res_sym[!duplicated(res_sym$gene_symbol), ]
  res_sym <- res_sym[order(-res_sym$t), ]

  safe_name <- gsub("[^A-Za-z0-9]", "_", tt)

  write.csv(res_sym[, c("gene_symbol", "ensembl_id", "logFC", "t", "P.Value", "adj.P.Val")],
            file.path("results/MB_DE_results", paste0(safe_name, "_vs_rest_limma.csv")),
            row.names = FALSE)

  rnk <- res_sym[, c("gene_symbol", "t")]
  write.table(rnk,
              file.path("results/MB_rnk_files", paste0(safe_name, "_vs_rest.rnk")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat(sprintf("Significant genes (adj.P.Val < 0.05): %d\n",
              sum(res_sym$adj.P.Val < 0.05, na.rm = TRUE)))

  # quick sanity check against known marker genes for WNT/SHH
  if (tt %in% names(marker_check)) {
    cat(sprintf("Marker check for %s:\n", tt))
    for (g in marker_check[[tt]]) {
      row <- res_sym[res_sym$gene_symbol == g, ]
      if (nrow(row) == 0) {
        cat(sprintf("  %s: not found in the mapped data\n", g))
      } else {
        direction <- ifelse(row$t[1] > 0, "UP", "DOWN")
        cat(sprintf("  %s: t = %.2f, logFC = %.2f, adj.P = %.2e  ->  %s\n",
                    g, row$t[1], row$logFC[1], row$adj.P.Val[1], direction))
      }
    }
  }

  rm(fit, res, res_sym, rnk, grp, design)
  gc()
}

cat("\nAll 4 subgroups done.\n")
cat("DE tables in results/MB_DE_results/, ranked lists in results/MB_rnk_files/\n")
