# Per-sample metabolite pathway scoring (ssGSEA/GSVA) + ROC/AUC
#
# This is basically trying to replicate what Gill et al. (2024) did in
# their Table 4 ("binary classification of tumours with discriminatory
# biomarkers... ROC curves... AUC, sensitivity, specificity"), except
# instead of using the directly-measured MRS metabolite concentration,
# I'm using the transcriptomic proxy - the metabolite pathway activity
# score from ssGSEA.
#
# Uses the same 17-metabolite gene-set panel I've already built and
# checked (MRS_Metabolite_GeneSets_suffix.gmt). The 138 individual gene
# sets get collapsed down into one union gene list per metabolite (so 17
# gene lists total), and every sample gets scored against all 17 with
# ssGSEA - giving one "how active is this metabolite's pathway in this
# patient" number per sample, which is the closest transcriptomic
# equivalent I have to Gill's per-patient mM concentration.
#
# Two cohorts:
#   1. GSE85217 (medulloblastoma subgroups, microarray, n=763)
#      -> WNT / SHH / Group3 / Group4 one-vs-rest ROC/AUC
#   2. E-MTAB-14038 (4 tumour types, RNA-seq, n=195)
#      -> the 6 pairwise ROC/AUC comparisons
#
# The GSE85217 expression matrix, symbol lookup and subgroup labels are
# read from the .rds file that the limma script (04) saves, so the dataset
# doesn't have to be downloaded a second time.
#
# Outputs go to results/ssGSEA_ROC/
# Run from the repository root:  Rscript scripts/07_ssgsea_roc_analysis.R

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GSVA", quietly = TRUE)) BiocManager::install("GSVA", update = FALSE, ask = FALSE)
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR", update = FALSE, ask = FALSE)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")

library(GSVA)
library(edgeR)
library(pROC)

dir.create("results/ssGSEA_ROC", recursive = TRUE, showWarnings = FALSE)

# ---- build the 17 metabolite -> union gene list from the gmt file ----
gmt_lines <- readLines("data/gene_sets/MRS_Metabolite_GeneSets_suffix.gmt")
gs_raw <- list()
for (ln in gmt_lines) {
  parts <- strsplit(ln, "\t")[[1]]
  if (length(parts) < 3) next
  name <- parts[1]
  genes <- parts[-c(1, 2)]
  genes <- genes[genes != ""]
  suf <- sub(".*_", "", name)          # metabolite suffix is on the end of the gene set name
  gs_raw[[suf]] <- union(gs_raw[[suf]], genes)
}
suffix_map <- c(CHO = "Choline", LIP = "Lipids", CR = "Creatine", PETH = "PEth", GLX = "Glx",
                 LAC = "Lactate", GLC = "Glucose", CIT = "Citrate", ALA = "Alanine", ASP = "Aspartate",
                 GLY = "Glycine", GABA = "GABA", TAU = "Taurine", GSH = "Glutathione",
                 MI = "Myo-inositol", NAA = "NAA", SI = "scyllo-Ins")
names(gs_raw) <- suffix_map[names(gs_raw)]
gs_raw <- gs_raw[!is.na(names(gs_raw))]
cat("Metabolite gene set sizes (union of every gene set sharing that suffix):\n")
print(sapply(gs_raw, length))

run_ssgsea <- function(expr_mat, gene_sets) {
  # newer GSVA versions want the ssgseaParam() object instead of the
  # old method="ssgsea" argument - this just handles either
  if (exists("ssgseaParam")) {
    p <- ssgseaParam(expr_mat, gene_sets)
    gsva(p)
  } else {
    gsva(expr_mat, gene_sets, method = "ssgsea", ssgsea.norm = TRUE)
  }
}

roc_table <- function(score_mat, groups, comparisons) {
  out <- list()
  for (comp in comparisons) {
    g1 <- comp[1]; g2 <- comp[2]
    if (g2 == "rest") {
      grp <- ifelse(groups == g1, g1, "rest")
      idx <- rep(TRUE, length(groups))
    } else {
      idx <- groups %in% c(g1, g2)
      grp <- groups[idx]
    }
    sub_scores <- score_mat[, idx, drop = FALSE]
    for (metab in rownames(sub_scores)) {
      x <- sub_scores[metab, ]
      y <- factor(grp, levels = c(g2, g1))
      if (length(unique(y)) < 2 || any(table(y) < 3)) next   # skip if a group's basically empty
      ro <- tryCatch(roc(y, x, quiet = TRUE), error = function(e) NULL)
      if (is.null(ro)) next
      auc_val <- as.numeric(auc(ro))
      best <- coords(ro, "best", best.method = "youden",
                      ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
      out[[length(out) + 1]] <- data.frame(
        comparison = paste(g1, "vs", g2), metabolite = metab, AUC = auc_val,
        cutoff = best$threshold[1], sensitivity = best$sensitivity[1], specificity = best$specificity[1],
        n1 = sum(y == g1), n2 = sum(y == g2)
      )
    }
  }
  res <- do.call(rbind, out)
  res[order(-res$AUC), ]   # best discriminators first, easier to scan
}

# ---- 1. GSE85217 medulloblastoma subgroups ----
cat("\n==== GSE85217: building a symbol-indexed expression matrix ====\n")
mb <- readRDS("data/raw/GSE85217/GSE85217_expr_symbols_subgroup.rds")   # saved by script 04
expr <- mb$expr
symbols <- mb$symbols
subgroup <- mb$subgroup
expr_sym <- expr
rownames(expr_sym) <- symbols[rownames(expr_sym)]
expr_sym <- expr_sym[!is.na(rownames(expr_sym)) & rownames(expr_sym) != "", ]
# where a symbol has multiple probes, keep whichever probe has the most variance
ord <- order(-apply(expr_sym, 1, var))
expr_sym <- expr_sym[ord, ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]
cat("expr_sym dim:", dim(expr_sym), "\n")

cat("Running ssGSEA on GSE85217 (763 samples - this takes a few minutes)...\n")
gsva_mb <- run_ssgsea(as.matrix(expr_sym), gs_raw)
write.csv(t(gsva_mb), "results/ssGSEA_ROC/gsva_scores_MB_subgroups.csv")
cat("Wrote gsva_scores_MB_subgroups.csv\n")

mb_groups <- setNames(as.character(subgroup), colnames(expr))
mb_groups <- mb_groups[colnames(gsva_mb)]
mb_comparisons <- list(c("WNT", "rest"), c("SHH", "rest"), c("Group3", "rest"), c("Group4", "rest"))
mb_roc <- roc_table(gsva_mb, mb_groups, mb_comparisons)
write.csv(mb_roc, "results/ssGSEA_ROC/ROC_AUC_MB_subgroups.csv", row.names = FALSE)
cat("\nMB subgroup results with AUC > 0.8:\n")
print(subset(mb_roc, AUC > 0.8))

# ---- 2. E-MTAB-14038 pairwise tumour types ----
cat("\n==== E-MTAB-14038: loading counts and normalising to log-CPM ====\n")
counts <- read.csv("data/raw/counts_matrix_195_fixed.csv", row.names = 1, check.names = FALSE)
design_full <- read.csv("data/sample_design.csv", check.names = FALSE)
stopifnot(all(colnames(counts) == design_full$sample_id))

dge <- DGEList(counts = counts, group = design_full$tumor_type)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE)
cat("logcpm dim:", dim(logcpm), "\n")

cat("Running ssGSEA on E-MTAB-14038...\n")
gsva_tt <- run_ssgsea(as.matrix(logcpm), gs_raw)
write.csv(t(gsva_tt), "results/ssGSEA_ROC/gsva_scores_pairwise_tumor_types.csv")
cat("Wrote gsva_scores_pairwise_tumor_types.csv\n")

tt_groups <- setNames(design_full$tumor_type, design_full$sample_id)
tt_groups <- tt_groups[colnames(gsva_tt)]
tt_pairs <- combn(c("Pilocytic astrocytoma", "Medulloblastoma", "Ependymoma", "DMG_H3K27_altered"), 2, simplify = FALSE)
tt_roc <- roc_table(gsva_tt, tt_groups, tt_pairs)
write.csv(tt_roc, "results/ssGSEA_ROC/ROC_AUC_pairwise_tumor_types.csv", row.names = FALSE)
cat("\nPairwise tumour-type results with AUC > 0.8:\n")
print(subset(tt_roc, AUC > 0.8))

cat("\nDone. Files written to results/ssGSEA_ROC/:\n")
cat("  gsva_scores_MB_subgroups.csv, ROC_AUC_MB_subgroups.csv\n")
cat("  gsva_scores_pairwise_tumor_types.csv, ROC_AUC_pairwise_tumor_types.csv\n")
