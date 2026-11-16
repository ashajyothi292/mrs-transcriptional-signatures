# Pairwise (one vs one) DE for the 4 tumour types in E-MTAB-14038 (n=195)
# Pilocytic astrocytoma, Medulloblastoma, Ependymoma, DMG_H3K27_altered
#
# My supervisor suggested doing direct pairwise comparisons instead of
# one-vs-rest, closer to how Gill et al. actually compared tumour types.
# That gives 6 pairs:
#   Pilocytic astrocytoma vs Medulloblastoma
#   Pilocytic astrocytoma vs Ependymoma
#   Pilocytic astrocytoma vs DMG_H3K27_altered
#   Medulloblastoma       vs Ependymoma
#   Medulloblastoma       vs DMG_H3K27_altered
#   Ependymoma            vs DMG_H3K27_altered
#
# Each pair is run as its own DESeq2 fit on just those two groups' samples
# (fresh dispersion estimate per pair - not reusing anything from the
# one-vs-rest run). Ranking statistic is the same DESeq2 Wald "stat" as
# before; positive = up in the first-named group of the pair.
#
# I added rm() + gc() after each loop iteration because the earlier
# one-vs-rest run crashed RStudio once on memory - better safe than sorry
# with 6 separate DESeq objects.
#
# Input:  data/raw/counts_matrix_195_fixed.csv (from script 02), data/sample_design.csv
# Run from the repository root:  Rscript scripts/03_deseq2_pairwise_tumour_types.R

library(DESeq2)

counts <- read.csv("data/raw/counts_matrix_195_fixed.csv", row.names = 1, check.names = FALSE)
design_full <- read.csv("data/sample_design.csv", check.names = FALSE)

stopifnot(all(colnames(counts) == design_full$sample_id))

cat("Sample counts per tumor type:\n")
print(table(design_full$tumor_type))

dir.create("results/pairwise_DE_results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pairwise_rnk_files", recursive = TRUE, showWarnings = FALSE)

tumor_types <- c("Pilocytic astrocytoma", "Medulloblastoma", "Ependymoma", "DMG_H3K27_altered")
pairs <- combn(tumor_types, 2, simplify = FALSE)

cat(sprintf("\nRunning %d pairwise comparisons...\n", length(pairs)))

for (pr in pairs) {

  g1 <- pr[1]
  g2 <- pr[2]
  cat(sprintf("\n==== %s vs %s ====\n", g1, g2))

  keep_samples <- design_full$tumor_type %in% c(g1, g2)
  sub_design <- design_full[keep_samples, ]
  sub_counts <- counts[, sub_design$sample_id]

  # re-filter genes within just this pair's samples, not the whole cohort
  keep_genes <- rowSums(sub_counts >= 10) >= 3
  cat(sprintf("Genes before filtering: %d | after filtering: %d | n samples: %d (%s=%d, %s=%d)\n",
              nrow(sub_counts), sum(keep_genes), nrow(sub_design),
              g1, sum(sub_design$tumor_type == g1),
              g2, sum(sub_design$tumor_type == g2)))
  sub_counts_filt <- sub_counts[keep_genes, ]

  grp <- factor(sub_design$tumor_type, levels = c(g2, g1))  # g1 = "positive" side of the comparison
  coldata <- data.frame(sample_id = sub_design$sample_id, group = grp)
  rownames(coldata) <- coldata$sample_id

  dds <- DESeqDataSetFromMatrix(countData = sub_counts_filt,
                                 colData = coldata,
                                 design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", g1, g2))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$stat), ]
  res_df <- res_df[order(-res_df$stat), ]

  safe1 <- gsub("[^A-Za-z0-9]", "_", g1)
  safe2 <- gsub("[^A-Za-z0-9]", "_", g2)
  pair_name <- paste0(safe1, "_vs_", safe2)

  write.csv(res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
            file.path("results/pairwise_DE_results", paste0(pair_name, "_DESeq2.csv")),
            row.names = FALSE)

  rnk <- res_df[, c("gene", "stat")]
  write.table(rnk,
              file.path("results/pairwise_rnk_files", paste0(pair_name, ".rnk")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat(sprintf("Significant genes (padj < 0.05): %d  (positive stat = up in %s)\n",
              sum(res_df$padj < 0.05, na.rm = TRUE), g1))

  # clean up before the next pair
  rm(dds, res, res_df, rnk, grp, coldata, sub_counts, sub_counts_filt, sub_design, keep_genes, keep_samples)
  gc()
}

cat("\nAll 6 pairwise comparisons done.\n")
cat("DE tables in results/pairwise_DE_results/, ranked lists in results/pairwise_rnk_files/\n")
