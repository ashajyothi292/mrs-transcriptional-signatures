# DESeq2 - one vs rest, four tumour types
#
# NOTE: this was my first attempt and is kept for the record only. The final
# analysis uses direct pairwise comparisons (scripts/03_deseq2_pairwise_tumour_types.R).
#
# Dataset: E-MTAB-14038 (Baker-Hernandez et al., Princess Maxima Center)
# Groups: Pilocytic astrocytoma (n=107), Medulloblastoma (n=39),
#         DMG H3K27-altered (n=26), Ependymoma (n=23)
#
# I'm only using this subset of E-MTAB-14038 (the four tumour types I need),
# not the full ~1700-sample pan-cancer cohort the dataset actually contains.
# Each type gets compared against the other three pooled together.
#
# Needs DESeq2 - install with:
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Run from the repository root:  Rscript scripts/exploratory/deseq2_one_vs_rest.R

library(DESeq2)

counts <- read.csv("data/raw/counts_matrix_195_fixed.csv", row.names = 1, check.names = FALSE)
design <- read.csv("data/sample_design.csv", check.names = FALSE)

# sanity check - sample order in the counts matrix has to line up with the design table
stopifnot(all(colnames(counts) == design$sample_id))

design$tumor_type <- factor(design$tumor_type)
cat("Sample counts per tumor type:\n")
print(table(design$tumor_type))

# standard DESeq2 pre-filter - drop genes that barely show up anywhere
keep <- rowSums(counts >= 10) >= 3
cat(sprintf("\nGenes before filtering: %d | after filtering: %d\n", nrow(counts), sum(keep)))
counts_filt <- counts[keep, ]

dir.create("results/exploratory/deseq2_results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/exploratory/rnk_files", recursive = TRUE, showWarnings = FALSE)

tumor_types <- levels(design$tumor_type)

for (tt in tumor_types) {

  cat(sprintf("\n==== %s vs rest ====\n", tt))

  grp <- ifelse(design$tumor_type == tt, tt, "Other")
  grp <- factor(grp, levels = c("Other", tt))
  coldata <- data.frame(sample_id = design$sample_id, group = grp)
  rownames(coldata) <- coldata$sample_id

  dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                                 colData = coldata,
                                 design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", tt, "Other"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$stat), ]
  res_df <- res_df[order(-res_df$stat), ]  # rank by Wald stat, highest (most up) first

  safe_name <- gsub("[^A-Za-z0-9]", "_", tt)

  # full results table
  write.csv(res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
            file.path("results/exploratory/deseq2_results", paste0(safe_name, "_vs_rest_DESeq2.csv")),
            row.names = FALSE)

  # ranked list for GSEA Preranked - gene + Wald stat, that's it, no header
  rnk <- res_df[, c("gene", "stat")]
  write.table(rnk,
              file.path("results/exploratory/rnk_files", paste0(safe_name, "_vs_rest.rnk")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cat(sprintf("Significant genes (padj < 0.05): %d\n", sum(res_df$padj < 0.05, na.rm = TRUE)))
}

cat("\nDone - load the .rnk files straight into GSEA Preranked.\n")
cat("Using the DESeq2 Wald statistic as the ranking metric (standard for DESeq2 -> GSEA).\n")
