# GSEA Preranked on the ten ranked gene lists, against my 138 metabolite gene sets
#
# Ranked lists:
#   results/pairwise_rnk_files/  6 tumour-type pairs      (DESeq2 Wald statistic)
#   results/MB_rnk_files/        4 subgroup-vs-rest lists (limma moderated t)
#
# Settings are the ones my supervisor asked for in the GSEA Preranked tool:
# 1,000 permutations, gene set size 15-500, no collapsing of gene symbols,
# and FDR q < 0.25 as the cut-off (the usual GSEA threshold).
#
# This script runs the same analysis in R with the fgsea package so the whole
# thing can be re-run from the command line. fgsea and the GSEA desktop tool use
# the same enrichment score, but both are permutation based, so NES and FDR
# values agree closely rather than to the last decimal.
#
# Outputs (results/GSEA/):
#   <comparison>_gsea.csv             every gene set tested in that comparison
#   GSEA_all_comparisons.csv          the ten tables stacked
#   GSEA_metabolite_summary.csv       one row per metabolite per comparison
#   GSEA_direction_counts.csv         how many metabolites go up / down per comparison
#
# Run from the repository root:  Rscript scripts/05_gsea_preranked.R

library(fgsea)

gmt_file <- "data/gene_sets/MRS_Metabolite_GeneSets_suffix.gmt"
gene_sets <- gmtPathways(gmt_file)
cat("Gene sets loaded:", length(gene_sets), "\n")

suffix_map <- c(MI = "Myo-inositol", CHO = "Choline", CR = "Creatine", NAA = "NAA", PETH = "PEth",
                GLX = "Glx", LIP = "Lipids", ALA = "Alanine", ASP = "Aspartate", CIT = "Citrate",
                GLY = "Glycine", GABA = "GABA", GSH = "Glutathione", GLC = "Glucose",
                LAC = "Lactate", TAU = "Taurine", SI = "scyllo-Ins")

rnk_files <- c(list.files("results/MB_rnk_files", pattern = "\\.rnk$", full.names = TRUE),
               list.files("results/pairwise_rnk_files", pattern = "\\.rnk$", full.names = TRUE))
if (length(rnk_files) != 10) {
  cat("Expected 10 ranked lists, found", length(rnk_files), "- run scripts 03 and 04 first\n")
}

dir.create("results/GSEA", recursive = TRUE, showWarnings = FALSE)
all_results <- list()

for (f in rnk_files) {

  comp <- sub("\\.rnk$", "", basename(f))
  rnk <- read.delim(f, header = FALSE, col.names = c("gene", "stat"))
  ranks <- setNames(rnk$stat, rnk$gene)

  set.seed(42)   # permutation test, so fix the seed to get the same numbers each run
  res <- fgseaSimple(pathways = gene_sets, stats = ranks, nperm = 1000,
                     minSize = 15, maxSize = 500)
  res <- as.data.frame(res)
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
  res$metabolite <- suffix_map[sub(".*_", "", res$pathway)]
  res$comparison <- comp
  res <- res[order(res$padj, -abs(res$NES)), ]

  write.csv(res[, c("comparison", "metabolite", "pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")],
            file.path("results/GSEA", paste0(comp, "_gsea.csv")), row.names = FALSE)
  all_results[[comp]] <- res

  cat(sprintf("%-45s %3d sets tested, %3d with FDR < 0.25\n", comp, nrow(res),
              sum(res$padj < 0.25, na.rm = TRUE)))
}

all_res <- do.call(rbind, all_results)
write.csv(all_res[, c("comparison", "metabolite", "pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")],
          "results/GSEA/GSEA_all_comparisons.csv", row.names = FALSE)

# ---- one line per metabolite per comparison ----
# A metabolite counts as changed if at least one of its gene sets has FDR < 0.25.
# I report its most significant gene set, and the sign of that NES gives the
# direction (positive = higher in the first-named group).
summary_rows <- list()
for (comp in names(all_results)) {
  res <- all_results[[comp]]
  for (m in suffix_map) {
    r <- res[res$metabolite == m, ]
    if (nrow(r) == 0) next
    sig <- r[!is.na(r$padj) & r$padj < 0.25, ]
    top <- if (nrow(sig) > 0) sig[1, ] else r[1, ]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      comparison = comp, metabolite = m,
      sets_tested = nrow(r), sets_significant = nrow(sig),
      sets_up = sum(sig$NES > 0), sets_down = sum(sig$NES < 0),
      top_gene_set = top$pathway, NES = top$NES, FDR = top$padj,
      direction = if (nrow(sig) == 0) "ns" else if (top$NES > 0) "up" else "down")
  }
}
met_summary <- do.call(rbind, summary_rows)
write.csv(met_summary, "results/GSEA/GSEA_metabolite_summary.csv", row.names = FALSE)

counts <- as.data.frame.matrix(table(factor(met_summary$comparison, levels = names(all_results)),
                                     factor(met_summary$direction, levels = c("up", "down", "ns"))))
counts <- cbind(comparison = rownames(counts), counts)
write.csv(counts, "results/GSEA/GSEA_direction_counts.csv", row.names = FALSE)

cat("\nMetabolites up / down per comparison (first-named group relative to the second):\n")
print(counts, row.names = FALSE)
cat("\nDone. Tables are in results/GSEA/\n")
