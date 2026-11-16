# Hallmark sanity check for the WNT and SHH medulloblastoma subgroups
#
# Before trusting my own metabolite gene sets, I checked that the pipeline finds
# what it should find. WNT tumours are defined by WNT/beta-catenin activation and
# SHH tumours by Hedgehog activation, so with the standard MSigDB Hallmark
# collection (used exactly as supplied, no curation by me):
#   WNT vs rest  should be enriched for HALLMARK_WNT_BETA_CATENIN_SIGNALING
#   SHH vs rest  should be enriched for HALLMARK_HEDGEHOG_SIGNALING
#
# The Hallmark file is not in this repository (MSigDB needs a free login).
# Download h.all.v2024.1.Hs.symbols.gmt from https://www.gsea-msigdb.org and put
# it in data/gene_sets/ - any recent version works, just change the name below.
#
# Outputs (results/GSEA_hallmark/): full Hallmark table per subgroup, the
# leading-edge genes with their rank, and the two enrichment plots.
#
# Run from the repository root:  Rscript scripts/06_hallmark_check.R

library(fgsea)
library(ggplot2)

hallmark_file <- "data/gene_sets/h.all.v2024.1.Hs.symbols.gmt"
if (!file.exists(hallmark_file)) stop("Hallmark GMT not found: ", hallmark_file)
hallmark <- gmtPathways(hallmark_file)

dir.create("results/GSEA_hallmark", recursive = TRUE, showWarnings = FALSE)

checks <- list(c("WNT", "HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
               c("SHH", "HALLMARK_HEDGEHOG_SIGNALING"))

for (chk in checks) {

  subgroup <- chk[1]
  pathway <- chk[2]
  cat(sprintf("\n==== %s vs rest ====\n", subgroup))

  rnk <- read.delim(sprintf("results/MB_rnk_files/%s_vs_rest.rnk", subgroup),
                    header = FALSE, col.names = c("gene", "stat"))
  ranks <- setNames(rnk$stat, rnk$gene)

  set.seed(42)
  res <- fgseaSimple(pathways = hallmark, stats = ranks, nperm = 1000, minSize = 15, maxSize = 500)
  res <- res[order(res$padj, -abs(res$NES)), ]

  hit <- res[res$pathway == pathway, ]
  cat(sprintf("%s: NES = %.2f, FDR = %.4f, %d genes\n", pathway, hit$NES, hit$padj, hit$size))

  # leading-edge genes = the ones driving the enrichment; rank 1 = top of the list
  le <- unlist(hit$leadingEdge)
  le_tab <- data.frame(gene = le, rank = match(le, names(ranks)), t = ranks[le])
  le_tab <- le_tab[order(le_tab$rank), ]
  cat("Leading-edge genes:", paste(head(le_tab$gene, 10), collapse = ", "), "\n")

  out <- as.data.frame(res)
  out$leadingEdge <- sapply(out$leadingEdge, paste, collapse = ";")
  write.csv(out, sprintf("results/GSEA_hallmark/%s_vs_rest_hallmark.csv", subgroup), row.names = FALSE)
  write.csv(le_tab, sprintf("results/GSEA_hallmark/%s_vs_rest_%s_leading_edge.csv", subgroup, pathway),
            row.names = FALSE)

  p <- plotEnrichment(hallmark[[pathway]], ranks) +
    labs(title = pathway,
         subtitle = sprintf("%s vs rest: NES %.2f, FDR q = %.3g, %d genes", subgroup, hit$NES, hit$padj, hit$size),
         x = "Rank in ordered gene list", y = "Enrichment score") +
    theme_bw()
  ggsave(sprintf("results/GSEA_hallmark/%s_vs_rest_%s.png", subgroup, pathway), p,
         width = 6, height = 4, dpi = 200)
}

cat("\nDone. Tables and plots are in results/GSEA_hallmark/\n")
