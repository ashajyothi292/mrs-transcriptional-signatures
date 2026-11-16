# Figures
#
# Draws the main figures from the result tables written by scripts 05, 07 and 08:
#   metabolite_profile_<cohort>.png   mean +/- SD ssGSEA score per group, drawn like an
#                                     MRS spectrum (after Gill et al. 2024, Figure 3)
#   gsea_bubble_plot.png              which metabolites go up / down in each comparison
#   top_metabolite_dotplots.png       per-patient scores for the best metabolite in each
#                                     comparison, with the ROC cut-off
#   lda_<cohort>.png                  PCA + LDA discriminant scores
#
# The medulloblastoma panels are skipped if scripts 04 and 07 have not been run.
# Run from the repository root:  Rscript scripts/09_figures.R

library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

metabolites <- c("Myo-inositol", "Choline", "Creatine", "NAA", "PEth", "Glx", "Lipids", "Alanine",
                 "Aspartate", "Citrate", "Glycine", "GABA", "Glutathione", "Glucose", "Lactate",
                 "Taurine", "scyllo-Ins")
tumour_types <- c("Pilocytic astrocytoma", "Medulloblastoma", "Ependymoma", "DMG_H3K27_altered")
subgroups <- c("WNT", "SHH", "Group3", "Group4")

# four colours that stay distinguishable for colour-blind readers
group_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")

my_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA))

# ---- load the score tables ----
design <- read.csv("data/sample_design.csv", check.names = FALSE)
tt_scores <- read.csv("results/ssGSEA_ROC/gsva_scores_pairwise_tumor_types.csv",
                      row.names = 1, check.names = FALSE)
tt_groups <- factor(design$tumor_type[match(rownames(tt_scores), design$sample_id)], levels = tumour_types)

have_mb <- file.exists("results/ssGSEA_ROC/gsva_scores_MB_subgroups.csv") &&
  file.exists("data/GSE85217_subgroups.csv")
if (have_mb) {
  mb_scores <- read.csv("results/ssGSEA_ROC/gsva_scores_MB_subgroups.csv", row.names = 1, check.names = FALSE)
  labels <- read.csv("data/GSE85217_subgroups.csv")
  mb_groups <- factor(labels$subgroup[match(rownames(mb_scores), labels$sample_id)], levels = subgroups)
}

# ---- 1. spectrum-style metabolite profiles ----
profile_plot <- function(scores, groups, title) {
  rows <- list()
  for (g in levels(groups)) {
    for (m in metabolites) {
      x <- scores[groups == g, m]
      rows[[length(rows) + 1]] <- data.frame(group = g, metabolite = m, mean = mean(x), sd = sd(x),
                                             panel = sprintf("%s (n = %d)", g, length(x)))
    }
  }
  d <- do.call(rbind, rows)
  d$metabolite <- factor(d$metabolite, levels = metabolites)
  d$panel <- factor(d$panel, levels = unique(d$panel))

  ggplot(d, aes(metabolite, mean, group = panel, colour = panel, fill = panel)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.3) +
    facet_wrap(~ panel, ncol = 2) +
    scale_colour_manual(values = group_cols, guide = "none") +
    scale_fill_manual(values = group_cols, guide = "none") +
    labs(title = title, subtitle = "Line = mean ssGSEA score, shaded band = SD",
         x = NULL, y = "ssGSEA pathway score") +
    my_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggsave("results/figures/metabolite_profile_tumour_types.png",
       profile_plot(tt_scores, tt_groups, "Metabolite-pathway activity by tumour type"),
       width = 9, height = 6.5, dpi = 200)
if (have_mb) {
  ggsave("results/figures/metabolite_profile_MB_subgroups.png",
         profile_plot(mb_scores, mb_groups, "Metabolite-pathway activity by medulloblastoma subgroup"),
         width = 9, height = 6.5, dpi = 200)
}

# ---- 2. GSEA bubble plot ----
gs <- read.csv("results/GSEA/GSEA_metabolite_summary.csv")
gs <- gs[gs$direction != "ns", ]
gs$block <- ifelse(grepl("_vs_rest$", gs$comparison), "Medulloblastoma subgroups", "Tumour types")
short <- gsub("Pilocytic_astrocytoma", "PA", gs$comparison)
short <- gsub("Medulloblastoma", "MB", short)
short <- gsub("Ependymoma", "EPN", short)
short <- gsub("DMG_H3K27_altered", "DMG", short)
gs$label <- gsub("_", " ", short)
label_order <- c(paste(subgroups, "vs rest"), "PA vs MB", "PA vs EPN", "PA vs DMG",
                 "MB vs EPN", "MB vs DMG", "EPN vs DMG")
gs$label <- factor(gs$label, levels = label_order)
gs$metabolite <- factor(gs$metabolite, levels = rev(metabolites))
gs$confidence <- -log10(pmax(gs$FDR, 1e-4))
lim <- max(abs(gs$NES))

bubble <- ggplot(gs, aes(label, metabolite)) +
  geom_point(aes(size = confidence, fill = NES), shape = 21, colour = "white", stroke = 0.5) +
  facet_grid(. ~ block, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                       limits = c(-lim, lim), name = "NES") +
  scale_size_continuous(range = c(2, 7), name = "-log10 FDR q") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "MRS metabolite pathway enrichment",
       subtitle = "GSEA Preranked, FDR q < 0.25. Red = higher in the first-named group, blue = lower, blank = not significant",
       x = NULL, y = NULL) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/gsea_bubble_plot.png", bubble, width = 9.5, height = 6.5, dpi = 200)

# ---- 3. dot plots for the top metabolite in each comparison ----
roc <- read.csv("results/ssGSEA_ROC/ROC_AUC_pairwise_tumor_types.csv")
if (have_mb && file.exists("results/ssGSEA_ROC/ROC_AUC_MB_subgroups.csv")) {
  roc <- rbind(read.csv("results/ssGSEA_ROC/ROC_AUC_MB_subgroups.csv"), roc)
}

dot_rows <- list()
for (comp in unique(roc$comparison)) {
  top <- roc[roc$comparison == comp, ][1, ]            # tables are already sorted by AUC
  g1 <- sub(" vs .*", "", comp)
  g2 <- sub(".* vs ", "", comp)
  if (g2 == "rest") {
    score <- mb_scores[, top$metabolite]
    side <- ifelse(mb_groups == g1, g1, "rest")
  } else {
    keep <- tt_groups %in% c(g1, g2)
    score <- tt_scores[keep, top$metabolite]
    side <- as.character(tt_groups[keep])
  }
  dot_rows[[comp]] <- data.frame(
    panel = sprintf("%s\n%s (AUC %.2f, cut-off %.2f)", gsub("_", " ", comp), top$metabolite, top$AUC, top$cutoff),
    side = gsub("_", " ", side), first = side == g1, score = score, cutoff = top$cutoff)
}
dots <- do.call(rbind, dot_rows)
dots$panel <- factor(dots$panel, levels = unique(dots$panel))

dotplot <- ggplot(dots, aes(side, score, colour = first)) +
  geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.6) +
  geom_hline(aes(yintercept = cutoff), linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  scale_colour_manual(values = c("TRUE" = group_cols[1], "FALSE" = group_cols[2]), guide = "none") +
  labs(title = "Best single metabolite pathway for each comparison",
       subtitle = "One dot per patient. Dashed line = Youden-optimal ROC cut-off",
       x = NULL, y = "ssGSEA pathway score") +
  my_theme +
  theme(strip.text = element_text(size = 8))
ggsave("results/figures/top_metabolite_dotplots.png", dotplot,
       width = 10, height = 3.2 * ceiling(length(dot_rows) / 3), dpi = 200)

# ---- 4. LDA scatter plots ----
lda_plot <- function(label, title) {
  d <- read.csv(sprintf("results/PCA_LDA/%s_lda_scores.csv", label))
  s <- read.csv("results/PCA_LDA/classifier_summary.csv")
  s <- s[s$cohort == label, ]
  d$group <- factor(gsub("_", " ", d$group), levels = gsub("_", " ", unique(d$group)))
  centres <- aggregate(cbind(LD1, LD2) ~ group, d, mean)

  ggplot(d, aes(LD1, LD2, colour = group)) +
    geom_point(size = 1.6, alpha = 0.65) +
    geom_label(data = centres, aes(label = group), colour = "black", size = 3, show.legend = FALSE) +
    scale_colour_manual(values = group_cols, name = NULL) +
    labs(title = title,
         subtitle = sprintf("Leave-one-out accuracy %.1f%% (%d of %d), %d principal components",
                            100 * s$loocv_accuracy, s$loocv_correct, s$n_samples, s$n_pcs)) +
    my_theme +
    theme(legend.position = "bottom")
}
ggsave("results/figures/lda_tumour_types.png",
       lda_plot("tumour_types", "PCA + LDA on 17 metabolite pathway scores: tumour types"),
       width = 6.5, height = 5.5, dpi = 200)
if (have_mb) {
  ggsave("results/figures/lda_MB_subgroups.png",
         lda_plot("MB_subgroups", "PCA + LDA on 17 metabolite pathway scores: medulloblastoma subgroups"),
         width = 6.5, height = 5.5, dpi = 200)
}

cat("Figures written to results/figures/\n")
