# PCA + LDA classifier on the 17 metabolite pathway scores
#
# Gill et al. (2024) classified tumours from their measured MRS metabolite
# concentrations. Here I ask the same question with the transcriptomic version:
# can the 17 ssGSEA metabolite-pathway scores on their own tell the tumour
# types (and the medulloblastoma subgroups) apart?
#
# Steps, the same for both cohorts:
#   1. standardise the 17 scores (mean 0, SD 1)
#   2. PCA, keep the first components that explain at least 95% of the variance
#   3. linear discriminant analysis on those components
#   4. leave-one-out cross-validation (each sample is predicted by a model
#      trained on all the other samples)
#
# Input:  the ssGSEA score tables from script 07
# Output: results/PCA_LDA/  accuracy summary, confusion matrix, LD scores, variance table
#
# Run from the repository root:  Rscript scripts/08_pca_lda_classifier.R

library(MASS)

dir.create("results/PCA_LDA", recursive = TRUE, showWarnings = FALSE)

run_pca_lda <- function(scores, groups, label) {

  cat(sprintf("\n==== %s: %d samples, %d classes ====\n", label, nrow(scores), nlevels(groups)))

  pca <- prcomp(scores, center = TRUE, scale. = TRUE)
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  n_pc <- which(cumsum(var_explained) >= 0.95)[1]
  cat(sprintf("Keeping %d PCs (%.1f%% of variance)\n", n_pc, 100 * cumsum(var_explained)[n_pc]))
  pcs <- pca$x[, 1:n_pc, drop = FALSE]

  # CV = TRUE makes lda() do the leave-one-out prediction for every sample
  cv <- lda(pcs, grouping = groups, CV = TRUE)
  correct <- sum(cv$class == groups)
  acc <- correct / length(groups)
  cat(sprintf("LOOCV accuracy: %.1f%% (%d of %d)\n", 100 * acc, correct, length(groups)))

  conf <- table(actual = groups, predicted = cv$class)
  print(conf)

  # per-class sensitivity, to see which groups get confused
  per_class <- diag(conf) / rowSums(conf)
  print(round(per_class, 3))

  # full-data fit, only used to get LD1/LD2 for plotting
  fit <- lda(pcs, grouping = groups)
  ld <- predict(fit)$x

  write.csv(data.frame(component = paste0("PC", seq_along(var_explained)),
                       variance_explained = var_explained,
                       cumulative = cumsum(var_explained)),
            sprintf("results/PCA_LDA/%s_pca_variance.csv", label), row.names = FALSE)
  write.csv(as.data.frame.matrix(conf), sprintf("results/PCA_LDA/%s_confusion_matrix.csv", label))
  write.csv(data.frame(sample_id = rownames(scores), group = groups,
                       loocv_prediction = cv$class, ld),
            sprintf("results/PCA_LDA/%s_lda_scores.csv", label), row.names = FALSE)

  data.frame(cohort = label, n_samples = length(groups), n_pcs = n_pc,
             loocv_correct = correct, loocv_accuracy = acc)
}

summary_rows <- list()

# ---- 1. E-MTAB-14038: four tumour types ----
tt_scores <- read.csv("results/ssGSEA_ROC/gsva_scores_pairwise_tumor_types.csv",
                      row.names = 1, check.names = FALSE)
design <- read.csv("data/sample_design.csv", check.names = FALSE)
tt_groups <- factor(design$tumor_type[match(rownames(tt_scores), design$sample_id)],
                    levels = c("Pilocytic astrocytoma", "Medulloblastoma", "Ependymoma", "DMG_H3K27_altered"))
stopifnot(!any(is.na(tt_groups)))
summary_rows[[1]] <- run_pca_lda(as.matrix(tt_scores), tt_groups, "tumour_types")

# ---- 2. GSE85217: four medulloblastoma subgroups ----
mb_file <- "results/ssGSEA_ROC/gsva_scores_MB_subgroups.csv"
label_file <- "data/GSE85217_subgroups.csv"          # written by script 04
if (file.exists(mb_file) && file.exists(label_file)) {
  mb_scores <- read.csv(mb_file, row.names = 1, check.names = FALSE)
  labels <- read.csv(label_file)
  mb_groups <- factor(labels$subgroup[match(rownames(mb_scores), labels$sample_id)],
                      levels = c("WNT", "SHH", "Group3", "Group4"))
  stopifnot(!any(is.na(mb_groups)))
  summary_rows[[2]] <- run_pca_lda(as.matrix(mb_scores), mb_groups, "MB_subgroups")
} else {
  cat("\nMedulloblastoma scores or subgroup labels not found - run scripts 04 and 07 first.\n")
}

write.csv(do.call(rbind, summary_rows), "results/PCA_LDA/classifier_summary.csv", row.names = FALSE)
cat("\nDone. Tables are in results/PCA_LDA/\n")
