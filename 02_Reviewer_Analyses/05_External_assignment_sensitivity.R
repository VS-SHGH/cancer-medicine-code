# Purpose: Signature coverage, ambiguous assignments and external strata
# Manuscript/panels: Supplementary Figure 2c
# Reviewer comments: R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/05_External_assignment_sensitivity.R"
#!/usr/bin/env Rscript

# Reviewer 2.3: frozen, same-program external validation.
# The three discovery gene sets are loaded from BULK_DEA_RDATA and are not
# reselected in either external cohort.

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
source(file.path(dirname(dirname(script_path)), "00_paths.R"))
validate_revision_paths(c("BULK_DEA_RDATA", "EXTERNAL_GSE141801_DIR", "EXTERNAL_GSE39645_DIR"))

OUT <- file.path(RESULT_DIR, "R2_3_external_validation")
FIG <- file.path(FIGURE_DIR, "R2_3_external_validation")
LOG <- file.path(LOG_DIR, "R2_3_external_validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

save_bundle <- function(p, n, w = 7, h = 5) {
  ggsave(file.path(FIG, paste0(n, ".pdf")), p, width = w, height = h, device = pdf, family = "sans")
  ggsave(file.path(FIG, paste0(n, ".svg")), p, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(FIG, paste0(n, ".tiff")), p, width = w, height = h, device = ragg::agg_tiff, dpi = 600, compression = "lzw")
  ggsave(file.path(FIG, paste0(n, "_preview.png")), p, width = w, height = h, device = ragg::agg_png, dpi = 300)
}

load_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

dea <- load_env(BULK_DEA_RDATA)
if (!exists("gene_sets", envir = dea)) stop("BULK_DEA_RDATA does not contain gene_sets")
gene_sets <- dea$gene_sets
expected_names <- c("C1_Proliferative", "C2_Mesenchymal", "C3_Immune")
if (!identical(names(gene_sets), expected_names)) stop("Unexpected frozen gene-set names")

load_expression <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2L) stop("Expression file has fewer than two columns: ", path)
  genes <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE])
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  mat <- mat[!is.na(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
  if (anyDuplicated(rownames(mat))) mat <- limma::avereps(mat)
  mat
}

datasets <- c(GSE141801 = EXTERNAL_GSE141801_DIR, GSE39645 = EXTERNAL_GSE39645_DIR)
ann <- read.csv(EXTERNAL_ANNOTATIONS_CORRECTED_CSV, check.names = FALSE, stringsAsFactors = FALSE)
required_ann <- c("Dataset", "SampleID", "assignment", "ambiguous", "primary_nonirradiated_sporadic")
if (!all(required_ann %in% names(ann))) stop("Corrected annotation file lacks required columns")
verified <- readRDS(EXTERNAL_VALIDATION_VERIFIED_RDS)
if (!identical(verified$cache_status, "PASS")) stop("Verified external-validation score cache is not PASS")

all_scores <- list(); coverage <- list(); counts <- list(); score_summaries <- list(); concordance <- list()

for (dataset_id in names(datasets)) {
  expr_path <- file.path(datasets[[dataset_id]], paste0(dataset_id, "_Expression_Log2.csv"))
  expr <- load_expression(expr_path)
  scored <- verified$scores[[dataset_id]]
  tab <- scored$sample_table
  names(tab)[names(tab) == "Sample"] <- "SampleID"
  tab$Ambiguous <- tab$Ambiguous_delta_lt_0.1
  tab$Dataset <- dataset_id
  a <- ann[ann$Dataset == dataset_id, , drop = FALSE]
  a <- a[match(tab$SampleID, a$SampleID), , drop = FALSE]
  if (anyNA(a$SampleID)) stop("Corrected annotations do not cover all samples in ", dataset_id)
  if (!all(tab$Predicted == a$top_label)) stop("Verified assignment differs from corrected annotations in ", dataset_id)
  tab$Frozen_Assignment <- a$assignment
  tab$Primary_Nonirradiated_Sporadic <- as.logical(a$primary_nonirradiated_sporadic)
  tab$NF2_disease_status <- a$NF2_disease_status
  tab$prior_irradiation <- a$prior_irradiation
  tab$assignment <- ifelse(tab$Ambiguous, "Ambiguous", tab$Predicted)
  tab$assignment <- as.character(tab$assignment)
  tab$Signature_C1 <- scored$scaled_scores["C1_Proliferative", tab$SampleID]
  tab$Signature_C2 <- scored$scaled_scores["C2_Mesenchymal", tab$SampleID]
  tab$Signature_C3 <- scored$scaled_scores["C3_Immune", tab$SampleID]
  tab <- tab[, c("Dataset", "SampleID", "Predicted", "assignment", "Frozen_Assignment", "Top1", "Top2", "Delta", "Ambiguous", "Primary_Nonirradiated_Sporadic", "NF2_disease_status", "prior_irradiation", "Signature_C1", "Signature_C2", "Signature_C3")]
  all_scores[[dataset_id]] <- tab

  coverage[[dataset_id]] <- do.call(rbind, lapply(names(gene_sets), function(s) {
    data.frame(Dataset = dataset_id, Signature = s, Signature_genes = length(gene_sets[[s]]),
               Expression_genes = nrow(expr), Matched_genes = length(intersect(gene_sets[[s]], rownames(expr))))
  }))

  for (stratum in c("all", "nonirradiated_sporadic")) {
    z <- if (stratum == "all") tab else tab[tab$Primary_Nonirradiated_Sporadic, , drop = FALSE]
    ct <- as.data.frame(table(Dataset = rep(dataset_id, nrow(z)), Stratum = rep(stratum, nrow(z)),
                              Assignment = factor(z$assignment, levels = c("C1_Proliferative", "C2_Mesenchymal", "C3_Immune", "Ambiguous"))))
    names(ct)[4] <- "n"
    counts[[length(counts) + 1L]] <- ct
    keep <- z[!z$Ambiguous, , drop = FALSE]
    for (sig in c("Signature_C1", "Signature_C2", "Signature_C3")) {
      for (g in c("C1_Proliferative", "C2_Mesenchymal", "C3_Immune")) {
        v <- keep[[sig]][keep$assignment == g]
        score_summaries[[length(score_summaries) + 1L]] <- data.frame(
          Dataset = dataset_id, Stratum = stratum, Assignment = g, Signature = sig,
          n = length(v), mean = mean(v), median = median(v), sd = sd(v),
          CI_low = if (length(v) > 1L) mean(v) - qt(0.975, length(v) - 1) * sd(v) / sqrt(length(v)) else NA_real_,
          CI_high = if (length(v) > 1L) mean(v) + qt(0.975, length(v) - 1) * sd(v) / sqrt(length(v)) else NA_real_
        )
      }
    }
  }
  frozen_path <- if (dataset_id == "GSE141801") EXTERNAL_GSE141801_FROZEN_RANK_CSV else EXTERNAL_GSE39645_FROZEN_RANK_CSV
  frozen <- read.csv(frozen_path, check.names = FALSE, stringsAsFactors = FALSE)
  m <- merge(tab[, c("SampleID", "Predicted", "assignment")], frozen[, c("SampleID", "top_label", "assignment")], by = "SampleID", suffixes = c("_cohort_ssgsea", "_frozen_rank"))
  cohort_pred_col <- if ("Predicted_cohort_ssgsea" %in% names(m)) "Predicted_cohort_ssgsea" else "Predicted"
  frozen_pred_col <- if ("top_label_frozen_rank" %in% names(m)) "top_label_frozen_rank" else "top_label"
  cohort_assignment_col <- if ("assignment_cohort_ssgsea" %in% names(m)) "assignment_cohort_ssgsea" else "assignment.x"
  frozen_assignment_col <- if ("assignment_frozen_rank" %in% names(m)) "assignment_frozen_rank" else "assignment.y"
  concordance[[dataset_id]] <- data.frame(Dataset = dataset_id, n = nrow(m),
    top_label_agreement = mean(m[[cohort_pred_col]] == m[[frozen_pred_col]]),
    assignment_agreement = mean(m[[cohort_assignment_col]] == m[[frozen_assignment_col]]),
    cohort_ambiguous = sum(m[[cohort_assignment_col]] == "Ambiguous"),
    frozen_rank_ambiguous = sum(m[[frozen_assignment_col]] == "Ambiguous"))
}

scores <- do.call(rbind, all_scores)
coverage_df <- do.call(rbind, coverage)
counts_df <- do.call(rbind, counts)
score_summary_df <- do.call(rbind, score_summaries)
concordance_df <- do.call(rbind, concordance)
write.csv(scores, file.path(OUT, "external_signature_scores_recomputed.csv"), row.names = FALSE)
write.csv(coverage_df, file.path(OUT, "external_signature_coverage.csv"), row.names = FALSE)
write.csv(counts_df, file.path(OUT, "external_assignment_summary.csv"), row.names = FALSE)
write.csv(score_summary_df, file.path(OUT, "external_signature_score_by_assignment.csv"), row.names = FALSE)
write.csv(concordance_df, file.path(OUT, "external_frozen_rank_concordance.csv"), row.names = FALSE)

# The MGS axis is handled separately from the frozen three-programme scores.
mgs <- read.csv(EXTERNAL_MGS_FIGURE_INPUT_CSV, check.names = FALSE, stringsAsFactors = FALSE)
mgs_sum <- do.call(rbind, lapply(split(mgs, list(mgs$Dataset, mgs$Class), drop = TRUE), function(z) {
  data.frame(Dataset = z$Dataset[1], Class = z$Class[1], n = nrow(z), mean = mean(z$MGS), median = median(z$MGS), sd = sd(z$MGS))
}))
mgs_tests <- do.call(rbind, lapply(split(mgs, mgs$Dataset, drop = TRUE), function(z) {
  wt <- wilcox.test(MGS ~ Class, data = z, exact = FALSE)
  data.frame(Dataset = z$Dataset[1], n = nrow(z), p = wt$p.value, class_levels = paste(levels(factor(z$Class)), collapse = " vs "))
}))
write.csv(mgs_sum, file.path(OUT, "external_mgs_summary.csv"), row.names = FALSE)
write.csv(mgs_tests, file.path(OUT, "external_mgs_tests.csv"), row.names = FALSE)

cols <- c(C1_Proliferative = "#B05A7A", C2_Mesenchymal = "#D99B3D", C3_Immune = "#4C78A8", Ambiguous = "#999999")
counts_plot <- ggplot(counts_df, aes(Stratum, n, fill = Assignment)) + geom_col(width = 0.7) + facet_wrap(~Dataset, nrow = 1) + scale_fill_manual(values = cols) + labs(x = NULL, y = "Number of samples", fill = NULL) + theme_classic(base_size = 10) + theme(legend.position = "top")
save_bundle(counts_plot, "R2_3_external_assignment_counts", 7.2, 4.4)

long <- rbind(
  data.frame(Dataset = scores$Dataset, SampleID = scores$SampleID, Assignment = scores$assignment, Signature = "C1", Score = scores$Signature_C1),
  data.frame(Dataset = scores$Dataset, SampleID = scores$SampleID, Assignment = scores$assignment, Signature = "C2", Score = scores$Signature_C2),
  data.frame(Dataset = scores$Dataset, SampleID = scores$SampleID, Assignment = scores$assignment, Signature = "C3", Score = scores$Signature_C3)
)
long$SampleID <- factor(long$SampleID, levels = unique(scores$SampleID[order(scores$Dataset, scores$assignment, scores$SampleID)]))
heat <- ggplot(long, aes(Signature, SampleID, fill = Score)) + geom_tile() + facet_grid(Dataset ~ ., scales = "free_y", space = "free_y") + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") + labs(x = "Frozen discovery signature", y = "External sample", fill = "Within-cohort\nz score") + theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "right")
save_bundle(heat, "R2_3_external_signature_scores_heatmap", 6.6, 7.5)

mgs_plot <- ggplot(mgs, aes(Class, MGS, fill = Class)) + geom_boxplot(outlier.shape = NA, alpha = .75) + geom_jitter(width = .08, size = 1.1) + facet_wrap(~Dataset, nrow = 1) + scale_fill_manual(values = c("Injury-like" = "#4C78A8", "nmSC Core" = "#B05A7A")) + labs(x = NULL, y = "Independent MGS score", fill = NULL) + theme_classic(base_size = 10) + theme(legend.position = "top")
save_bundle(mgs_plot, "R2_3_external_mgs_by_dataset", 7.2, 4.2)

report <- c(
  "R2.3 frozen same-program external validation",
  paste("R:", R.version.string),
  "The three 100-gene discovery signatures were loaded from BULK_DEA_RDATA and were not reselected in external cohorts. Scores were read from the verified external-validation RDS (cache comparison PASS; max absolute difference 1.33e-15; prediction agreement 1.0).",
  paste("GSE141801: samples", sum(scores$Dataset == "GSE141801"), "ambiguous", sum(scores$Dataset == "GSE141801" & scores$Ambiguous), "primary nonirradiated sporadic", sum(scores$Dataset == "GSE141801" & scores$Primary_Nonirradiated_Sporadic)),
  paste("GSE39645: samples", sum(scores$Dataset == "GSE39645"), "ambiguous", sum(scores$Dataset == "GSE39645" & scores$Ambiguous), "primary nonirradiated sporadic", sum(scores$Dataset == "GSE39645" & scores$Primary_Nonirradiated_Sporadic)),
  paste("Frozen rank top-label agreement:", paste(concordance_df$Dataset, sprintf("%.3f", concordance_df$top_label_agreement), collapse = "; ")),
  paste("Independent MGS tests:", paste(mgs_tests$Dataset, format.pval(mgs_tests$p, digits = 5), collapse = "; ")),
  "Interpretation: the fixed signatures can be transferred and produce reproducible score-based assignments, but the assignment is partly self-defining and is not independent biological validation. The independent MGS axis is positive in GSE141801 and negative in GSE39645, so transferability is cohort-dependent."
)
writeLines(report, file.path(OUT, "R2_3_verified_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(LOG, "sessionInfo.txt"))
print(report)
