# Purpose: Clustering after removing composition-associated expression
# Manuscript/panels: Supplementary Figure 2b
# Reviewer comments: R2 Comment 3
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/04_Composition_residual_clustering.R"
#!/usr/bin/env Rscript

# Reviewer 2.2: residualized expression clustering sensitivity.
# This is a sensitivity analysis, not a claim that ESTIMATE scores are true
# cellular fractions. The primary clustering procedure is kept unchanged.

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(cluster)
  library(pheatmap)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
source(file.path(dirname(dirname(script_path)), "00_paths.R"))
validate_revision_paths(c("BULK_STEP1_RDATA", "BULK_CLUSTER_RDATA"))

OUT <- file.path(RESULT_DIR, "R2_2_composition")
FIG <- file.path(FIGURE_DIR, "R2_2_composition")
LOG <- file.path(LOG_DIR, "R2_2_composition")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

load_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}
adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  a <- sum(choose2(tab)); b <- sum(choose2(rowSums(tab))); c <- sum(choose2(colSums(tab)))
  expected <- b * c / choose2(n)
  denom <- (b + c) / 2 - expected
  if (denom == 0) return(NA_real_)
  (a - expected) / denom
}
map_to_reference <- function(labels, reference) {
  ids <- names(labels)
  ref <- as.integer(sub("^C", "", reference[ids]))
  labs <- sort(unique(as.integer(labels)))
  perms <- as.matrix(expand.grid(rep(list(labs), length(labs))))
  perms <- perms[apply(perms, 1, function(z) length(unique(z)) == length(labs)), , drop = FALSE]
  agreements <- apply(perms, 1, function(p) mean(p[as.integer(labels)] == ref))
  p <- perms[which.max(agreements), ]
  setNames(paste0("C", p[as.integer(labels)]), ids)
}
consensus_metrics <- function(cc, x, k) {
  cm <- cc[[k]]$consensusMatrix
  lab <- cc[[k]]$consensusClass
  vals <- cm[lower.tri(cm)]
  d <- as.dist(1 - cor(x, method = "pearson"))
  sil <- silhouette(as.integer(lab), d)
  data.frame(k = k,
             PAC_0.1_0.9 = mean(vals > 0.1 & vals < 0.9),
             silhouette_mean = mean(sil[, "sil_width"]),
             silhouette_min = min(sil[, "sil_width"]),
             cluster_sizes = paste(as.integer(table(lab)), collapse = "/"))
}
patient_metrics <- function(cc, k, mapped) {
  cm <- cc[[k]]$consensusMatrix
  ids <- names(cc[[k]]$consensusClass)
  if (is.null(rownames(cm))) dimnames(cm) <- list(ids, ids)
  mapped <- mapped[ids]
  do.call(rbind, lapply(ids, function(id) {
    same <- names(mapped)[mapped == mapped[id] & names(mapped) != id]
    other_groups <- setdiff(unique(mapped), mapped[id])
    within <- if (length(same)) mean(cm[id, same]) else NA_real_
    between <- if (length(other_groups)) max(vapply(other_groups, function(g) {
      mean(cm[id, names(mapped)[mapped == g]])
    }, numeric(1))) else NA_real_
    data.frame(SampleID = id, Subtype = mapped[id],
               within_cluster_consensus = within,
               strongest_other_cluster_consensus = between,
               consensus_margin = within - between)
  }))
}

s1 <- load_env(BULK_STEP1_RDATA)
s2 <- load_env(BULK_CLUSTER_RDATA)
est <- read.csv(ESTIMATE_OFFICIAL_SCORES_CSV, check.names = FALSE)
est$SampleID <- sub("^tpm[.]", "", est$SampleID)
reference <- setNames(as.character(s2$meta_step2$Subtype), s2$meta_step2$SampleID)
ids <- names(reference)
est <- est[match(ids, est$SampleID), , drop = FALSE]
stopifnot(identical(est$SampleID, ids), length(ids) == 38L)

tpm <- as.matrix(s1$clean_tpm)[, ids, drop = FALSE]
expr <- log2(tpm + 1)
design <- model.matrix(~ Stromal_Score + Immune_Score,
                       data = data.frame(Stromal_Score = est$Stromal_Score,
                                         Immune_Score = est$Immune_Score))
resid_expr <- t(qr.resid(qr(design), t(expr)))
vmad <- apply(resid_expr, 1, mad, na.rm = TRUE)
genes <- names(sort(vmad[is.finite(vmad) & vmad > 0], decreasing = TRUE))[seq_len(min(5000L, sum(is.finite(vmad) & vmad > 0)))]
x <- resid_expr[genes, , drop = FALSE]
x <- sweep(x, 1, apply(x, 1, median), "-")

cc_dir <- file.path(LOG, "ccp_cache", "residualized_5000_log2_median_center_pearson_average")
dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123456)
cc <- ConsensusClusterPlus(x, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
                           clusterAlg = "hc", distance = "pearson",
                           innerLinkage = "average", finalLinkage = "average",
                           seed = 123456, title = cc_dir, plot = NULL,
                           writeTable = FALSE, verbose = FALSE)
metrics <- do.call(rbind, lapply(2:6, function(k) consensus_metrics(cc, x, k)))
write.csv(metrics, file.path(OUT, "residualized_k2_to_k6_metrics.csv"), row.names = FALSE)

mapped <- map_to_reference(cc[[3]]$consensusClass, reference)
labels <- data.frame(SampleID = ids,
                     original_Subtype = unname(reference[ids]),
                     residualized_Subtype = unname(mapped[ids]))
labels$matches_original <- labels$original_Subtype == labels$residualized_Subtype
write.csv(labels, file.path(OUT, "residualized_k3_patient_labels.csv"), row.names = FALSE)
patient <- patient_metrics(cc, 3, mapped)
write.csv(patient, file.path(OUT, "residualized_k3_patient_consensus_stability.csv"), row.names = FALSE)

summary <- data.frame(
  n_genes = length(genes),
  residualization = "log2(TPM+1) ~ official Stromal + Immune scores",
  k3_ARI_vs_original = adjusted_rand(labels$residualized_Subtype, labels$original_Subtype),
  exact_label_agreement = mean(labels$matches_original),
  changed_samples = paste(labels$SampleID[!labels$matches_original], collapse = ";"),
  k3_PAC = metrics$PAC_0.1_0.9[metrics$k == 3],
  k3_silhouette = metrics$silhouette_mean[metrics$k == 3],
  k3_cluster_sizes = metrics$cluster_sizes[metrics$k == 3],
  median_consensus_margin = median(patient$consensus_margin, na.rm = TRUE),
  min_consensus_margin = min(patient$consensus_margin, na.rm = TRUE)
)
write.csv(summary, file.path(OUT, "residualized_clustering_summary.csv"), row.names = FALSE)

cm <- cc[[3]]$consensusMatrix
cm_ids <- names(cc[[3]]$consensusClass)
if (is.null(rownames(cm))) dimnames(cm) <- list(cm_ids, cm_ids)
ord <- order(reference[cm_ids], cm_ids)
ann <- data.frame(Original_Subtype = reference[cm_ids][ord])
rownames(ann) <- cm_ids[ord]
pdf(file.path(FIG, "R2_2_residualized_consensus_matrix.pdf"), width = 7.2, height = 6.8)
pheatmap(cm[ord, ord], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = ann, annotation_col = ann,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         color = colorRampPalette(c("white", "#FDD0A2", "#D94801"))(100),
         main = "Composition-residualized consensus matrix (k = 3)")
dev.off()
png(file.path(FIG, "R2_2_residualized_consensus_matrix_preview.png"), width = 1800, height = 1700, res = 300)
pheatmap(cm[ord, ord], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = ann, annotation_col = ann,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         color = colorRampPalette(c("white", "#FDD0A2", "#D94801"))(100),
         main = "Composition-residualized consensus matrix (k = 3)")
dev.off()

writeLines(capture.output(sessionInfo()), file.path(LOG, "sessionInfo_residualized_clustering.txt"))
writeLines(c(
  paste0("n genes = ", length(genes)),
  paste0("k3 ARI vs original = ", sprintf("%.4f", summary$k3_ARI_vs_original)),
  paste0("exact label agreement = ", sprintf("%.4f", summary$exact_label_agreement)),
  paste0("changed samples = ", summary$changed_samples),
  paste0("k3 PAC = ", sprintf("%.4f", summary$k3_PAC)),
  paste0("k3 mean silhouette = ", sprintf("%.4f", summary$k3_silhouette)),
  paste0("k3 cluster sizes = ", summary$k3_cluster_sizes),
  "Interpretation: residualized clustering is a sensitivity analysis for composition-associated expression; it does not establish tumour-cell-intrinsic states."
), file.path(OUT, "residualized_clustering_verified_summary.txt"))

print(summary)
