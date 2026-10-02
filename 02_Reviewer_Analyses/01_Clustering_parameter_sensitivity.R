# Purpose: 36-setting clustering sensitivity and consensus margins
# Manuscript/panels: Supplementary Figure 1
# Reviewer comments: R2 Comment 2
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/01_Clustering_parameter_sensitivity.R"
#!/usr/bin/env Rscript

# Reviewer 2, Major Comment 1: reproducible audit of the primary clustering.
# Inputs are read through 00_paths.R. All outputs stay in the active CM workspace.

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(cluster)
  library(ggplot2)
  library(pheatmap)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
source(file.path(dirname(dirname(script_path)), "00_paths.R"))
validate_revision_paths(c("BULK_STEP1_RDATA", "BULK_CLUSTER_RDATA"))

OUT <- file.path(RESULT_DIR, "R2_1_clustering")
FIG <- file.path(FIGURE_DIR, "R2_1_clustering")
LOG <- file.path(LOG_DIR, "R2_1_clustering")
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
  a <- sum(choose2(tab))
  b <- sum(choose2(rowSums(tab)))
  c <- sum(choose2(colSums(tab)))
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

consensus_metrics <- function(cc, expr, k) {
  cm <- cc[[k]]$consensusMatrix
  lab <- cc[[k]]$consensusClass
  vals <- cm[lower.tri(cm)]
  d <- switch(attr(cc, "distance_used"),
              pearson = as.dist(1 - cor(expr, method = "pearson")),
              spearman = as.dist(1 - cor(expr, method = "spearman")),
              euclidean = dist(t(expr)),
              as.dist(1 - cor(expr, method = "pearson")))
  sil <- silhouette(as.integer(lab), d)
  auc <- 1 - mean(vals) # integral of the empirical CDF over [0,1]
  data.frame(
    k = k,
    PAC_0.1_0.9 = mean(vals > 0.1 & vals < 0.9),
    consensus_CDF_area = auc,
    silhouette_mean = mean(sil[, "sil_width"]),
    silhouette_min = min(sil[, "sil_width"]),
    cluster_sizes = paste(as.integer(table(lab)), collapse = "/")
  )
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

step1 <- load_env(BULK_STEP1_RDATA)
step2 <- load_env(BULK_CLUSTER_RDATA)
tpm <- as.matrix(step1$clean_tpm)
reference <- setNames(as.character(step2$meta_step2$Subtype), step2$meta_step2$SampleID)
tpm <- tpm[, names(reference), drop = FALSE]
stopifnot(ncol(tpm) == 38L, !anyNA(tpm), all(tpm >= 0))

log_tpm <- log2(tpm + 1)
mads <- apply(log_tpm, 1, mad, na.rm = TRUE)
rank_normal <- apply(log_tpm, 2, function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x)))
rownames(rank_normal) <- rownames(log_tpm)

prepare_input <- function(n_genes, transform) {
  source <- if (transform == "rank_normal_z") rank_normal else log_tpm
  vmad <- apply(source, 1, mad, na.rm = TRUE)
  genes <- names(sort(vmad[is.finite(vmad) & vmad > 0], decreasing = TRUE))[seq_len(n_genes)]
  x <- source[genes, , drop = FALSE]
  if (transform == "log2_median_center") {
    x <- sweep(x, 1, apply(x, 1, median), "-")
  } else {
    s <- apply(x, 1, sd)
    x <- sweep(sweep(x, 1, rowMeans(x), "-"), 1, s, "/")
  }
  x
}

# Confirm what the saved primary object actually contains.
saved_check <- data.frame(
  evidence = c("00_config.R parameter", "saved Step2 clust_input", "saved Step2 samples",
               "saved Step2 subtype sizes", "CCP default linkage when unspecified"),
  value = c("n_top_mad_genes = 5000", nrow(step2$clust_input), ncol(step2$clust_input),
            paste(as.integer(table(reference)), collapse = "/"), "average")
)
write.csv(saved_check, file.path(OUT, "primary_input_evidence.csv"), row.names = FALSE)

# Exact baseline reconstruction check before any sensitivity analysis.
baseline_input <- prepare_input(5000L, "log2_median_center")
common <- intersect(rownames(baseline_input), rownames(step2$clust_input))
reconstruction <- data.frame(
  metric = c("saved_rows", "reconstructed_rows", "same_gene_set", "same_sample_order",
             "max_absolute_difference_on_common_entries"),
  value = c(nrow(step2$clust_input), nrow(baseline_input),
            setequal(rownames(step2$clust_input), rownames(baseline_input)),
            identical(colnames(step2$clust_input), colnames(baseline_input)),
            max(abs(step2$clust_input[common, ] - baseline_input[common, ])))
)
write.csv(reconstruction, file.path(OUT, "baseline_input_reconstruction.csv"), row.names = FALSE)

configs <- expand.grid(
  n_genes = c(2000L, 3000L, 5000L, 7500L),
  transform = c("log2_median_center", "log2_gene_z", "rank_normal_z"),
  parameter = c("pearson_average", "spearman_average", "euclidean_ward.D2"),
  stringsAsFactors = FALSE
)
# Full 36-condition grid; fixed seed and 500 resamples per condition.
configs$config_id <- sprintf("G%04d_%s_%s", configs$n_genes, configs$transform, configs$parameter)

all_summary <- list()
all_labels <- list()
baseline_cc <- NULL
baseline_patient <- NULL

for (i in seq_len(nrow(configs))) {
  z <- configs[i, ]
  x <- prepare_input(z$n_genes, z$transform)
  distance <- sub("_.*$", "", z$parameter)
  linkage <- sub("^[^_]+_", "", z$parameter)
  ccp_dir <- file.path(LOG, "ccp_cache", z$config_id)
  dir.create(ccp_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(123456)
  cc <- ConsensusClusterPlus(
    x, maxK = 6, reps = 500, pItem = 0.8, pFeature = 1,
    clusterAlg = "hc", distance = distance,
    innerLinkage = linkage, finalLinkage = linkage,
    seed = 123456, title = ccp_dir, plot = NULL,
    writeTable = FALSE, verbose = FALSE
  )
  attr(cc, "distance_used") <- distance
  met <- do.call(rbind, lapply(2:6, function(k) consensus_metrics(cc, x, k)))
  met$config_id <- z$config_id
  met$n_genes <- z$n_genes
  met$transform <- z$transform
  met$parameter <- z$parameter
  met$delta_CDF_area <- c(NA, diff(met$consensus_CDF_area))
  all_summary[[i]] <- met

  mapped <- map_to_reference(cc[[3]]$consensusClass, reference)
  all_labels[[i]] <- data.frame(config_id = z$config_id, SampleID = names(mapped),
                                reference_Subtype = reference[names(mapped)],
                                sensitivity_Subtype = mapped,
                                matches_reference = mapped == reference[names(mapped)])
  if (z$n_genes == 5000L && z$transform == "log2_median_center" &&
      z$parameter == "pearson_average") {
    baseline_cc <- cc
    baseline_patient <- patient_metrics(cc, 3, mapped)
    saveRDS(cc, file.path(OUT, "baseline_consensus_k2_to_k6.rds"), compress = FALSE)
  }
  message(sprintf("[%d/%d] %s", i, nrow(configs), z$config_id))
}

summary_df <- do.call(rbind, all_summary)
labels_df <- do.call(rbind, all_labels)
write.csv(summary_df, file.path(OUT, "clustering_k2_to_k6_all_conditions.csv"), row.names = FALSE)
write.csv(labels_df, file.path(OUT, "clustering_k3_all_patient_labels.csv"), row.names = FALSE)

ari_df <- do.call(rbind, lapply(split(labels_df, labels_df$config_id), function(z) {
  data.frame(config_id = z$config_id[1],
             ARI_vs_reference = adjusted_rand(z$sensitivity_Subtype, z$reference_Subtype),
             exact_label_agreement = mean(z$matches_reference),
             changed_samples = paste(z$SampleID[!z$matches_reference], collapse = ";"))
}))
ari_df <- merge(configs, ari_df, by = "config_id", sort = FALSE)
write.csv(ari_df, file.path(OUT, "clustering_k3_parameter_stability.csv"), row.names = FALSE)
write.csv(baseline_patient, file.path(OUT, "baseline_k3_patient_consensus_stability.csv"), row.names = FALSE)

baseline_metrics <- subset(summary_df,
  n_genes == 5000 & transform == "log2_median_center" & parameter == "pearson_average")
write.csv(baseline_metrics, file.path(OUT, "baseline_k2_to_k6_metrics.csv"), row.names = FALSE)

# Baseline consensus matrix ordered by the locked C1/C2/C3 labels.
cm <- baseline_cc[[3]]$consensusMatrix
cm_ids <- names(baseline_cc[[3]]$consensusClass)
if (is.null(rownames(cm))) dimnames(cm) <- list(cm_ids, cm_ids)
ord <- order(reference[cm_ids], cm_ids)
ann <- data.frame(Subtype = reference[cm_ids][ord])
rownames(ann) <- cm_ids[ord]
pdf(file.path(FIG, "R2_1_baseline_k3_consensus_matrix.pdf"), width = 7.2, height = 6.8)
pheatmap(cm[ord, ord], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = ann, annotation_col = ann,
         show_rownames = FALSE, show_colnames = FALSE, border_color = NA,
         color = colorRampPalette(c("white", "#CFE8F3", "#2166AC"))(100),
         main = "Consensus matrix (k = 3; 5,000 MAD genes)")
dev.off()

p_metrics <- ggplot(baseline_metrics, aes(k)) +
  geom_line(aes(y = PAC_0.1_0.9, colour = "PAC"), linewidth = 0.8) +
  geom_point(aes(y = PAC_0.1_0.9, colour = "PAC"), size = 2.5) +
  geom_line(aes(y = silhouette_mean, colour = "Mean silhouette"), linewidth = 0.8) +
  geom_point(aes(y = silhouette_mean, colour = "Mean silhouette"), size = 2.5) +
  scale_x_continuous(breaks = 2:6) +
  scale_colour_manual(values = c("PAC" = "#B2182B", "Mean silhouette" = "#2166AC")) +
  labs(x = "Number of clusters (k)", y = "Metric value", colour = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(FIG, "R2_1_baseline_PAC_silhouette.pdf"), p_metrics, width = 5.8, height = 4.2)

p_ari <- ggplot(ari_df, aes(factor(n_genes), parameter, fill = ARI_vs_reference)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI_vs_reference)), size = 3) +
  facet_wrap(~transform, ncol = 1) +
  scale_fill_gradient2(limits = c(-1, 1), low = "#B2182B", mid = "white", high = "#2166AC") +
  labs(x = "Number of selected genes", y = NULL, fill = "ARI") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
ggsave(file.path(FIG, "R2_1_k3_sensitivity_ARI.pdf"), p_ari, width = 7.4, height = 8.5)

p_patient <- ggplot(baseline_patient, aes(reorder(SampleID, consensus_margin), consensus_margin,
                                           fill = Subtype)) +
  geom_col(width = 0.8) + coord_flip() +
  labs(x = "Patient", y = "Within-cluster consensus minus strongest alternative") +
  theme_classic(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(FIG, "R2_1_patient_consensus_margin.pdf"), p_patient, width = 6.8, height = 7.0)

writeLines(capture.output(sessionInfo()), file.path(LOG, "sessionInfo.txt"))
writeLines(c(
  paste("Input:", normalizePath(BULK_STEP1_RDATA)),
  paste("Saved clustering object:", normalizePath(BULK_CLUSTER_RDATA)),
  paste("Rows/columns:", nrow(tpm), ncol(tpm)),
  paste("Sensitivity conditions:", nrow(configs)),
  "Consensus resamples per condition: 500",
  "Baseline: 5,000 top-MAD genes; log2(TPM+1); gene median-centering; Pearson distance; average linkage."
), file.path(LOG, "run_summary.txt"))

message("Completed R2.1 clustering robustness analysis: ", OUT)
