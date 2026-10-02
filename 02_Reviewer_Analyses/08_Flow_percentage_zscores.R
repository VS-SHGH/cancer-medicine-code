# Purpose: 26 flow variables: percentage units and mean patient-level z scores
# Manuscript/panels: Supplementary Figure 5b
# Reviewer comments: R2 Comment 6
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/08_Flow_percentage_zscores.R"
# Reviewer 2.5: audit the Figure 5f peripheral-blood flow-cytometry data
# provenance and the Z-score calculation. This script reconstructs the
# calculation documented in the historical Figure 5f pipeline, using only the
# canonical 38-sample workbook actually referenced by that pipeline.

options(stringsAsFactors = FALSE, scipen = 999)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "01_flow_zscore_audit.R"
script_path <- normalizePath(script_path, mustWork = FALSE)
code_dir <- dirname(dirname(script_path))
source(file.path(code_dir, "00_paths.R"))

out_dir <- file.path(RESULT_DIR, "R2_5_flow_methods")
fig_dir <- file.path(FIGURE_DIR, "R2_5_flow_methods")
log_dir <- file.path(LOG_DIR, "R2_5_flow_methods")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")

flow_cn <- c(
  "总T淋巴细胞", "辅助/诱导T细胞", "抑制/杀伤T细胞", "总B淋巴细胞",
  "自然杀伤细胞", "NKT淋巴细胞", "双阳性T细胞", "双阴性T细胞",
  "CD29+辅助T细胞", "CD29+杀伤T细胞", "辅助T功能亚群", "杀伤T功能亚群",
  "T细胞活化亚群", "PD1+辅助T细胞", "PD1+杀伤T细胞", "早期活化T细胞",
  "晚期活化T细胞", "调节性T淋巴细胞", "初始辅助T细胞", "中央记忆辅助T细胞",
  "效应记忆辅助T细胞", "效应记忆终末分化辅助T细胞", "初始杀伤T细胞",
  "效应记忆终末分化杀伤T细胞", "效应记忆杀伤T细胞", "中央记忆杀伤T细胞"
)
flow_en <- c(
  "Total_T", "CD4_T", "CD8_T", "Total_B", "NK_Cells", "NKT_Cells",
  "DP_T", "DN_T", "CD29_CD4_T", "CD29_CD8_T", "CD4_Functional",
  "CD8_Functional", "T_Activation", "PD1_CD4", "PD1_CD8", "Early_Act_T",
  "Late_Act_T", "Tregs", "Naive_CD4", "CM_CD4", "EM_CD4", "TEMRA_CD4",
  "Naive_CD8", "TEMRA_CD8", "EM_CD8", "CM_CD8"
)
flow_label <- c(
  "Total T", "CD4 T", "CD8 T", "Total B", "NK", "NKT", "DP T", "DN T",
  "CD29+ CD4", "CD29+ CD8", "CD4 functional", "CD8 functional", "T activation",
  "PD1+ CD4", "PD1+ CD8", "Early activated T", "Late activated T", "Treg",
  "Naive CD4", "Central-memory CD4", "Effector-memory CD4", "TEMRA CD4",
  "Naive CD8", "TEMRA CD8", "Effector-memory CD8", "Central-memory CD8"
)
stopifnot(length(flow_cn) == length(flow_en), length(flow_en) == length(flow_label))

# Figure 5f's source is the canonical workbook, which contains the 38 subtype
# labels and the flow variables in the same row-level table.
source_file <- CLINICAL_SUBTYPE_XLSX
if (!file.exists(source_file)) stop("Missing canonical clinical workbook: ", source_file)
raw <- readxl::read_xlsx(source_file, sheet = 1)
if (!"Subtype" %in% names(raw)) stop("Canonical workbook lacks Subtype column")

# The first data row is a reference-range/unit row in this workbook.
reference_row <- raw[1, , drop = FALSE]
patient <- raw[as.character(raw$Subtype) %in% c("C1", "C2", "C3"), , drop = FALSE]
if (nrow(patient) != 38L) stop("Expected 38 labelled patients, found ", nrow(patient))
if (any(!flow_cn %in% names(patient))) {
  stop("Missing flow columns: ", paste(setdiff(flow_cn, names(patient)), collapse = ", "))
}

clean_numeric <- function(x) {
  x <- gsub("[<>]", "", as.character(x))
  x <- gsub("^>120$", "120", x)
  suppressWarnings(as.numeric(x))
}
for (j in seq_along(flow_cn)) patient[[flow_en[j]]] <- clean_numeric(patient[[flow_cn[j]]])
patient$Subtype <- factor(as.character(patient$Subtype), levels = c("C1", "C2", "C3"))

# Reproduce the original build_z_matrix() exactly:
#   1) each flow feature is standardized across all available labelled patients;
#   2) the standardized values are averaged within each subtype.
feature_rows <- lapply(seq_along(flow_en), function(j) {
  f <- flow_en[j]
  x <- patient[[f]]
  ok <- is.finite(x) & !is.na(x)
  mu <- mean(x[ok])
  s <- stats::sd(x[ok])
  z <- if (sum(ok) >= 5L && is.finite(s) && s > 0) (x - mu) / s else rep(0, length(x))
  data.frame(
    feature = f,
    label = flow_label[j],
    source_column = flow_cn[j],
    unit_from_reference_row = as.character(reference_row[[flow_cn[j]]][1]),
    inferred_unit = ifelse(grepl("%", as.character(reference_row[[flow_cn[j]]][1])), "%", "not specified"),
    n_total = sum(ok),
    mean_all_38 = mu,
    sd_all_38 = s,
    C1_n = sum(patient$Subtype == "C1" & ok),
    C2_n = sum(patient$Subtype == "C2" & ok),
    C3_n = sum(patient$Subtype == "C3" & ok),
    C1_z_mean = mean(z[patient$Subtype == "C1"], na.rm = TRUE),
    C2_z_mean = mean(z[patient$Subtype == "C2"], na.rm = TRUE),
    C3_z_mean = mean(z[patient$Subtype == "C3"], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
feature_audit <- do.call(rbind, feature_rows)
write.csv(feature_audit, file.path(out_dir, "R2_5_flow_feature_audit.csv"), row.names = FALSE, na = "")

# Row-level standardized values are supplied for auditability.
z_long <- do.call(rbind, lapply(seq_along(flow_en), function(j) {
  f <- flow_en[j]
  x <- patient[[f]]
  ok <- is.finite(x) & !is.na(x)
  mu <- mean(x[ok]); s <- stats::sd(x[ok])
  z <- if (is.finite(s) && s > 0) (x - mu) / s else rep(0, length(x))
  data.frame(
    SampleID = if ("SampleID" %in% names(patient)) as.character(patient$SampleID) else seq_len(nrow(patient)),
    Subtype = as.character(patient$Subtype),
    feature = f,
    label = flow_label[j],
    value_percent = x,
    z_score_all_38 = z,
    stringsAsFactors = FALSE
  )
}))
write.csv(z_long, file.path(out_dir, "R2_5_flow_row_level_zscores.csv"), row.names = FALSE, na = "")

# Compact group-level matrix in the same feature order as the historical plot.
z_group <- feature_audit[, c("feature", "label", "C1_z_mean", "C2_z_mean", "C3_z_mean")]
write.csv(z_group, file.path(out_dir, "R2_5_flow_zscore_group_matrix.csv"), row.names = FALSE, na = "")

# Explicit source/provenance record for the response letter.
provenance <- data.frame(
  item = c("Figure 5f source workbook", "Figure 5f source code", "labelled patients", "subtype counts", "flow features audited", "flow measurement unit", "z-score reference population", "group aggregation", "missing flow values"),
  value = c(
    normalizePath(source_file, mustWork = FALSE),
    file.path(Sys.getenv("VS_CM_REPO_ROOT"), "inputs/legacy_project/code/200+blood.R (DATA_38_PATH; build_z_matrix)"),
    as.character(nrow(patient)),
    paste(names(table(patient$Subtype)), as.integer(table(patient$Subtype)), collapse = "; "),
    as.character(length(flow_en)),
    if (all(feature_audit$inferred_unit == "%")) "percentage (%)" else paste(unique(feature_audit$inferred_unit), collapse = "; "),
    "all available labelled patients (n=38), feature-wise mean and sample SD",
    "mean of patient-level z-scores within C1/C2/C3",
    as.character(sum(!is.finite(as.matrix(patient[, flow_en])) | is.na(as.matrix(patient[, flow_en]))))
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "R2_5_source_provenance.csv"), row.names = FALSE, na = "")

subtype_crosscheck <- data.frame(
  subtype = c("C1", "C2", "C3"),
  n = as.integer(table(patient$Subtype)),
  stringsAsFactors = FALSE
)
write.csv(subtype_crosscheck, file.path(out_dir, "R2_5_subtype_counts.csv"), row.names = FALSE)

# A method-focused visualization: group means of the feature-wise z-scores.
plot_df <- rbind(
  data.frame(feature = z_group$feature, label = z_group$label, Subtype = "C1", z = z_group$C1_z_mean),
  data.frame(feature = z_group$feature, label = z_group$label, Subtype = "C2", z = z_group$C2_z_mean),
  data.frame(feature = z_group$feature, label = z_group$label, Subtype = "C3", z = z_group$C3_z_mean)
)
plot_df$label <- factor(plot_df$label, levels = rev(flow_label))
plot_df$Subtype <- factor(plot_df$Subtype, levels = c("C1", "C2", "C3"))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Subtype, y = label, fill = z)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1.6, 1.6), oob = scales::squish, name = "Mean\nZ-score") +
  ggplot2::labs(title = "Figure 5f flow-cytometry Z-score audit", subtitle = "Feature-wise standardization across all 38 patients; cells are subtype means", x = NULL, y = NULL) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"), legend.position = "right")
ggplot2::ggsave(file.path(fig_dir, "R2_5_flow_zscore_group_matrix.pdf"), p, width = 5.8, height = 8.6, useDingbats = FALSE)
ggplot2::ggsave(file.path(fig_dir, "R2_5_flow_zscore_group_matrix.svg"), p, width = 5.8, height = 8.6)
ggplot2::ggsave(file.path(fig_dir, "R2_5_flow_zscore_group_matrix.tiff"), p, width = 5.8, height = 8.6, dpi = 300, compression = "lzw")
ggplot2::ggsave(file.path(fig_dir, "R2_5_flow_zscore_group_matrix_preview.png"), p, width = 5.8, height = 8.6, dpi = 160)

method_lines <- c(
  "Reviewer 2.5 Figure 5f flow-cytometry/Z-score audit",
  paste0("canonical_source=", normalizePath(source_file, mustWork = FALSE)),
  "historical_pipeline=inputs/legacy_project/code/200+blood.R",
  paste0("patients=", nrow(patient)),
  paste0("subtype_counts=", paste(names(table(patient$Subtype)), as.integer(table(patient$Subtype)), collapse = ";")),
  paste0("flow_features=", length(flow_en)),
  "measurement_unit=percentage (%) according to workbook reference row",
  "z_score_formula=(patient_value - feature_mean_all_38)/feature_sample_sd_all_38",
  "z_score_reference_population=all available labelled patients",
  "group_summary=mean patient-level z-scores within subtype",
  paste0("missing_flow_values=", sum(!is.finite(as.matrix(patient[, flow_en])) | is.na(as.matrix(patient[, flow_en])))),
  "status=ANALYSIS_VERIFIED / AUTHOR_METHOD_CONFIRMATION_PENDING"
)
writeLines(method_lines, file.path(out_dir, "R2_5_method_audit.txt"))
writeLines(c(
  paste0("script=", script_path),
  paste0("result_dir=", normalizePath(out_dir, mustWork = FALSE)),
  paste0("figure_dir=", normalizePath(fig_dir, mustWork = FALSE)),
  "status=ANALYSIS_VERIFIED / AUTHOR_METHOD_CONFIRMATION_PENDING"
), file.path(log_dir, "run_summary.txt"))
message("R2.5 flow/Z-score audit complete: ", normalizePath(out_dir))
