# Purpose: CIBERSORT-pass comparisons, cross-method correlations and adjusted clinical analyses
# Manuscript/panels: Supplementary Figure 3a; Supplementary Tables 3,4
# Reviewer comments: R1 Comments 5.3–5.7
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/07_Deconvolution_and_clinical_sensitivity.R"
#!/usr/bin/env Rscript

# Reviewer 1: complete audit using the frozen Cancer Medicine inputs.
# This script creates analysis assets only; it does not edit the manuscript or
# response-letter drafts.

suppressPackageStartupMessages({
  library(readxl)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(FSA)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript so that --file is available.")
script_file <- sub("^--file=", "", script_arg[1])
code_dir <- normalizePath(dirname(script_file), mustWork = TRUE)
work_dir <- normalizePath(file.path(code_dir, "..", ".."), mustWork = TRUE)
result_dir <- file.path(work_dir, "03_Results", "R1_reviewer1_audit")
figure_dir <- file.path(work_dir, "04_Figures", "R1_reviewer1_audit")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

clinical_file <- file.path(work_dir, "02_Code", "R2_6_variable_selection", "input_cl_bulk_with_Subtype.xlsx")
surgery_file <- file.path(work_dir, "02_Code", "R2_6_variable_selection", "input_VS_surgery_analysis.xlsx")
cibersort_file <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "inputs/legacy_project/Figures/step8_Immune/CIBERSORT_Result.csv")
mcp_file <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "inputs/prior_project/VS_revision/output/mcp_counter_results.rds")
score_file <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "inputs/prior_project/supplement/all/step8_Immune/VS_score_summary.csv")
estimate_file <- file.path(work_dir, "..", "02_Imported_Reference_Assets", "03_Existing_Results_Reference", "verified_outputs", "estimate_official", "estimate_official_scores.csv")
model_dir <- file.path(work_dir, "..", "02_Imported_Reference_Assets", "03_Existing_Results_Reference", "verified_outputs", "clinical_models")
external_dir <- file.path(work_dir, "03_Results", "R2_3_external_validation")

stopifnot(file.exists(clinical_file), file.exists(surgery_file), file.exists(cibersort_file),
          file.exists(mcp_file), file.exists(score_file), file.exists(estimate_file))

to_num <- function(x) {
  x <- trimws(as.character(x))
  x[grepl("^$|^NA$|^--$|^—$", x, ignore.case = TRUE)] <- NA_character_
  x <- gsub("[<>]", "", x)
  suppressWarnings(as.numeric(x))
}

write_csv_safe <- function(x, name) {
  write.csv(x, file.path(result_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

# -------------------------------------------------------------------------
# 1. Frozen clinical table and recurrence audit
# -------------------------------------------------------------------------
raw <- read_excel(clinical_file)
dat <- raw %>% filter(Subtype %in% c("C1", "C2", "C3"))
stopifnot(nrow(dat) == 38, !anyDuplicated(dat$SampleID))
surg <- read_excel(surgery_file)
stopifnot(!anyDuplicated(surg$SampleID), setequal(dat$SampleID, surg$SampleID))

surg_note_cols <- names(surg)[vapply(surg, function(x) is.character(x) || is.factor(x), logical(1))]
surg_note_cols <- setdiff(surg_note_cols, "SampleID")
surg_note <- surg %>%
  mutate(Surgical_note_text = apply(as.data.frame(across(all_of(surg_note_cols))), 1,
                                    function(z) paste(z[!is.na(z)], collapse = " | "))) %>%
  select(SampleID, Surgical_note_text)

rec <- dat %>%
  select(SampleID, Subtype, Name, Recurrence, Diagnosis, SurgeryDate, Age, TumorSize_raw) %>%
  left_join(surg_note, by = "SampleID") %>%
  mutate(
    Recurrence_norm = case_when(
      grepl("^yes$|^是$|^1$", as.character(Recurrence), ignore.case = TRUE) ~ "Yes",
      grepl("^no$|^否$|^0$", as.character(Recurrence), ignore.case = TRUE) ~ "No",
      TRUE ~ NA_character_
    ),
    Narrative_explicit_recurrence_flag = grepl("复发", Surgical_note_text),
    Narrative_prior_treatment_flag = grepl("复发|再次|二次手术|前次|上次|曾行|伽玛刀术后|[0-9]+年手术后|[0-9]+月前做过一次手术", Surgical_note_text),
    Max_Dimension_mm_raw = vapply(strsplit(gsub("×", "*", as.character(TumorSize_raw)), "\\*"), function(z) {
      v <- suppressWarnings(as.numeric(z)); if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    }, numeric(1)),
    Max_Dimension_mm_corrected = ifelse(SampleID == "P32", 22, Max_Dimension_mm_raw),
    Recurrence_note_agreement = ifelse(Recurrence_norm == "Yes" & Narrative_explicit_recurrence_flag, "Yes/Yes",
                                       ifelse(Recurrence_norm == "No" & !Narrative_explicit_recurrence_flag, "No/No", "discordant_or_unclear")),
    Prior_treatment_field_note_status = ifelse(Recurrence_norm == "Yes" & Narrative_prior_treatment_flag, "Yes/possible_prior_treatment",
                                               ifelse(Recurrence_norm == "No" & Narrative_prior_treatment_flag, "No/possible_prior_treatment", "no_prior_treatment_flag"))
  )
write_csv_safe(rec, "R1_recurrence_sample_audit.csv")

rec_by_subtype <- rec %>%
  count(Subtype, Recurrence_norm, name = "n") %>%
  group_by(Subtype) %>% mutate(total = sum(n), percent = 100 * n / total) %>% ungroup()
write_csv_safe(rec_by_subtype, "R1_recurrence_by_subtype.csv")
rec_tab <- table(rec$Subtype, rec$Recurrence_norm, useNA = "ifany")
rec_fisher <- fisher.test(rec_tab)
write_csv_safe(data.frame(test = "Subtype_vs_recurrence_Fisher_exact", n = nrow(rec),
                          recurrence_positive = sum(rec$Recurrence_norm == "Yes", na.rm = TRUE),
                          recurrence_positive_percent = 100 * mean(rec$Recurrence_norm == "Yes", na.rm = TRUE),
                          fisher_p = rec_fisher$p.value), "R1_recurrence_test_summary.csv")

# The record-level audit is deliberately descriptive because only three cases
# are recurrence-positive; recurrence is not used as a prognostic endpoint.

# -------------------------------------------------------------------------
# 2. CIBERSORT quality, retention sensitivity, and MCP-counter concordance
# -------------------------------------------------------------------------
cib <- read.csv(cibersort_file, check.names = FALSE, stringsAsFactors = FALSE)
names(cib)[1] <- "SampleID_full"
cib <- cib %>% mutate(SampleID = sub("^tpm\\.", "", SampleID_full)) %>%
  left_join(dat %>% select(SampleID, Subtype), by = "SampleID") %>%
  mutate(Pass_Permutation_P_lt_0_05 = `P-value` < 0.05)
stopifnot(nrow(cib) == 38, all(!is.na(cib$Subtype)))

cib_quality <- cib %>% select(SampleID, Subtype, `P-value`, Correlation, RMSE,
                              Pass_Permutation_P_lt_0_05)
write_csv_safe(cib_quality, "R1_CIBERSORT_sample_quality.csv")
write_csv_safe(data.frame(
  Sample_set = c("All discovery samples", "Permutation-P<0.05 subset"),
  n = c(nrow(cib), sum(cib$Pass_Permutation_P_lt_0_05)),
  percentage_of_38 = c(100, 100 * mean(cib$Pass_Permutation_P_lt_0_05))
), "R1_CIBERSORT_quality_summary.csv")

cib_by_subtype <- cib %>% group_by(Subtype) %>% summarise(
  n_all = n(), n_pass = sum(Pass_Permutation_P_lt_0_05), pass_percent = 100 * n_pass / n_all,
  median_p = median(`P-value`), median_correlation = median(Correlation), median_RMSE = median(RMSE),
  .groups = "drop"
)
write_csv_safe(cib_by_subtype, "R1_CIBERSORT_quality_by_subtype.csv")
pass_subtype_fisher <- fisher.test(table(cib$Subtype, cib$Pass_Permutation_P_lt_0_05))
write_csv_safe(data.frame(test = "Pass_status_by_subtype_Fisher_exact", p = pass_subtype_fisher$p.value),
               "R1_CIBERSORT_pass_status_test.csv")

cov_for_pass <- rec %>% select(SampleID, Age, Max_Dimension_mm_corrected, Recurrence_norm) %>%
  left_join(cib_quality %>% select(SampleID, Pass_Permutation_P_lt_0_05), by = "SampleID")
cov_tests <- lapply(c("Age", "Max_Dimension_mm_corrected"), function(v) {
  x <- to_num(cov_for_pass[[v]]); g <- cov_for_pass$Pass_Permutation_P_lt_0_05
  wt <- wilcox.test(x ~ g, exact = FALSE)
  data.frame(variable = v, test = "Wilcoxon_pass_vs_fail", n = sum(is.finite(x) & !is.na(g)), p = wt$p.value)
}) %>% bind_rows()
cov_tests <- bind_rows(cov_tests, data.frame(variable = "Recurrence_norm", test = "Fisher_pass_vs_fail",
                                             n = nrow(cov_for_pass),
                                             p = fisher.test(table(cov_for_pass$Pass_Permutation_P_lt_0_05,
                                                                  cov_for_pass$Recurrence_norm))$p.value))
write_csv_safe(cov_tests, "R1_CIBERSORT_pass_fail_covariate_tests.csv")

lm22_cols <- setdiff(names(cib), c("SampleID_full", "SampleID", "Subtype", "P-value",
                                   "Correlation", "RMSE", "Pass_Permutation_P_lt_0_05"))
run_kw_dunn <- function(df, value_col, label, sample_set) {
  x <- df[[value_col]]
  g <- droplevels(factor(df$Subtype, levels = c("C1", "C2", "C3")))
  ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 3 || nlevels(droplevels(g[ok])) < 2) return(NULL)
  kt <- kruskal.test(x[ok] ~ g[ok])
  dd <- data.frame(Value = x[ok], Subtype = g[ok])
  dr <- tryCatch(FSA::dunnTest(Value ~ Subtype, data = dd, method = "none")$res,
                 error = function(e) NULL)
  med <- tapply(x[ok], g[ok], median, na.rm = TRUE)
  out <- data.frame(CellType = label, Sample_set = sample_set, n = nrow(dd), KW_p = kt$p.value,
                    Median_C1 = unname(med["C1"]), Median_C2 = unname(med["C2"]),
                    Median_C3 = unname(med["C3"]), stringsAsFactors = FALSE)
  if (!is.null(dr)) {
    for (i in seq_len(nrow(dr))) {
      nm <- gsub(" ", "_", dr$Comparison[i])
      out[[paste0("P_raw_", nm)]] <- dr$P.unadj[i]
    }
  }
  out
}

for (set_name in c("All38", "Filtered22")) {
  dd <- if (set_name == "All38") cib else cib %>% filter(Pass_Permutation_P_lt_0_05)
  rows <- lapply(lm22_cols, function(v) run_kw_dunn(dd, v, v, set_name)) %>% bind_rows()
  if (nrow(rows)) {
    rows$KW_BH_22 <- p.adjust(rows$KW_p, method = "BH")
    pair_cols <- grep("^P_raw_", names(rows), value = TRUE)
    if (length(pair_cols)) {
      pvec <- unlist(rows[pair_cols], use.names = FALSE)
      qvec <- p.adjust(pvec, method = "BH")
      for (j in seq_along(pair_cols)) {
        rows[[sub("^P_raw_", "Q_global66_", pair_cols[j])]] <-
          qvec[((j - 1) * nrow(rows) + 1):(j * nrow(rows))]
      }
    }
  }
  write_csv_safe(rows, paste0("R1_CIBERSORT_LM22_", set_name, ".csv"))
}

# MCP-counter subtype screen, including the Fibroblasts endpoint mentioned by
# Reviewer 1. It is kept separate from the CIBERSORT LM22 screen.
mcp <- readRDS(mcp_file)
mcp_df <- data.frame(SampleID_full = colnames(mcp), t(mcp), check.names = FALSE,
                     stringsAsFactors = FALSE) %>%
  mutate(SampleID = sub("^tpm\\.", "", SampleID_full)) %>%
  left_join(dat %>% select(SampleID, Subtype), by = "SampleID")
mcp_rows <- lapply(rownames(mcp), function(cell) {
  x <- mcp_df[[cell]]; ok <- is.finite(x) & !is.na(mcp_df$Subtype)
  dd <- data.frame(Value = x[ok], Subtype = droplevels(factor(mcp_df$Subtype[ok], levels = c("C1", "C2", "C3"))))
  kt <- kruskal.test(Value ~ Subtype, data = dd)
  dr <- FSA::dunnTest(Value ~ Subtype, data = dd, method = "none")$res
  med <- tapply(dd$Value, dd$Subtype, median)
  data.frame(CellType = cell, n = sum(ok), KW_p = kt$p.value,
             P_C1_C2 = dr$P.unadj[1], P_C1_C3 = dr$P.unadj[2], P_C2_C3 = dr$P.unadj[3],
             Median_C1 = unname(med["C1"]), Median_C2 = unname(med["C2"]),
             Median_C3 = unname(med["C3"]), stringsAsFactors = FALSE)
}) %>% bind_rows()
mcp_rows$KW_BH_10 <- p.adjust(mcp_rows$KW_p, method = "BH")
all_mcp_p <- unlist(mcp_rows[c("P_C1_C2", "P_C1_C3", "P_C2_C3")], use.names = FALSE)
all_mcp_q <- p.adjust(all_mcp_p, method = "BH")
mcp_rows$Q_C1_C2_global30 <- all_mcp_q[seq_len(nrow(mcp_rows))]
mcp_rows$Q_C1_C3_global30 <- all_mcp_q[nrow(mcp_rows) + seq_len(nrow(mcp_rows))]
mcp_rows$Q_C2_C3_global30 <- all_mcp_q[2 * nrow(mcp_rows) + seq_len(nrow(mcp_rows))]
write_csv_safe(mcp_rows, "R1_MCP_counter_celltype_subtype_screen.csv")

# Harmonized CIBERSORT–MCP-counter correlations. The manuscript discrepancy is
# retained as an audit: Pearson, Spearman asymptotic, and requested exact
# Spearman calculations are all exported from the same 38 samples.
cib_cor <- cib %>% transmute(SampleID_full,
                             CIBERSORT_CD8 = `T cells CD8`,
                             # These mappings reproduce the historical
                             # cross-method comparison exactly; they are not
                             # redefined as aggregate biological compartments.
                             CIBERSORT_B = `B cells naive`,
                             CIBERSORT_Monocytic = `Macrophages M2`,
                             CIBERSORT_NK = `NK cells resting`)
mcp_cor <- data.frame(SampleID_full = colnames(mcp), t(mcp), check.names = FALSE,
                      stringsAsFactors = FALSE) %>%
  transmute(SampleID_full, MCP_CD8 = `CD8 T cells`, MCP_B = `B lineage`,
            MCP_Monocytic = `Monocytic lineage`, MCP_NK = `NK cells`)
cor_map <- data.frame(
  Population = c("CD8 T cells", "B lineage", "Monocytic lineage", "NK cells"),
  MCP_counter_component = c("CD8 T cells", "B lineage", "Monocytic lineage", "NK cells"),
  CIBERSORT_component = c("T cells CD8", "B cells naive", "Macrophages M2", "NK cells resting"),
  Mapping_status = c("direct nominal match", rep("historical component proxy", 3)),
  mcp = c("MCP_CD8", "MCP_B", "MCP_Monocytic", "MCP_NK"),
  cib = c("CIBERSORT_CD8", "CIBERSORT_B", "CIBERSORT_Monocytic", "CIBERSORT_NK"),
  stringsAsFactors = FALSE
)
cor_rows <- lapply(seq_len(nrow(cor_map)), function(i) {
  z <- inner_join(mcp_cor, cib_cor, by = "SampleID_full")
  x <- z[[cor_map$cib[i]]]; y <- z[[cor_map$mcp[i]]]
  ok <- complete.cases(x, y)
  has_ties <- anyDuplicated(x[ok]) > 0L || anyDuplicated(y[ok]) > 0L
  s_asym <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  s_exact <- suppressWarnings(cor.test(x, y, method = "spearman", exact = TRUE))
  p <- cor.test(x, y, method = "pearson")
  data.frame(
    Cell_Population = cor_map$Population[i],
    MCP_counter_component = cor_map$MCP_counter_component[i],
    CIBERSORT_component = cor_map$CIBERSORT_component[i],
    Mapping_status = cor_map$Mapping_status[i],
    N = sum(complete.cases(x, y)),
    Spearman_Rho = unname(s_asym$estimate),
    P_Spearman_asymptotic = s_asym$p.value,
    P_Spearman_exact_TRUE_request = s_exact$p.value,
    Pearson_R = unname(p$estimate), P_Pearson = p$p.value,
    Mean_MCP = mean(y, na.rm = TRUE), Mean_CIBERSORT = mean(x, na.rm = TRUE),
    Spearman_exact_TRUE_method = if (has_ties) "asymptotic fallback because of ties" else if (sum(ok) < 10L) "exact" else "AS89/Edgeworth approximation"
  )
}) %>% bind_rows()
write_csv_safe(cor_rows, "R1_CIBERSORT_MCP_correlation_reanalysis.csv")
write_csv_safe(cor_rows %>% filter(Cell_Population == "CD8 T cells"), "R1_CD8_pvalue_conventions.csv")

# -------------------------------------------------------------------------
# 3. Official ESTIMATE and custom fibrosis/ECM score audit
# -------------------------------------------------------------------------
est <- read_csv(estimate_file, show_col_types = FALSE) %>%
  mutate(Subtype = factor(Subtype, levels = c("C1", "C2", "C3")))
est_metrics <- c("Stromal_Score", "Immune_Score", "ESTIMATE_Score", "Tumour_Purity")
est_summary <- lapply(est_metrics, function(v) {
  x <- est[[v]]; kt <- kruskal.test(x ~ est$Subtype)
  med <- tapply(x, est$Subtype, median); means <- tapply(x, est$Subtype, mean)
  data.frame(Metric = v, n = sum(is.finite(x)), KW_p = kt$p.value,
             Mean_C1 = means["C1"], Mean_C2 = means["C2"], Mean_C3 = means["C3"],
             Median_C1 = med["C1"], Median_C2 = med["C2"], Median_C3 = med["C3"])
}) %>% bind_rows() %>% mutate(KW_BH_4 = p.adjust(KW_p, method = "BH"))
write_csv_safe(est_summary, "R1_official_ESTIMATE_subtype_summary.csv")

est_pair <- lapply(est_metrics, function(v) {
  dd <- data.frame(Value = est[[v]], Subtype = est$Subtype)
  dr <- FSA::dunnTest(Value ~ Subtype, data = dd, method = "none")$res
  data.frame(Metric = v, Comparison = dr$Comparison, P_raw = dr$P.unadj,
             P_BH_within_metric = p.adjust(dr$P.unadj, "BH"))
}) %>% bind_rows()
est_pair$Q_global12 <- p.adjust(est_pair$P_raw, "BH")
write_csv_safe(est_pair, "R1_official_ESTIMATE_pairwise_tests.csv")

score <- read_csv(score_file, show_col_types = FALSE) %>%
  mutate(Subtype = factor(Subtype, levels = c("C1", "C2", "C3")))
custom_metrics <- intersect(c("Fibrosis_Score", "ECM_Score", "Stromal_Score", "Immune_Score"), names(score))
custom_summary <- lapply(custom_metrics, function(v) {
  x <- score[[v]]; kt <- kruskal.test(x ~ score$Subtype)
  med <- tapply(x, score$Subtype, median); means <- tapply(x, score$Subtype, mean)
  data.frame(Metric = v, n = sum(is.finite(x)), KW_p = kt$p.value,
             Mean_C1 = means["C1"], Mean_C2 = means["C2"], Mean_C3 = means["C3"],
             Median_C1 = med["C1"], Median_C2 = med["C2"], Median_C3 = med["C3"])
}) %>% bind_rows() %>% mutate(KW_BH = p.adjust(KW_p, "BH"))
write_csv_safe(custom_summary, "R1_custom_score_subtype_summary.csv")

custom_pair <- lapply(custom_metrics, function(v) {
  dd <- data.frame(Value = score[[v]], Subtype = score$Subtype)
  dr <- FSA::dunnTest(Value ~ Subtype, data = dd, method = "none")$res
  data.frame(Metric = v, Comparison = dr$Comparison, P_raw = dr$P.unadj,
             P_BH_within_metric = p.adjust(dr$P.unadj, "BH"))
}) %>% bind_rows()
custom_pair$Q_global <- p.adjust(custom_pair$P_raw, "BH")
write_csv_safe(custom_pair, "R1_custom_score_pairwise_tests.csv")

# -------------------------------------------------------------------------
# 4. Size-adjusted surgical phenotype audit
# -------------------------------------------------------------------------
clinical <- dat %>%
  mutate(
    Age_num = to_num(Age),
    Max_Dimension_mm_raw = vapply(strsplit(gsub("×", "*", as.character(TumorSize_raw)), "\\*"), function(z) {
      v <- suppressWarnings(as.numeric(z)); if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    }, numeric(1)),
    Max_Dimension_mm = ifelse(SampleID == "P32", 22, Max_Dimension_mm_raw),
    Size_cm = Max_Dimension_mm / 10
  ) %>%
  left_join(surg %>% select(SampleID, `肿瘤质地`, `脑干面粘连`), by = "SampleID") %>%
  mutate(
    Texture = case_when(grepl("软|Soft", `肿瘤质地`, ignore.case = TRUE) ~ "Soft",
                        grepl("中|Medium", `肿瘤质地`, ignore.case = TRUE) ~ "Medium",
                        grepl("硬|Hard", `肿瘤质地`, ignore.case = TRUE) ~ "Hard",
                        TRUE ~ NA_character_),
    Texture_ord = ordered(Texture, levels = c("Soft", "Medium", "Hard")),
    Adhesion = case_when(grepl("有|Yes|Y|1|丰富|富", `脑干面粘连`, ignore.case = TRUE) ~ 1,
                         grepl("无|No|N|0", `脑干面粘连`, ignore.case = TRUE) ~ 0,
                         TRUE ~ NA_real_)
  )
write_csv_safe(clinical %>% select(SampleID, Subtype, Age_num, TumorSize_raw,
                                   Max_Dimension_mm_raw, Max_Dimension_mm, Size_cm,
                                   Texture, Adhesion, Recurrence),
               "R1_clinical_model_input_audit.csv")

model_rows <- list()
for (model_name in c("texture_unadjusted", "texture_size_adjusted", "texture_size_age_adjusted")) {
  rhs <- switch(model_name,
                texture_unadjusted = "Subtype",
                texture_size_adjusted = "Subtype + Size_cm",
                texture_size_age_adjusted = "Subtype + Size_cm + Age_num")
  dd <- clinical %>% drop_na(Texture_ord, Subtype, Size_cm, Age_num)
  if (model_name == "texture_unadjusted") dd <- clinical %>% drop_na(Texture_ord, Subtype)
  fit <- MASS::polr(as.formula(paste("Texture_ord ~", rhs)), data = dd,
                    Hess = TRUE, method = "logistic")
  co <- coef(summary(fit)); ci <- suppressWarnings(confint.default(fit))
  keep <- grep("^Subtype", rownames(co), value = TRUE)
  if (length(keep)) model_rows[[model_name]] <- data.frame(
    model = model_name, outcome = "Texture_ordinal", n = nrow(dd), term = keep,
    OR = exp(co[keep, "Value"]), CI_low_wald = exp(ci[keep, 1]),
    CI_high_wald = exp(ci[keep, 2]),
    P = 2 * pnorm(abs(co[keep, "t value"]), lower.tail = FALSE), stringsAsFactors = FALSE)
}

# Ordinary logistic sensitivities are exported for transparency. The verified
# Firth coefficients are copied below as the official small-sample estimates.
for (model_name in c("adhesion_unadjusted", "adhesion_size_adjusted", "adhesion_size_age_adjusted")) {
  rhs <- switch(model_name,
                adhesion_unadjusted = "Subtype",
                adhesion_size_adjusted = "Subtype + Size_cm",
                adhesion_size_age_adjusted = "Subtype + Size_cm + Age_num")
  dd <- clinical %>% drop_na(Adhesion, Subtype, Size_cm, Age_num)
  if (model_name == "adhesion_unadjusted") dd <- clinical %>% drop_na(Adhesion, Subtype)
  fit <- suppressWarnings(glm(as.formula(paste("Adhesion ~", rhs)), data = dd, family = binomial()))
  co <- coef(summary(fit)); keep <- grep("^Subtype", rownames(co), value = TRUE)
  if (length(keep)) model_rows[[model_name]] <- data.frame(
    model = paste0(model_name, "_ordinary_logistic_sensitivity"),
    outcome = "Brainstem_adhesion", n = nrow(dd), term = keep,
    OR = exp(co[keep, "Estimate"]),
    CI_low_wald = exp(co[keep, "Estimate"] - 1.96 * co[keep, "Std. Error"]),
    CI_high_wald = exp(co[keep, "Estimate"] + 1.96 * co[keep, "Std. Error"]),
    P = co[keep, "Pr(>|z|)"], stringsAsFactors = FALSE)
}
write_csv_safe(bind_rows(model_rows), "R1_size_adjusted_models_reanalysis.csv")

if (dir.exists(model_dir)) {
  official_coef <- read.csv(file.path(model_dir, "clinical_model_coefficients.csv"), stringsAsFactors = FALSE)
  official_omni <- read.csv(file.path(model_dir, "clinical_omnibus_tests.csv"), stringsAsFactors = FALSE)
  write_csv_safe(official_coef, "R1_official_verified_clinical_model_coefficients.csv")
  write_csv_safe(official_omni, "R1_official_verified_clinical_omnibus_tests.csv")
}

size_kw <- kruskal.test(Max_Dimension_mm ~ Subtype, data = clinical)
write_csv_safe(data.frame(test = "Tumor_size_by_subtype_KW", n = nrow(clinical), p = size_kw$p.value,
                          median_C1 = median(clinical$Max_Dimension_mm[clinical$Subtype == "C1"]),
                          median_C2 = median(clinical$Max_Dimension_mm[clinical$Subtype == "C2"]),
                          median_C3 = median(clinical$Max_Dimension_mm[clinical$Subtype == "C3"])),
               "R1_tumor_size_by_subtype.csv")

# -------------------------------------------------------------------------
# 5. External validation reconciliation
# -------------------------------------------------------------------------
# R2.3 already produced a frozen, same-signature transfer analysis. For R1 we
# audit those files and copy only the verified summary assets; we do not create
# a second external-data version or silently replace the negative cohort result.
external_files <- c("external_mgs_tests.csv", "external_assignment_summary.csv",
                    "external_signature_coverage.csv", "external_frozen_rank_concordance.csv")
external_audit <- lapply(external_files, function(fn) {
  f <- file.path(external_dir, fn)
  if (!file.exists(f)) return(data.frame(file = fn, exists = FALSE, rows = NA, columns = NA))
  x <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(file = fn, exists = TRUE, rows = nrow(x), columns = ncol(x))
}) %>% bind_rows()
write_csv_safe(external_audit, "R1_external_validation_asset_audit.csv")
for (fn in external_files) {
  f <- file.path(external_dir, fn)
  if (file.exists(f)) file.copy(f, file.path(result_dir, paste0("R1_", fn)), overwrite = TRUE)
}

# -------------------------------------------------------------------------
# 6. Diagnostic figures and audit summary
# -------------------------------------------------------------------------
ggplot(cib_quality, aes(x = Subtype, fill = Pass_Permutation_P_lt_0_05)) +
  geom_bar(position = "dodge", colour = "black") +
  scale_fill_manual(values = c(`TRUE` = "#4DAF4A", `FALSE` = "#BDBDBD"),
                    name = "CIBERSORT\npermutation P<0.05") +
  labs(x = NULL, y = "Number of samples", title = "CIBERSORT quality by discovery group") +
  theme_classic(base_size = 12) -> p_cib
ggsave(file.path(figure_dir, "R1_CIBERSORT_quality_by_subtype.pdf"), p_cib, width = 6.5, height = 4.5)
ggsave(file.path(figure_dir, "R1_CIBERSORT_quality_by_subtype.png"), p_cib, width = 6.5, height = 4.5, dpi = 300)

est_long <- est %>% select(SampleID, Subtype, all_of(est_metrics)) %>%
  pivot_longer(all_of(est_metrics), names_to = "Metric", values_to = "Value")
ggplot(est_long, aes(Subtype, Value, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.12, size = 1.3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) + theme_classic(base_size = 11) +
  theme(legend.position = "none") -> p_est
ggsave(file.path(figure_dir, "R1_official_ESTIMATE_by_subtype.pdf"), p_est, width = 8, height = 6)
ggsave(file.path(figure_dir, "R1_official_ESTIMATE_by_subtype.png"), p_est, width = 8, height = 6, dpi = 300)

summary_lines <- c(
  "Reviewer 1 analysis audit (frozen Cancer Medicine workspace)",
  paste0("Date: ", Sys.Date()),
  "No manuscript or response-letter file was modified by this script.",
  "",
  "R1.m1 recurrence audit:",
  paste0("  n=", nrow(rec), "; recurrence-positive=", sum(rec$Recurrence_norm == "Yes", na.rm = TRUE),
         " (", sprintf("%.1f", 100 * mean(rec$Recurrence_norm == "Yes", na.rm = TRUE)), "%)."),
  paste0("  Subtype-vs-recurrence Fisher P=", formatC(rec_fisher$p.value, format = "e", digits = 3), "."),
  paste0("  Narrative/field discordant-or-unclear cases=", sum(rec$Recurrence_note_agreement == "discordant_or_unclear", na.rm = TRUE), "."),
  "",
  "R1.M4 CIBERSORT reliability:",
  paste0("  permutation P<0.05: ", sum(cib$Pass_Permutation_P_lt_0_05), "/", nrow(cib),
         " (", sprintf("%.1f", 100 * mean(cib$Pass_Permutation_P_lt_0_05)), "%)."),
  paste0("  pass-status-by-subtype Fisher P=", formatC(pass_subtype_fisher$p.value, format = "e", digits = 3), "."),
  "  CIBERSORT LM22 all-sample and pass-only tables written.",
  "  MCP-counter fibroblast and other cell-type screen written.",
  "",
  "R1.M5 CIBERSORT-MCP concordance:",
  "  Pearson and Spearman (asymptotic and AS89/Edgeworth exact=TRUE) values written from the same 38 samples.",
  "  Only CD8 is a direct nominal component match; B, monocytic and NK comparisons are historical component proxies.",
  "",
  "R1.M3 ESTIMATE:",
  "  Official Stromal/Immune/ESTIMATE/Purity means, omnibus and pairwise tests written.",
  "  Custom Fibrosis/ECM score sensitivity table written separately.",
  "",
  "R1.M7 size adjustment:",
  "  Texture ordinal models and ordinary-logistic adhesion sensitivities rerun from canonical clinical inputs.",
  "  Verified Firth coefficients and omnibus tests copied into the R1 audit directory for comparison.",
  "",
  "R1.M2 external validation:",
  "  Frozen R2.3 external outputs audited and copied; no alternative version silently substituted.",
  "",
  "Remaining non-computational confirmations: author/affiliation changes, recurrence field meaning if needed,",
  "flow/laboratory metadata, data deposition, funding, graphical-abstract provenance, and literature additions."
)
writeLines(summary_lines, file.path(result_dir, "R1_complete_audit_summary.txt"), useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")
