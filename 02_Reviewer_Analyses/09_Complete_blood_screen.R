# Purpose: 83 group tests, 249 Dunn comparisons and 168 operative-phenotype correlations
# Manuscript/panels: Supplementary Figure 5c,d; Figure 5g–i interpretation
# Reviewer comments: R1 Comment 5.6; R2 Comment 7
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/09_Complete_blood_screen.R"
#!/usr/bin/env Rscript

# Reviewer 2, Major comment 6
# Complete audit of variable selection for Figure 5g-h/i.
#
# This script intentionally starts from the canonical 38-patient workbook and
# re-runs the full subtype and surgical-phenotype screens. It also reproduces
# the historical Figure 5 candidate screen (seven variables; min P < 0.1)
# without treating that screen as an a priori analysis.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(FSA)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
code_dir <- normalizePath(dirname(script_file), mustWork = TRUE)
work_dir <- normalizePath(file.path(code_dir, "..", ".."), mustWork = TRUE)
result_dir <- file.path(work_dir, "03_Results", "R2_6_variable_selection")
figure_dir <- file.path(work_dir, "04_Figures", "R2_6_variable_selection")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_clinical <- file.path(code_dir, "input_cl_bulk_with_Subtype.xlsx")
input_surgery <- file.path(code_dir, "input_VS_surgery_analysis.xlsx")
stopifnot(file.exists(input_clinical), file.exists(input_surgery))

to_num <- function(x) {
  x <- trimws(as.character(x))
  x[grepl("^$|^NA$|^--$|^—$", x, ignore.case = TRUE)] <- NA_character_
  x <- gsub("[<>]", "", x)
  suppressWarnings(as.numeric(x))
}

# Map the canonical Chinese workbook columns to the names used in the
# historical analysis. The six derived indices are calculated below.
feature_map <- c(
  "白细胞计数" = "WBC", "红细胞计数" = "RBC", "中性粒细胞%" = "Neutrophil_pct",
  "血红蛋白量" = "Hemoglobin", "淋巴细胞%" = "Lymphocyte_pct", "红细胞比积" = "Hematocrit",
  "单核细胞%" = "Monocyte_pct", "平均红细胞体积" = "MCV", "嗜酸性粒细胞%" = "Eosinophil_pct",
  "平均红细胞血红蛋白含量" = "MCH", "嗜碱性粒细胞%" = "Basophil_pct",
  "中性粒细胞绝对数" = "ANC", "红细胞体积分布宽度变异系数" = "RDW_CV",
  "淋巴细胞绝对数" = "ALC", "血小板计数" = "Platelets", "单核细胞绝对数" = "AMC",
  "平均血小板体积" = "MPV", "嗜酸性粒细胞绝对数" = "AEC", "嗜碱性粒细胞绝对数" = "ABC",
  "血小板比积" = "Plateletcrit", "血小板体积分布宽度变异系数" = "PDW",
  "平均红细胞血红蛋白浓度" = "MCHC", "CRP" = "CRP",
  "肌酐" = "Creatinine", "白蛋白" = "Albumin", "碱性磷酸酶" = "ALP",
  "丙氨酸氨基转移酶" = "ALT", "天冬氨酸氨基转移酶" = "AST", "直接胆红素" = "D_Bilirubin",
  "γ-谷氨酰转肽酶" = "GGT", "葡萄糖" = "Glucose", "钾" = "Potassium",
  "乳酸脱氢酶" = "LDH", "钠" = "Sodium", "前白蛋白" = "Prealbumin",
  "视黄醇结合蛋白" = "RBP", "总胆汁酸" = "TBA", "总胆红素" = "T_Bilirubin",
  "甘油三酯" = "Triglycerides", "总蛋白" = "Total_Protein", "尿酸" = "Uric_Acid",
  "尿素" = "Urea", "线粒体-AST" = "mAST", "谷氨酸脱氢酶" = "GLDH",
  "白/球比" = "AG_Ratio", "球蛋白" = "Globulin", "氯" = "Chloride",
  "总胆固醇" = "Total_Cholesterol", "AST/ALT比值" = "AST_ALT_Ratio",
  "间接胆红素" = "I_Bilirubin", "估算肾小球滤过率(CKD-EPI)" = "eGFR",
  "总T淋巴细胞" = "Total_T", "辅助/诱导T细胞" = "CD4_T", "抑制/杀伤T细胞" = "CD8_T",
  "总B淋巴细胞" = "Total_B", "自然杀伤细胞" = "NK_Cells", "NKT淋巴细胞" = "NKT_Cells",
  "双阳性T细胞" = "DP_T", "双阴性T细胞" = "DN_T", "CD29+辅助T细胞" = "CD29_CD4_T",
  "CD29+杀伤T细胞" = "CD29_CD8_T", "辅助T功能亚群" = "CD4_Functional",
  "杀伤T功能亚群" = "CD8_Functional", "T细胞活化亚群" = "T_Activation",
  "PD1+辅助T细胞" = "PD1_CD4", "PD1+杀伤T细胞" = "PD1_CD8",
  "早期活化T细胞" = "Early_Act_T", "晚期活化T细胞" = "Late_Act_T",
  "调节性T淋巴细胞" = "Tregs", "初始辅助T细胞" = "Naive_CD4",
  "中央记忆辅助T细胞" = "CM_CD4", "效应记忆辅助T细胞" = "EM_CD4",
  "效应记忆终末分化辅助T细胞" = "TEMRA_CD4", "初始杀伤T细胞" = "Naive_CD8",
  "效应记忆终末分化杀伤T细胞" = "TEMRA_CD8", "效应记忆杀伤T细胞" = "EM_CD8",
  "中央记忆杀伤T细胞" = "CM_CD8"
)

raw <- read_excel(input_clinical) %>% filter(Subtype %in% c("C1", "C2", "C3"))
stopifnot(nrow(raw) == 38, !anyDuplicated(raw$SampleID))

surg <- read_excel(input_surgery) %>% select(SampleID, all_of(c("肿瘤质地", "脑干面粘连")))
stopifnot(!anyDuplicated(surg$SampleID), setequal(raw$SampleID, surg$SampleID))

dat <- raw %>%
  left_join(surg, by = "SampleID") %>%
  mutate(
    Subtype = factor(Subtype, levels = c("C1", "C2", "C3")),
    Age = to_num(Age),
    Texture = case_when(
      grepl("软|Soft", .data[["肿瘤质地"]], ignore.case = TRUE) ~ "Soft",
      grepl("中|Medium", .data[["肿瘤质地"]], ignore.case = TRUE) ~ "Medium",
      grepl("硬|Hard", .data[["肿瘤质地"]], ignore.case = TRUE) ~ "Hard",
      TRUE ~ NA_character_
    ),
    Texture_Score = match(Texture, c("Soft", "Medium", "Hard")),
    Adhesion = case_when(
      grepl("有|Yes|Y|1|丰富|富", .data[["脑干面粘连"]], ignore.case = TRUE) ~ "Yes",
      grepl("无|No|N|0", .data[["脑干面粘连"]], ignore.case = TRUE) ~ "No",
      TRUE ~ NA_character_
    ),
    Adhesion_Score = ifelse(Adhesion == "Yes", 1, ifelse(Adhesion == "No", 0, NA_real_))
  )

for (cn in names(feature_map)) {
  if (cn %in% names(dat)) dat[[feature_map[[cn]]]] <- to_num(dat[[cn]])
}

dat <- dat %>% mutate(
  NLR = ANC / ALC,
  PLR = Platelets / ALC,
  MLR = AMC / ALC,
  SII = Platelets * ANC / ALC,
  CD4_CD8_Ratio = CD4_T / CD8_T,
  SIRI = ANC * AMC / ALC
)

feature_order <- c(
  "CRP", "Total_T", "CD4_T", "CD8_T", "Total_B", "NK_Cells", "NKT_Cells", "DP_T", "DN_T",
  "CD29_CD4_T", "CD29_CD8_T", "CD4_Functional", "CD8_Functional", "T_Activation", "PD1_CD4", "PD1_CD8",
  "Early_Act_T", "Late_Act_T", "Tregs", "Naive_CD4", "CM_CD4", "EM_CD4", "TEMRA_CD4", "Naive_CD8",
  "TEMRA_CD8", "EM_CD8", "CM_CD8", "WBC", "RBC", "Neutrophil_pct", "Hemoglobin", "Lymphocyte_pct",
  "Hematocrit", "Monocyte_pct", "MCV", "Eosinophil_pct", "MCH", "Basophil_pct", "ANC", "RDW_CV",
  "ALC", "Platelets", "AMC", "MPV", "AEC", "ABC", "Plateletcrit", "PDW", "MCHC",
  "Creatinine", "Albumin", "ALP", "ALT", "AST", "D_Bilirubin", "GGT", "Glucose", "Potassium", "LDH",
  "Sodium", "Prealbumin", "RBP", "TBA", "T_Bilirubin", "Triglycerides", "Total_Protein", "Uric_Acid",
  "Urea", "mAST", "GLDH", "AG_Ratio", "Globulin", "Chloride", "Total_Cholesterol", "AST_ALT_Ratio",
  "I_Bilirubin", "eGFR", "NLR", "PLR", "MLR", "SII", "CD4_CD8_Ratio", "SIRI"
)
stopifnot(length(feature_order) == 83, all(feature_order %in% names(dat)))

category <- function(x) {
  case_when(
    x %in% c("Total_B", "NK_Cells", "NKT_Cells") ~ "NK/B cell",
    x %in% c("Total_T", "CD4_T", "CD8_T", "DP_T", "DN_T", "Tregs", "CD4_CD8_Ratio") ~ "T-cell core",
    x %in% c("CD29_CD4_T", "CD29_CD8_T", "CD4_Functional", "CD8_Functional", "T_Activation", "PD1_CD4", "PD1_CD8", "Early_Act_T", "Late_Act_T") ~ "T-cell function",
    x %in% c("Naive_CD4", "CM_CD4", "EM_CD4", "TEMRA_CD4", "Naive_CD8", "TEMRA_CD8", "EM_CD8", "CM_CD8") ~ "T-cell memory",
    x %in% c("CRP", "WBC", "RBC", "Neutrophil_pct", "Hemoglobin", "Lymphocyte_pct", "Hematocrit", "Monocyte_pct", "MCV", "Eosinophil_pct", "MCH", "Basophil_pct", "ANC", "RDW_CV", "ALC", "Platelets", "AMC", "MPV", "AEC", "ABC", "Plateletcrit", "PDW", "MCHC", "NLR", "PLR", "MLR", "SII", "SIRI") ~ "CBC/inflammation",
    TRUE ~ "Biochemistry"
  )
}

safe_kw <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if (length(x) < 3 || nlevels(g) < 2) return(c(statistic = NA_real_, p = NA_real_))
  kt <- tryCatch(kruskal.test(x ~ g), error = function(e) NULL)
  if (is.null(kt)) c(statistic = NA_real_, p = NA_real_) else c(statistic = unname(kt$statistic), p = kt$p.value)
}

# 1. Complete subtype screen: all 83 variables, with BH across all omnibus tests.
kw_rows <- lapply(feature_order, function(f) {
  x <- dat[[f]]
  z <- safe_kw(x, dat$Subtype)
  ok <- is.finite(x)
  med <- tapply(x[ok], dat$Subtype[ok], median, na.rm = TRUE)
  data.frame(
    Feature = f, Category = category(f), N_nonmiss = sum(ok),
    KW_chi2 = unname(z[["statistic"]]), KW_pvalue = unname(z[["p"]]),
    Epsilon2 = ifelse(is.na(z[["statistic"]]), NA_real_, (z[["statistic"]] - 3 + 1) / (sum(ok) - 3)),
    Median_C1 = unname(med["C1"]), Median_C2 = unname(med["C2"]), Median_C3 = unname(med["C3"]),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>% mutate(FDR_BH_83 = p.adjust(KW_pvalue, method = "BH")) %>% arrange(KW_pvalue)
write.csv(kw_rows, file.path(result_dir, "R2_6_complete_subtype_screen_83_features.csv"), row.names = FALSE, na = "")

# 2. Complete pairwise Dunn screen. P.adj_within_feature is the ordinary
# three-comparison BH correction; q_global_249 additionally treats all 249
# pairwise tests as one family for transparent selective-reporting control.
dunn_rows <- lapply(feature_order, function(f) {
  x <- dat[[f]]
  ok <- is.finite(x) & !is.na(dat$Subtype)
  dd <- data.frame(Value = x[ok], Subtype = droplevels(dat$Subtype[ok]))
  if (nrow(dd) < 3 || nlevels(dd$Subtype) < 2) return(NULL)
  out <- tryCatch({
    invisible(capture.output(tmp_dunn <- FSA::dunnTest(Value ~ Subtype, data = dd, method = "none")))
    tmp_dunn$res
  }, error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out %>% transmute(
    Feature = f, Comparison = Comparison, Z = Z, P_raw = P.unadj,
    P_adj_within_feature_BH = p.adjust(P.unadj, method = "BH")
  )
}) %>% bind_rows() %>% mutate(Q_global_249_BH = p.adjust(P_raw, method = "BH"))
write.csv(dunn_rows, file.path(result_dir, "R2_6_complete_pairwise_dunn_249_tests.csv"), row.names = FALSE, na = "")

# 3. Complete surgical-phenotype correlation screen: age plus 83 blood
# variables x two
# prespecified endpoints. Texture is the ordinal Soft < Medium < Hard score;
# adhesion is No=0/Yes=1. Both endpoint-wise and all-166-test BH values are
# exported; the latter is the conservative family for selective reporting.
corr_features <- c("Age", feature_order)
corr_rows <- lapply(corr_features, function(f) {
  x <- dat[[f]]
  one <- function(y, endpoint) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
      return(data.frame(Feature = f, Endpoint = endpoint, N = sum(ok), Spearman_rho = NA_real_, P_raw = NA_real_))
    }
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.frame(Feature = f, Endpoint = endpoint, N = sum(ok), Spearman_rho = unname(ct$estimate), P_raw = ct$p.value)
  }
  bind_rows(one(dat$Texture_Score, "Tumor_texture_ordinal"), one(dat$Adhesion_Score, "Brainstem_adhesion"))
}) %>% bind_rows() %>%
  group_by(Endpoint) %>% mutate(Q_endpoint_BH = p.adjust(P_raw, method = "BH")) %>% ungroup() %>%
  mutate(Q_global_168_BH = p.adjust(P_raw, method = "BH")) %>%
  arrange(Endpoint, P_raw)
write.csv(corr_rows, file.path(result_dir, "R2_6_complete_clinical_correlation_screen_168_tests.csv"), row.names = FALSE, na = "")

# 4. Exact historical Figure 5 candidate screen reproduced from the original
# code: seven variables, complete cases for each variable, and min(P_texture,
# P_adhesion) < 0.1. This is an audit of the old rule, not the recommended
# inferential analysis.
historical_candidates <- c("NK_Cells", "Albumin", "Age", "MPV", "Total_T", "eGFR", "AG_Ratio")
hist_rows <- lapply(historical_candidates, function(f) {
  x <- dat[[f]]
  ok <- is.finite(x) & is.finite(dat$Texture_Score) & is.finite(dat$Adhesion_Score)
  pt <- suppressWarnings(cor.test(dat$Texture_Score[ok], x[ok], method = "spearman", exact = FALSE))
  pa <- suppressWarnings(cor.test(dat$Adhesion_Score[ok], x[ok], method = "spearman", exact = FALSE))
  data.frame(
    Candidate = f, N_complete = sum(ok),
    rho_texture = unname(pt$estimate), P_texture = pt$p.value,
    rho_adhesion = unname(pa$estimate), P_adhesion = pa$p.value,
    Min_P_historical_rule = min(pt$p.value, pa$p.value),
    Selected_by_historical_MinP_lt_0_1 = min(pt$p.value, pa$p.value) < 0.1,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
all_candidate_p <- c(hist_rows$P_texture, hist_rows$P_adhesion)
all_candidate_q <- p.adjust(all_candidate_p, method = "BH")
hist_rows$Q_texture_over_14_candidate_tests <- all_candidate_q[seq_len(nrow(hist_rows))]
hist_rows$Q_adhesion_over_14_candidate_tests <- all_candidate_q[nrow(hist_rows) + seq_len(nrow(hist_rows))]
write.csv(hist_rows, file.path(result_dir, "R2_6_historical_figure5_candidate_screen.csv"), row.names = FALSE, na = "")

# 5. Record which variables are visibly retained in the submitted Figure 5
# panels and where they came from in the historical code.
panel_provenance <- data.frame(
  Panel = c("Figure5g_correlogram", "Figure5h_subgroup_lab", "Figure5i_subgroup_lab"),
  Variables_visible = c("Tumor Texture; Adhesion; NK Cells; Age; MPV; Total T Cells", "NK Cells by tumour size and subtype", "NK Cells and Age by adhesion; Total T Cells by texture"),
  Historical_rule_or_source = c("candidate blood_dict + Min_P < 0.1 screen in figure5/术中 code; eGFR passed the reproduced rule but is absent from the displayed panel", "manually specified subgroup panel in figure5 code; not an a priori analysis plan", "manually specified subgroup panel in figure5 code; not an a priori analysis plan"),
  Reviewer2_6_issue = c("selective subset of a larger candidate set", "selective endpoint/marker display", "selective endpoint/marker display"),
  stringsAsFactors = FALSE
)
write.csv(panel_provenance, file.path(result_dir, "R2_6_figure5_panel_provenance.csv"), row.names = FALSE, na = "")

subtype_counts <- as.data.frame(table(dat$Subtype)); names(subtype_counts) <- c("Subtype", "N")
write.csv(subtype_counts, file.path(result_dir, "R2_6_subtype_counts.csv"), row.names = FALSE)

missingness <- data.frame(Feature = feature_order, N_nonmiss = vapply(feature_order, function(f) sum(is.finite(dat[[f]])), integer(1))) %>% mutate(N_missing = nrow(dat) - N_nonmiss, Missing_fraction = N_missing / nrow(dat))
write.csv(missingness, file.path(result_dir, "R2_6_feature_missingness.csv"), row.names = FALSE)

capture.output({
  cat("R2.6 complete variable-selection audit\n")
  cat("Clinical source copy:", normalizePath(input_clinical), "\n")
  cat("Surgery source copy:", normalizePath(input_surgery), "\n")
  cat("Patients:", nrow(dat), "\n")
  cat("Subtype counts:\n"); print(table(dat$Subtype))
  cat("Texture counts:\n"); print(table(dat$Texture, useNA = "ifany"))
  cat("Adhesion counts:\n"); print(table(dat$Adhesion, useNA = "ifany"))
  cat("Complete blood features:", length(feature_order), "\n")
  cat("Historical candidates:", paste(historical_candidates, collapse = ", "), "\n")
  cat("Historical rule: min(P_texture, P_adhesion) < 0.1\n")
  cat("Historical selected:", paste(hist_rows$Candidate[hist_rows$Selected_by_historical_MinP_lt_0_1], collapse = ", "), "\n")
  cat("Features with subtype KW FDR < 0.05:", sum(kw_rows$FDR_BH_83 < 0.05, na.rm = TRUE), "\n")
  cat("Features with surgical correlation global q < 0.05:", sum(corr_rows$Q_global_168_BH < 0.05, na.rm = TRUE), "of 168 tests\n")
}, file = file.path(result_dir, "R2_6_audit_summary.txt"))

# 6. Transparent overview figures for supplementary use.
corr_plot <- corr_rows %>%
  mutate(Feature = factor(Feature, levels = rev(corr_features)), Endpoint = factor(Endpoint, levels = c("Tumor_texture_ordinal", "Brainstem_adhesion")))
p_corr <- ggplot(corr_plot, aes(x = Endpoint, y = Feature, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(aes(size = -log10(pmax(Q_global_168_BH, 1e-300))), color = "black", alpha = 0.65) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  scale_size_continuous(name = "-log10(global q)", range = c(0.5, 5)) +
  labs(x = NULL, y = NULL, fill = "Spearman rho", title = "Complete surgical-phenotype screen: age + 83 blood variables") +
  theme_classic(base_size = 8) + theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 25, hjust = 1), plot.title = element_text(face = "bold", size = 10), legend.position = "right")
ggsave(file.path(figure_dir, "R2_6_complete_surgical_correlation_heatmap.pdf"), p_corr, width = 6.4, height = 15, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(figure_dir, "R2_6_complete_surgical_correlation_heatmap.png"), p_corr, width = 6.4, height = 15, dpi = 300)

kw_plot <- kw_rows %>% mutate(Feature = factor(Feature, levels = rev(feature_order)))
p_kw <- ggplot(kw_plot, aes(x = "Subtype KW", y = Feature, fill = -log10(pmax(FDR_BH_83, 1e-300)))) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(FDR_BH_83 < 0.05, "q<0.05", sprintf("q=%.3f", FDR_BH_83))), size = 1.7) +
  scale_fill_viridis_c(option = "C", name = "-log10 BH q") +
  labs(x = NULL, y = NULL, title = "Complete subtype screen: 83 variables") +
  theme_classic(base_size = 8) + theme(axis.text.y = element_text(size = 5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(face = "bold", size = 10), legend.position = "right")
ggsave(file.path(figure_dir, "R2_6_complete_subtype_screen_heatmap.pdf"), p_kw, width = 4.8, height = 15, device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(figure_dir, "R2_6_complete_subtype_screen_heatmap.png"), p_kw, width = 4.8, height = 15, dpi = 300)

cat("R2.6 complete audit finished. Results:", result_dir, "\n")
