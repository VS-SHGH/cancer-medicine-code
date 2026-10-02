# Purpose: Recruitment/checkpoint programme scores and C3 purity-adjusted correlation
# Manuscript/panels: Figure 2d,e; Supplementary Table 6
# Reviewer comments: R1 Comment 5.6; R2 Comment 3
# Execution requires the inputs and working-directory layout described in README.md; "01_Main_Analyses/05_Immune_programs_C3_correlation.R"
# Reviewer 1 focused audit: C3-only immune-program correlations
# This script reproduces the scoring convention used in VS.4/code/step7_immune.R
# and does not modify the manuscript or response letter.

suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

rdata_file <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "inputs/legacy_project/RData/Step2_Clustering.RData")
purity_file <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "workspace/02_Imported_Reference_Assets/03_Existing_Results_Reference/verified_outputs/estimate_official/estimate_official_scores.csv")
score_file <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "inputs/prior_project/supplement/all/step8_Immune/VS_score_summary.csv")
work_dir <- file.path(Sys.getenv("VS_CM_REPO_ROOT"), "workspace/03_Active_Revision_Work")
result_dir <- file.path(work_dir, "03_Results", "R1_reviewer1_audit")
figure_dir <- file.path(work_dir, "04_Figures", "R1_reviewer1_audit")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

load(rdata_file)
stopifnot(exists("clean_tpm"), exists("meta_step2"), all(c("SampleID", "Subtype") %in% names(meta_step2)))
meta_step2$Subtype <- factor(meta_step2$Subtype, levels = c("C1", "C2", "C3"))
common <- intersect(colnames(clean_tpm), meta_step2$SampleID)
stopifnot(length(common) == 38L, !anyDuplicated(common))
clean_tpm <- clean_tpm[, common, drop = FALSE]
meta_step2 <- meta_step2[match(common, meta_step2$SampleID), , drop = FALSE]

recruitment_genes <- c("CCL4", "CCL5", "CXCL9", "CXCL10")
suppression_genes <- c("PDCD1", "CD274", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "VSIR", "CD96", "SIGLEC7")
checkpoint_genes <- suppression_genes
recruitment_present <- intersect(recruitment_genes, rownames(clean_tpm))
suppression_present <- intersect(suppression_genes, rownames(clean_tpm))
stopifnot(length(recruitment_present) >= 3, length(suppression_present) >= 3)

# Reproduce Section 32 exactly: recruitment and suppression were scored together
# in a separate two-gene-set ssGSEA run. Checkpoint_Score was generated earlier
# in the historical multi-gene-set run, so it is read from the frozen score file.
gene_sets <- list(
  Immune_Recruitment_Score = recruitment_present,
  Immune_Suppression_Score = suppression_present
)
score_mat <- GSVA::gsva(as.matrix(clean_tpm), gene_sets, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
score_df <- data.frame(SampleID = colnames(score_mat), t(score_mat), check.names = FALSE)
score_df <- score_df %>% left_join(meta_step2 %>% select(SampleID, Subtype), by = "SampleID")

purity <- read_csv(purity_file, show_col_types = FALSE) %>% select(SampleID, Tumour_Purity)
historical_scores <- read_csv(score_file, show_col_types = FALSE) %>%
  select(SampleID, Subtype_historical = Subtype, Checkpoint_Score)
stopifnot(nrow(historical_scores) == 38L, setequal(historical_scores$SampleID, score_df$SampleID))
score_df <- score_df %>% left_join(purity, by = "SampleID") %>%
  left_join(historical_scores, by = "SampleID")
stopifnot(!anyNA(score_df$Tumour_Purity), !anyNA(score_df$Checkpoint_Score),
          as.character(score_df$Subtype) == score_df$Subtype_historical)
write_csv(score_df, file.path(result_dir, "R1_immune_program_scores_all38.csv"))

# Program summaries by subtype (descriptive; all 38 samples)
sub_sum <- score_df %>% group_by(Subtype) %>% summarise(
  n = n(),
  Recruitment_mean = mean(Immune_Recruitment_Score), Recruitment_median = median(Immune_Recruitment_Score),
  Suppression_mean = mean(Immune_Suppression_Score), Suppression_median = median(Immune_Suppression_Score),
  Checkpoint_mean = mean(Checkpoint_Score), Checkpoint_median = median(Checkpoint_Score),
  .groups = "drop"
)
write_csv(sub_sum, file.path(result_dir, "R1_immune_program_subtype_summary.csv"))

# Correlation helper: the manuscript's original convention is Spearman, exact=FALSE.
cor_one <- function(df, x, y, scope, purity_adjusted = FALSE) {
  d <- df %>% select(SampleID, all_of(c(x, y)), any_of("Tumour_Purity")) %>% filter(complete.cases(across(all_of(c(x, y)))))
  if (purity_adjusted) {
    d <- d %>% filter(!is.na(Tumour_Purity))
    d$x_adj <- resid(lm(reformulate("Tumour_Purity", response = x), data = d))
    d$y_adj <- resid(lm(reformulate("Tumour_Purity", response = y), data = d))
    xx <- d$x_adj; yy <- d$y_adj
  } else { xx <- d[[x]]; yy <- d[[y]] }
  ok <- complete.cases(xx, yy)
  has_ties <- anyDuplicated(xx[ok]) > 0L || anyDuplicated(yy[ok]) > 0L
  s_asym <- suppressWarnings(cor.test(xx, yy, method = "spearman", exact = FALSE))
  s_exact <- suppressWarnings(cor.test(xx, yy, method = "spearman", exact = TRUE))
  p <- suppressWarnings(cor.test(xx, yy, method = "pearson"))
  data.frame(Scope = scope, X = x, Y = y, Purity_adjusted = purity_adjusted,
             N = length(xx), Spearman_Rho = unname(s_asym$estimate),
             P_Spearman_asymptotic = s_asym$p.value,
             P_Spearman_exact_TRUE_request = s_exact$p.value,
             Pearson_R = unname(p$estimate), P_Pearson = p$p.value,
             Spearman_exact_TRUE_method = if (has_ties) "asymptotic fallback because of ties" else if (sum(ok) < 10L) "exact" else "AS89/Edgeworth approximation",
             stringsAsFactors = FALSE)
}

c3 <- score_df %>% filter(Subtype == "C3")
rows <- list(
  cor_one(c3, "Immune_Recruitment_Score", "Immune_Suppression_Score", "C3_only", FALSE),
  cor_one(c3, "Immune_Recruitment_Score", "Checkpoint_Score", "C3_only", FALSE),
  cor_one(c3, "Immune_Suppression_Score", "Checkpoint_Score", "C3_only", FALSE),
  cor_one(c3, "Immune_Recruitment_Score", "Immune_Suppression_Score", "C3_only", TRUE),
  cor_one(c3, "Immune_Recruitment_Score", "Checkpoint_Score", "C3_only", TRUE),
  cor_one(c3, "Immune_Suppression_Score", "Checkpoint_Score", "C3_only", TRUE)
)
cor_out <- bind_rows(rows) %>% mutate(
  Endpoint_role = ifelse(X == "Immune_Suppression_Score" & Y == "Checkpoint_Score",
                         "same-gene-set concordance check; not an independent hypothesis",
                         "exploratory association")
)
write_csv(cor_out, file.path(result_dir, "R1_C3_internal_program_correlations.csv"))

# Median-split quadrants exactly as the historical exploratory output.
quad <- c3 %>% mutate(
  Recruitment_group = ifelse(Immune_Recruitment_Score >= median(Immune_Recruitment_Score), "High", "Low"),
  Suppression_group = ifelse(Immune_Suppression_Score >= median(Immune_Suppression_Score), "High", "Low"),
  Quadrant = paste(Recruitment_group, "recruit /", Suppression_group, "suppress")
) %>% count(Quadrant)
write_csv(quad, file.path(result_dir, "R1_C3_recruitment_suppression_quadrants.csv"))

# Multiplicity inventory for Reviewer 1's statistical reporting audit.
inv <- data.frame(
  Analysis_family = c("R2.6 subtype omnibus screen", "R2.6 clinical correlations", "R2.6 pairwise Dunn tests", "CIBERSORT LM22 subtype screen", "MCP-counter subtype screen", "CIBERSORT-MCP correlations", "C3 internal program correlations (2 executed; 1 unique)", "Surgical phenotype Fisher tests"),
  Number_of_primary_tests = c(83, 168, 249, 22, 10, 4, 2, 2),
  Primary_adjustment_or_convention = c("BH across 83 omnibus tests", "BH across 168 correlations", "BH across 249 pairwise tests", "BH across 22 cell types; filtered sensitivity reported separately", "BH across 10 cell types; pairwise global BH across 30", "Correlation methods harmonized; no family-wide claim", "Two executed C3-only Spearman comparisons reduce to one unique association because suppression and checkpoint use the same genes; n=13; purity-adjusted sensitivity", "Exact Fisher; model-based size-adjusted confirmation"),
  stringsAsFactors = FALSE
)
write_csv(inv, file.path(result_dir, "R1_multiplicity_inventory.csv"))

# Compact reviewer-facing scatterplot for the key C3 checkpoint-recruitment result.
p <- ggplot(c3, aes(Immune_Recruitment_Score, Checkpoint_Score)) +
  geom_point(size = 2.5, colour = "#3C5488") +
  geom_smooth(method = "lm", se = TRUE, colour = "#3C5488", linewidth = 0.7) +
  labs(x = "Immune recruitment ssGSEA score", y = "Checkpoint ssGSEA score", title = "C3-only recruitment–checkpoint association") +
  theme_classic(base_size = 11)
ggsave(file.path(figure_dir, "R1_C3_recruitment_checkpoint_scatter.pdf"), p, width = 5.2, height = 4.2, useDingbats = FALSE)

txt <- c(
  "Reviewer 1 focused audit: C3-only immune-program correlations",
  paste0("Recruitment genes present: ", paste(recruitment_present, collapse = ", ")), 
  paste0("Suppression/checkpoint genes present: ", paste(suppression_present, collapse = ", ")), 
  "Suppression and checkpoint use the same genes but come from separate ssGSEA runs; they are not independent endpoints.",
  paste0("C3 samples: ", nrow(c3)),
  paste0("Key recruitment–checkpoint Spearman rho (exact=FALSE): ", sprintf("%.6f", cor_out$Spearman_Rho[cor_out$X == "Immune_Recruitment_Score" & cor_out$Y == "Checkpoint_Score" & !cor_out$Purity_adjusted]), "; P=", format(cor_out$P_Spearman_asymptotic[cor_out$X == "Immune_Recruitment_Score" & cor_out$Y == "Checkpoint_Score" & !cor_out$Purity_adjusted], scientific = TRUE)),
  "Interpretation: this is an exploratory within-C3 correlation with n=13; it should not be presented as independent confirmation without explicitly stating the small sample and multiplicity context."
)
writeLines(txt, file.path(result_dir, "R1_C3_program_correlation_summary.txt"))
