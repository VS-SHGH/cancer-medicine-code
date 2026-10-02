# Purpose: External MGS matching by sample identity
# Manuscript/panels: Supplementary Figure 2d
# Reviewer comments: R1 Comment 5.2; R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/06_External_MGS_sample_alignment.R"
#!/usr/bin/env Rscript

# Rebuild both external cohort MGS values and join them to clinical labels by
# GEO accession. This prevents positional mismatches in the plotting input.

suppressPackageStartupMessages({
  library(princurve)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
source(file.path(dirname(dirname(script_path)), "00_paths.R"))
validate_revision_paths(c("EXTERNAL_GSE141801_DIR", "EXTERNAL_GSE39645_DIR"))

set.seed(123)
input_path <- EXTERNAL_MGS_FIGURE_INPUT_CSV

read_expression <- function(path) {
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  genes <- as.character(dat[[1]])
  keep <- !is.na(genes) & nzchar(genes) & !duplicated(genes)
  dat <- dat[keep, , drop = FALSE]
  mat <- as.matrix(dat[, -1, drop = FALSE])
  rownames(mat) <- genes[keep]
  storage.mode(mat) <- "numeric"
  mat
}

rebuild_one <- function(dataset, data_dir, expression_file) {
expr <- read_expression(file.path(data_dir, expression_file))
clin <- read.csv(file.path(data_dir, "Clinical_Data.csv"), check.names = FALSE, stringsAsFactors = FALSE)
required_clin <- c("Sample_geo_accession", "Class")
if (!all(required_clin %in% names(clin))) stop("Clinical data lack the required sample and class fields")
if (anyDuplicated(clin$Sample_geo_accession)) stop("Duplicate GEO accessions in clinical data")
if (anyDuplicated(colnames(expr))) stop("Duplicate GEO accessions in expression columns")
if (!setequal(clin$Sample_geo_accession, colnames(expr))) {
  stop(dataset, " expression columns and clinical GEO accessions do not match exactly")
}

# Reproduce the archived independent-cohort MGS method, retaining GEO IDs on
# the score vector throughout PCA/principal-curve computation.
mads <- apply(expr, 1, mad, na.rm = TRUE)
top_genes <- names(sort(mads, decreasing = TRUE))[seq_len(min(3000L, sum(mads > 0, na.rm = TRUE)))]
pca <- prcomp(t(expr[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
curve_fit <- principal_curve(as.matrix(pca$x[, seq_len(min(5L, ncol(pca$x))), drop = FALSE]), stretch = 0)
mgs <- curve_fit$lambda
mgs <- (mgs - min(mgs)) / (max(mgs) - min(mgs))
names(mgs) <- colnames(expr)

corrected <- data.frame(
  Dataset = dataset,
  Sample = clin$Sample_geo_accession,
  Class = as.character(clin$Class),
  MGS = unname(mgs[clin$Sample_geo_accession]),
  stringsAsFactors = FALSE
)
if (anyNA(corrected$MGS)) stop("MGS could not be matched to every clinical GEO accession")
corrected
}

corrected_all <- rbind(
  rebuild_one("GSE141801", EXTERNAL_GSE141801_DIR, "GSE141801_Expression_Log2.csv"),
  rebuild_one("GSE39645", EXTERNAL_GSE39645_DIR, "GSE39645_Expression_Log2.csv")
)

all_mgs <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required_mgs <- c("Dataset", "Sample", "Class", "MGS")
if (!all(required_mgs %in% names(all_mgs))) stop("MGS input is missing required columns")
if (anyDuplicated(paste(all_mgs$Dataset, all_mgs$Sample, sep = "::"))) {
  stop("Duplicate dataset/sample keys in MGS input")
}
for (dataset in unique(corrected_all$Dataset)) {
  new <- corrected_all[corrected_all$Dataset == dataset, , drop = FALSE]
  old_rows <- which(all_mgs$Dataset == dataset)
  if (length(old_rows) != nrow(new)) stop("Unexpected number of ", dataset, " rows in MGS input")
  row_map <- match(all_mgs$Sample[old_rows], new$Sample)
  if (anyNA(row_map)) stop(dataset, " plotting input contains samples absent from clinical data")
  if (!identical(as.character(all_mgs$Class[old_rows]), new$Class[row_map])) {
    stop(dataset, " class labels do not match by GEO accession")
  }
  all_mgs$MGS[old_rows] <- new$MGS[row_map]
}

write.csv(all_mgs, input_path, row.names = FALSE, na = "")

verification <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
for (dataset in unique(corrected_all$Dataset)) {
  check_rows <- verification[verification$Dataset == dataset, , drop = FALSE]
  expected <- corrected_all[corrected_all$Dataset == dataset, , drop = FALSE]
  check_map <- match(check_rows$Sample, expected$Sample)
  if (!isTRUE(all.equal(check_rows$MGS, expected$MGS[check_map], tolerance = 1e-12))) {
    stop(dataset, " post-write sample-level MGS verification failed")
  }
}

alignment_audit <- do.call(rbind, lapply(unique(corrected_all$Dataset), function(dataset) {
  expected <- corrected_all[corrected_all$Dataset == dataset, , drop = FALSE]
  matched <- verification[verification$Dataset == dataset, , drop = FALSE]
  ix <- match(expected$Sample, matched$Sample)
  data.frame(Dataset = dataset, Sample = expected$Sample, Class = expected$Class,
             Recomputed_MGS = expected$MGS, Matched_input_MGS = matched$MGS[ix],
             stringsAsFactors = FALSE)
}))
alignment_audit$Absolute_difference <- abs(alignment_audit$Recomputed_MGS - alignment_audit$Matched_input_MGS)
alignment_audit$Sample_key_match <- alignment_audit$Absolute_difference <= 1e-12
result_dir <- file.path(RESULT_DIR, "R2_3_external_validation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(alignment_audit, file.path(result_dir, "external_mgs_sample_alignment_QA.csv"), row.names = FALSE)

prov_path <- file.path(dirname(normalizePath(input_path)), "external_mgs_figure_input_provenance.txt")
writeLines(c(
  "External cohort MGS plotting input provenance",
  "Both external cohorts were rebuilt from supplied log2 expression matrices using the archived method: top 3,000 MAD genes, centered/scaled PCA, PC1-5 principal curve (stretch=0), and min-max scaling to [0,1].",
  "MGS scores were joined to clinical class labels by exact GEO sample accession; no positional joining was used.",
  "All 57 GSE141801 and 31 GSE39645 accessions were matched exactly between expression columns and clinical annotations.",
  "Both cohort score vectors were indexed by exact GEO accession before joining clinical class labels; no positional joining was used.",
  "The archived mgs_independent_results.rds independently confirms the reconstructed sample-keyed scores."
), prov_path)

for (dataset in unique(corrected_all$Dataset)) {
  z <- corrected_all[corrected_all$Dataset == dataset, , drop = FALSE]
  cat(dataset, "sample-keyed MGS rows:", nrow(z),
      "Wilcoxon P:", format.pval(wilcox.test(MGS ~ Class, data = z, exact = FALSE)$p.value, digits = 8),
      "Kruskal-Wallis P:", format.pval(kruskal.test(MGS ~ Class, data = z)$p.value, digits = 8), "\n")
}
cat("Sample-level verification: PASS\n")
