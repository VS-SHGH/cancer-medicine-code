# Purpose: Correctly sample-matched external MGS plot
# Manuscript/panels: Supplementary Figure 2d
# Reviewer comments: R1 Comment 5.2; R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "03_Final_Figure_Code/05_External_MGS_panel.R"
#!/usr/bin/env Rscript

# Refresh only the external-MGS replacement panel; all other replacement panels
# and manuscript/reviewer-response files are left untouched.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(svglite)
  library(ragg)
})

root <- file.path(Sys.getenv('VS_CM_REPO_ROOT'), 'workspace/03_Active_Revision_Work')
imported <- file.path(Sys.getenv('VS_CM_REPO_ROOT'), 'workspace/02_Imported_Reference_Assets')
input <- file.path(imported, "03_Existing_Results_Reference/verified_outputs/external_validation/external_mgs_figure_input.csv")
out <- file.path(root, "04_Figures/Replacement_Panels_Original_Style")
data_out <- file.path(out, "source_data")
dir.create(data_out, recursive = TRUE, showWarnings = FALSE)

ext_mgs <- read_csv(input, show_col_types = FALSE) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("GSE141801", "GSE39645")),
    Class = factor(Class, levels = c("Injury-like", "nmSC Core"))
  )
stopifnot(nrow(ext_mgs) == 88L,
          !anyDuplicated(paste(ext_mgs$Dataset, ext_mgs$Sample, sep = "::")),
          all(is.finite(ext_mgs$MGS)),
          all(ext_mgs$MGS >= 0 & ext_mgs$MGS <= 1))
write_csv(ext_mgs, file.path(data_out, "Fig4_external_MGS.csv"))

external_mgs_p <- vapply(levels(ext_mgs$Dataset), function(ds) {
  wilcox.test(MGS ~ Class, data = filter(ext_mgs, Dataset == ds), exact = FALSE)$p.value
}, numeric(1))
format_p <- function(p) {
  if (p < 0.001) {
    if (abs(p - 3.017523e-7) < 1e-12) return("Wilcoxon P = 3.02×10⁻⁷")
    return(sprintf("Wilcoxon P = %.2e", p))
  }
  sprintf("Wilcoxon P = %.2f", p)
}
pdat <- data.frame(
  Dataset = factor(c("GSE141801", "GSE39645"), levels = levels(ext_mgs$Dataset)),
  label = vapply(external_mgs_p, format_p, character(1))
)
class_cols <- c("Injury-like" = "#CD5C5C", "nmSC Core" = "#4682B4")
p <- ggplot(ext_mgs, aes(Class, MGS, fill = Class)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = "black", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 1.3,
             colour = "black", alpha = 0.82) +
  facet_wrap(~Dataset, nrow = 1) +
  geom_text(data = pdat, aes(x = 1.5, y = 1.04, label = label),
            inherit.aes = FALSE, size = 2.5, family = "Helvetica") +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_y_continuous(limits = c(-0.03, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "Independently fitted MGS (0–1)",
       title = "External MGS sensitivity") +
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text = element_text(face = "plain", size = 8),
    plot.title = element_text(face = "bold", size = 9, hjust = 0),
    panel.grid = element_blank()
  )

stem <- file.path(out, "Fig4_external_MGS_sensitivity")
w <- 146 / 25.4
h <- 72 / 25.4
svglite::svglite(paste0(stem, ".svg"), width = w, height = h,
                 system_fonts = list(sans = "Arial"))
print(p); dev.off()
write_vector_pdf <- function(path, plot) {
  cairo_ok <- tryCatch({
    grDevices::cairo_pdf(path, width = w, height = h, family = "Helvetica")
    print(plot)
    grDevices::dev.off()
    TRUE
  }, error = function(e) {
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
    if (file.exists(path)) unlink(path)
    FALSE
  })
  if (!cairo_ok) {
    grDevices::pdf(path, width = w, height = h, family = "Helvetica")
    print(plot)
    grDevices::dev.off()
  }
}
write_vector_pdf(paste0(stem, ".pdf"), p)
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in",
               res = 600, background = "white")
print(p); dev.off()
ragg::agg_png(paste0(stem, "_preview.png"), width = w, height = h, units = "in",
              res = 300, background = "white")
print(p); dev.off()
writeLines(c(
  "External MGS replacement panel refresh",
  "Style inherited from the existing replacement panel; only the corrected sample-keyed MGS values and data-derived Wilcoxon annotations changed.",
  sprintf("GSE141801 Wilcoxon P = %.8g", external_mgs_p[["GSE141801"]]),
  sprintf("GSE39645 Wilcoxon P = %.8g", external_mgs_p[["GSE39645"]]),
  "Sample unit: one GEO accession per point; boxplot centre is the median and hinges are the interquartile range."
), file.path(out, "Fig4_external_MGS_sensitivity_QA.txt"))
cat("Updated external MGS replacement panel.\n")
