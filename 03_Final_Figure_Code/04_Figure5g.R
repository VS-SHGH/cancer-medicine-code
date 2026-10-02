# Purpose: Final pairwise-complete Spearman rho heatmap without stars or q values
# Manuscript/panels: Figure 5g
# Reviewer comments: R2 Comment 7
# Execution requires the inputs and working-directory layout described in README.md; "03_Final_Figure_Code/04_Figure5g.R"
#!/usr/bin/env Rscript

# Re-create the Figure 5g lower-triangle Spearman correlogram from the
# canonical 38-patient workbooks. This script writes to its own folder and
# never overwrites the submitted Figure 5.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_path) || !nzchar(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  stopifnot(length(script_arg) == 1L)
  script_path <- sub("^--file=", "", script_arg)
}
out_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
work_dir <- normalizePath(file.path(out_dir, ".."), mustWork = TRUE)
input_dir <- file.path(work_dir, "02_Code", "R2_6_variable_selection")
clinical_file <- file.path(input_dir, "input_cl_bulk_with_Subtype.xlsx")
surgery_file <- file.path(input_dir, "input_VS_surgery_analysis.xlsx")
stopifnot(file.exists(clinical_file), file.exists(surgery_file))

to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "--", "—")] <- NA_character_
  suppressWarnings(as.numeric(gsub("[<>]", "", x)))
}

clean_texture <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("软|Soft", x, ignore.case = TRUE) ~ "Soft",
    grepl("中|Medium", x, ignore.case = TRUE) ~ "Medium",
    grepl("硬|Hard", x, ignore.case = TRUE) ~ "Hard",
    TRUE ~ NA_character_
  )
}

clean_adhesion <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("无|否|No|^N$|^0$", x, ignore.case = TRUE) ~ 0,
    grepl("有|是|Yes|^Y$|^1$|丰富|富", x, ignore.case = TRUE) ~ 1,
    TRUE ~ NA_real_
  )
}

clinical <- read_excel(clinical_file, .name_repair = "minimal") |>
  filter(Subtype %in% c("C1", "C2", "C3"))
surgery <- read_excel(surgery_file, .name_repair = "minimal") |>
  select(SampleID, all_of(c("肿瘤质地", "脑干面粘连")))
stopifnot(nrow(clinical) == 38L, !anyDuplicated(clinical$SampleID),
          !anyDuplicated(surgery$SampleID), setequal(clinical$SampleID, surgery$SampleID))

dat <- clinical |>
  select(SampleID, Subtype, Age, `自然杀伤细胞`, `平均血小板体积`, `总T淋巴细胞`) |>
  left_join(surgery, by = "SampleID") |>
  transmute(
    SampleID,
    `Tumor Texture` = match(clean_texture(`肿瘤质地`), c("Soft", "Medium", "Hard")),
    Adhesion = clean_adhesion(`脑干面粘连`),
    `NK Cells` = to_numeric(`自然杀伤细胞`),
    Age = to_numeric(Age),
    MPV = to_numeric(`平均血小板体积`),
    `Total T Cells` = to_numeric(`总T淋巴细胞`)
  )

variable_order <- c("Tumor Texture", "Adhesion", "NK Cells", "Age", "MPV", "Total T Cells")
stopifnot(nrow(dat) == 38L, all(variable_order %in% names(dat)))

# Pairwise-complete Spearman tests; keep the pairwise N in the source-data table.
pair_rows <- list()
for (i in seq_along(variable_order)) {
  if (i == 1L) next
  for (j in seq_len(i - 1L)) {
    x_name <- variable_order[j]
    y_name <- variable_order[i]
    x <- dat[[x_name]]
    y <- dat[[y_name]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
      rho <- NA_real_
    } else {
      rho <- unname(suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$estimate))
    }
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      Row_variable = y_name,
      Column_variable = x_name,
      N_pairwise = sum(ok),
      Spearman_rho = rho,
      stringsAsFactors = FALSE
    )
  }
}
source_pairs <- bind_rows(pair_rows)
write.csv(source_pairs, file.path(out_dir, "Fig5g_pairwise_source_data.csv"), row.names = FALSE, na = "")

plot_data <- source_pairs |>
  mutate(
    x = factor(Column_variable, levels = variable_order),
    y = factor(Row_variable, levels = rev(variable_order)),
    rho_label = ifelse(is.na(Spearman_rho), "NA", sprintf("%.2f", Spearman_rho))
  )

p <- ggplot(plot_data, aes(x = x, y = y, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.45, width = 0.98, height = 0.98) +
  geom_text(aes(label = rho_label), family = "Helvetica", size = 3.0, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
    name = expression("Spearman " * rho), na.value = "grey90"
  ) +
  scale_x_discrete(expand = expansion(add = 0.02)) +
  scale_y_discrete(expand = expansion(add = 0.02)) +
  labs(
    tag = "g",
    x = NULL, y = NULL,
    title = "Exploratory correlations with surgical phenotypes",
    subtitle = "Pairwise-complete observations (n = 37–38)"
  ) +
  coord_fixed(clip = "off") +
  theme_minimal(base_size = 8.5, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8.0, angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(size = 8.0, colour = "black"),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 10.2, face = "bold", hjust = 0.5, colour = "black"),
    plot.tag = element_text(size = 11.0, face = "bold", hjust = 0, vjust = 1, colour = "black"),
    plot.tag.position = c(0, 1),
    plot.tag.location = "plot",
    plot.subtitle = element_text(size = 8.2, hjust = 0.5, colour = "black", margin = margin(t = 1, b = 4)),
    legend.title = element_text(size = 8.0),
    legend.text = element_text(size = 7.2),
    legend.key.height = grid::unit(9, "mm"),
    legend.position = "right",
    plot.margin = margin(8, 7, 6, 10)
  )

# Keep the panel proportions close to the Figure 5g footprint; text remains
# editable in SVG/PDF, with a 600-dpi PNG for quick review.
pdf_device <- function(filename, width, height, bg, ...) {
  # The local macOS R build advertises Cairo but its X11-linked Cairo DLL is
  # unavailable, so use base PDF there; use cairo_pdf on other supported hosts.
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    grDevices::pdf(filename, width = width, height = height, bg = bg,
                   family = "Helvetica", useDingbats = FALSE)
  } else {
    grDevices::cairo_pdf(filename, width = width, height = height, bg = bg,
                         family = "Helvetica")
  }
}
svg_device <- function(filename, width, height, bg, ...) {
  svglite::svglite(filename, width = width, height = height, bg = bg, ...)
}
ggsave(file.path(out_dir, "Figure5g_replacement.svg"), p, width = 4.45, height = 3.55,
       device = svg_device, bg = "white")
ggsave(file.path(out_dir, "Figure5g_replacement.pdf"), p, width = 4.45, height = 3.55,
       device = pdf_device, bg = "white")
ggsave(file.path(out_dir, "Figure5g_replacement.png"), p, width = 4.45, height = 3.55,
       dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(out_dir, "Figure5g_replacement.tiff"), p, width = 4.45, height = 3.55,
       dpi = 600, device = ragg::agg_tiff, bg = "white", compression = "lzw")

capture.output({
  cat("Rows in canonical cohort:", nrow(dat), "\n")
  cat("Missing values per plotted variable:\n")
  print(vapply(dat[variable_order], function(x) sum(!is.finite(x)), integer(1)))
  cat("Pairwise correlations (pairwise complete):\n")
  print(source_pairs, row.names = FALSE)
}, file = file.path(out_dir, "Fig5g_QA_notes.txt"))

message("Figure 5g replacement panel and source data exported to: ", out_dir)
