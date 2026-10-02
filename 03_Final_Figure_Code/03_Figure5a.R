# Purpose: Final palette-aligned clinical annotation matrix
# Manuscript/panels: Figure 5a
# Reviewer comments: not a separate reviewer analysis
# Execution requires the inputs and working-directory layout described in README.md; "03_Final_Figure_Code/03_Figure5a.R"
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(svglite)
  library(ragg)
  library(patchwork)
})

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
out_dir <- dirname(normalizePath(script_arg, mustWork = TRUE))
data_file <- file.path(out_dir, "source_data", "Fig5a_clinical_matrix.csv")
if (!file.exists(data_file)) stop("Missing Fig5a_clinical_matrix.csv")

group_cols <- c(C1 = "#00A087", C2 = "#4DBBD5", C3 = "#E64B35")
clinical_cols <- c(
  # Match the Small/Large hues used in the Figure 5h size comparison.
  "Tumour size::Small" = "#91A3CE", "Tumour size::Large" = "#D28C8E",
  "KOOS::I" = "#E2E9F3", "KOOS::II" = "#B5C8DF", "KOOS::III" = "#82A0C2", "KOOS::IV" = "#51779F",
  "Texture::Soft" = "#B2E2E2", "Texture::Medium" = "#FEB24C", "Texture::Hard" = "#E31A1C",
  "Blood supply::No" = "#AABCAE", "Blood supply::Yes" = "#C6878F",
  # Match the No/Yes adhesion palette used in the Figure 5 adhesion panels.
  "Brainstem adhesion::No" = "#4DBBD5", "Brainstem adhesion::Yes" = "#E64B35",
  "Resection::GTR" = "#EAEAEA", "Resection::NTR" = "#D0D0D0", "Resection::STR" = "#A6A6A6", "Resection::PR" = "#7D7D7D"
)

dat <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(Group = factor(Group, levels = c("C1", "C2", "C3"))) %>%
  arrange(Group, SampleID) %>%
  mutate(order = row_number())
if (nrow(dat) != 38L || any(table(dat$Group) != c(13L, 12L, 13L))) {
  stop("Expected 38 patients, with group counts C1=13, C2=12, C3=13")
}
if (anyDuplicated(dat$SampleID)) stop("Duplicate sample IDs in source data")

long <- dat %>%
  select(Group, order, `Tumour size`, KOOS, Texture, `Blood supply`, `Brainstem adhesion`, Resection) %>%
  pivot_longer(-c(Group, order), names_to = "Feature", values_to = "Value") %>%
  mutate(
    Feature = factor(Feature, levels = c("Resection", "Brainstem adhesion", "Blood supply", "Texture", "KOOS", "Tumour size")),
    FillKey = paste(Feature, Value, sep = "::")
  )

group_bar <- dat %>%
  group_by(Group) %>%
  summarise(xmin = min(order) - 0.5, xmax = max(order) + 0.5, n = n(), .groups = "drop") %>%
  mutate(xmid = (xmin + xmax) / 2, label = paste0(Group, " (n=", n, ")"))

base_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(1, 2, 0, 2),
    text = element_text(colour = "black")
  )

p_group <- ggplot(group_bar) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = Group), colour = "white", linewidth = 0.35) +
  geom_text(aes(x = xmid, y = 0.5, label = label), family = "Arial", size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  scale_x_continuous(limits = c(0.5, 38.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text = element_blank(), legend.position = "none", plot.margin = margin(0, 2, 0, 2))

p_matrix <- ggplot(long, aes(x = order, y = Feature, fill = FillKey)) +
  geom_tile(width = 0.96, height = 0.94, colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = c(13.5, 25.5), colour = "white", linewidth = 1) +
  scale_fill_manual(values = clinical_cols, na.value = "#F5F5F5", guide = "none") +
  scale_x_continuous(limits = c(0.5, 38.5), expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(labels = c(
    "Resection" = "Resection", "Brainstem adhesion" = "Brainstem adhesion",
    "Blood supply" = "Blood supply", "Texture" = "Tumour texture",
    "KOOS" = "KOOS grade", "Tumour size" = "Tumour size"
  ), expand = expansion(add = 0.05)) +
  labs(x = "Patients ordered by transcriptional group", y = NULL) +
  base_theme +
  theme(
    axis.text.y = element_text(size = 6.6, margin = margin(r = 4)),
    axis.title.x = element_text(size = 6.5, margin = margin(t = 4)),
    plot.margin = margin(0, 2, 2, 2)
  )

panel <- p_group / p_matrix + plot_layout(heights = c(0.17, 1))

legend_data <- tibble::tribble(
  ~Feature, ~Value,
  "Tumour size", "Small", "Tumour size", "Large",
  "KOOS grade", "I", "KOOS grade", "II", "KOOS grade", "III", "KOOS grade", "IV",
  "Tumour texture", "Soft", "Tumour texture", "Medium", "Tumour texture", "Hard",
  "Blood supply", "No", "Blood supply", "Yes",
  "Brainstem adhesion", "No", "Brainstem adhesion", "Yes",
  "Resection", "GTR", "Resection", "NTR", "Resection", "STR", "Resection", "PR"
) %>%
  mutate(
    Feature = factor(Feature, levels = c("Tumour size", "KOOS grade", "Tumour texture", "Blood supply", "Brainstem adhesion", "Resection")),
    row = as.numeric(Feature),
    col = ave(as.numeric(Feature), Feature, FUN = seq_along),
    FillKey = paste(recode(as.character(Feature), "KOOS grade" = "KOOS", "Tumour texture" = "Texture"), Value, sep = "::")
  )
legend <- ggplot(legend_data, aes(col, row)) +
  geom_point(aes(fill = FillKey), shape = 22, size = 3, colour = "black", stroke = 0.2) +
  geom_text(aes(x = col + 0.18, label = Value), hjust = 0, family = "Arial", size = 2.1) +
  scale_fill_manual(values = clinical_cols, guide = "none") +
  scale_x_continuous(limits = c(0.7, 4.8), breaks = NULL, expand = c(0, 0)) +
  scale_y_reverse(breaks = 1:6, labels = levels(legend_data$Feature), expand = expansion(add = 0.45)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(axis.text.y = element_text(colour = "black", size = 6.3), plot.margin = margin(2, 5, 2, 2))

stem <- "Fig5a_patient_clinical_matrix_paletteAligned_20260930"
write_plot <- function(plot, name, width_mm, height_mm) {
  svg_file <- file.path(out_dir, paste0(name, ".svg"))
  pdf_file <- file.path(out_dir, paste0(name, ".pdf"))
  tiff_file <- file.path(out_dir, paste0(name, ".tiff"))
  png_file <- file.path(out_dir, paste0(name, "_preview.png"))
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite(svg_file, width = w, height = h, system_fonts = list(sans = "Arial")); print(plot); dev.off()
  # The local Cairo text backend cannot resolve Arial; Quartz preserves text in the PDF.
  grDevices::quartz(type = "pdf", file = pdf_file, width = w, height = h, family = "Arial")
  print(plot); dev.off()
  ragg::agg_tiff(tiff_file, width = w, height = h, units = "in", res = 600, background = "white"); print(plot); dev.off()
  ragg::agg_png(png_file, width = w, height = h, units = "in", res = 300, background = "white"); print(plot); dev.off()
}
write_plot(panel, stem, 178, 76)
write_plot(legend, "Fig5a_clinical_matrix_legend_paletteAligned_20260930", 142, 54)
write_csv(dat %>% select(-order), file.path(out_dir, "source_data", "Fig5a_clinical_matrix_paletteAligned_20260930.csv"))
writeLines(c(
  "Figure 5a palette-aligned redraw",
  "All 38 patients retained; group counts C1=13, C2=12, C3=13.",
  "Subtype strip uses locked colors C1 #00A087, C2 #4DBBD5, C3 #E64B35.",
  "Tumour-size colors match Figure 5h (Small #91A3CE; Large #D28C8E).",
  "Texture colors match Figures 5b/c; adhesion colors match Figures 5d/e.",
  "KOOS, blood-supply, and resection colors retain their original category palettes because no other main panel repeats these variables.",
  "No values or categories were changed; patient identifiers are not printed.",
  "SVG, PDF, TIFF, and PNG preview exported using R."
), file.path(out_dir, paste0(stem, "_QA.txt")))
