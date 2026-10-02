# Purpose: Discovery one-versus-rest limma-trend contrasts and 100-gene signatures
# Manuscript/panels: Figure 1c; Supplementary Table 2
# Reviewer comments: not a separate reviewer analysis
# Execution requires the inputs and working-directory layout described in README.md; "01_Main_Analyses/02_Differential_expression.R"
# ==============================================================================
# Step 3: 差异表达分析 (limma-trend, One-vs-Rest)
# ==============================================================================
# 注意: 当前数据仅有 TPM, 无原始 counts, 故使用 limma-trend
# 理想方案为对 counts 使用 limma-voom/DESeq2, 请在文章 Methods 中注明

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

source(Sys.getenv("VS_CONFIG_R"))

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
})

setwd(PROJECT_DIR)

# ==============================================================================
# 1. 加载 Step 2 数据
# ==============================================================================
load(file.path(RDATA_DIR, "Step2_Clustering.RData"))

# 确保对齐
clean_tpm <- clean_tpm[, meta_step2$SampleID]
stopifnot(all(colnames(clean_tpm) == meta_step2$SampleID))

message("✅ 数据加载: ", ncol(clean_tpm), " 样本")

# ==============================================================================
# 2. limma-trend 差异分析
# ==============================================================================
log_tpm <- log2(clean_tpm + 1)

# 设计矩阵
group <- factor(meta_step2$Subtype, levels = c("C1", "C2", "C3"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- meta_step2$SampleID

# 拟合
fit <- lmFit(log_tpm, design)

# 对比矩阵: One-vs-Rest
contrasts_mat <- makeContrasts(
  C1_vs_Others = C1 - (C2 + C3) / 2,
  C2_vs_Others = C2 - (C1 + C3) / 2,
  C3_vs_Others = C3 - (C1 + C2) / 2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2, trend = TRUE)  # trend=TRUE 适配连续值数据

# ==============================================================================
# 3. 提取结果
# ==============================================================================
dea_dir <- file.path(TABLE_DIR, "Step3_DEA")
if (!dir.exists(dea_dir)) dir.create(dea_dir, recursive = TRUE)

extract_deg <- function(fit_obj, coef_name) {
  res <- topTable(fit_obj, coef = coef_name, number = Inf, adjust.method = "BH")
  res$Gene <- rownames(res)
  res$Change <- case_when(
    res$adj.P.Val < PARAMS$dea_fdr_cutoff & res$logFC >  PARAMS$dea_logfc_cutoff ~ "Up",
    res$adj.P.Val < PARAMS$dea_fdr_cutoff & res$logFC < -PARAMS$dea_logfc_cutoff ~ "Down",
    TRUE ~ "Stable"
  )
  return(res)
}

deg_c1 <- extract_deg(fit2, "C1_vs_Others")
deg_c2 <- extract_deg(fit2, "C2_vs_Others")
deg_c3 <- extract_deg(fit2, "C3_vs_Others")

# 统计
for (nm in c("C1", "C2", "C3")) {
  deg <- get(paste0("deg_c", gsub("C", "", nm)))
  message(nm, ": Up=", sum(deg$Change == "Up"),
          " Down=", sum(deg$Change == "Down"),
          " Stable=", sum(deg$Change == "Stable"))
}

# 保存表格
write.csv(deg_c1, file.path(dea_dir, "DEG_C1_vs_Others.csv"), row.names = FALSE)
write.csv(deg_c2, file.path(dea_dir, "DEG_C2_vs_Others.csv"), row.names = FALSE)
write.csv(deg_c3, file.path(dea_dir, "DEG_C3_vs_Others.csv"), row.names = FALSE)

# ==============================================================================
# 4. 提取签名基因集
# ==============================================================================
n <- PARAMS$dea_top_n_markers

top_c1 <- deg_c1 %>% filter(Change == "Up") %>%
  slice_max(order_by = logFC, n = n) %>% pull(Gene)
top_c2 <- deg_c2 %>% filter(Change == "Up") %>%
  slice_max(order_by = logFC, n = n) %>% pull(Gene)
top_c3 <- deg_c3 %>% filter(Change == "Up") %>%
  slice_max(order_by = logFC, n = n) %>% pull(Gene)

signature_genes <- c(top_c1, top_c2, top_c3)
# 若有重复, 保留但不 unique (保持分组顺序)

message("🧬 签名基因: C1=", length(top_c1),
        " C2=", length(top_c2), " C3=", length(top_c3))

# ==============================================================================
# 5. 可视化: 签名基因热图
# ==============================================================================
fig3_dir <- file.path(FIG_DIR, "Step3_DEA")

meta_ordered <- meta_step2[order(meta_step2$Subtype), ]
expr_plot    <- log_tpm[signature_genes, meta_ordered$SampleID]

ann_col <- data.frame(Subtype = meta_ordered$Subtype)
rownames(ann_col) <- meta_ordered$SampleID
ann_colors <- list(Subtype = SUBTYPE_COLORS[sort(unique(meta_ordered$Subtype))])

gap_cols <- cumsum(table(meta_ordered$Subtype))
gap_cols <- gap_cols[-length(gap_cols)]
gap_rows <- c(length(top_c1), length(top_c1) + length(top_c2))

break_list <- seq(-2, 2, length.out = 101)

pdf(file.path(fig3_dir, "Signature_Heatmap.pdf"), width = 8, height = 10)
pheatmap(expr_plot,
         scale            = "row",
         breaks           = break_list,
         color            = HEATMAP_COLORS,
         cluster_cols     = FALSE,
         cluster_rows     = FALSE,
         gaps_col         = gap_cols,
         gaps_row         = gap_rows,
         show_colnames    = FALSE,
         show_rownames    = TRUE,
         fontsize_row     = 7,
         annotation_col   = ann_col,
         annotation_colors = ann_colors,
         main             = paste0("Top ", n, " Marker Genes per Subtype"),
         border_color     = NA)
dev.off()

# ==============================================================================
# 6. 可视化: 火山图 (可选, 用于补充材料)
# ==============================================================================
plot_volcano <- function(deg_df, subtype_name, colour) {
  ggplot(deg_df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(colour = Change), size = 0.8, alpha = 0.6) +
    scale_colour_manual(
      values = c("Up" = colour, "Down" = "steelblue", "Stable" = "grey80")
    ) +
    geom_hline(yintercept = -log10(PARAMS$dea_fdr_cutoff),
               linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-PARAMS$dea_logfc_cutoff, PARAMS$dea_logfc_cutoff),
               linetype = "dashed", colour = "grey40") +
    labs(
      x     = expression(log[2]~Fold~Change),
      y     = expression(-log[10]~FDR),
      title = paste0(subtype_name, " vs Others")
    ) +
    theme_publication(base_size = 12) +
    theme(legend.position = "none")
}

p_vol_c1 <- plot_volcano(deg_c1, "C1", SUBTYPE_COLORS["C1"])
p_vol_c2 <- plot_volcano(deg_c2, "C2", SUBTYPE_COLORS["C2"])
p_vol_c3 <- plot_volcano(deg_c3, "C3", SUBTYPE_COLORS["C3"])

p_vol_all <- cowplot::plot_grid(p_vol_c1, p_vol_c2, p_vol_c3,
                                ncol = 3, labels = c("C1", "C2", "C3"))

save_half_vector(p_vol_all,
                 file.path(fig3_dir, "Volcano_All_Subtypes.pdf"),
                 width = 14, height = 4.5)

# ==============================================================================
# 7. 保存
# ==============================================================================
save(deg_c1, deg_c2, deg_c3, signature_genes, top_c1, top_c2, top_c3,
     file = file.path(RDATA_DIR, "Step3_DEA_Results.RData"))

message("🎉 Step 3 完成!")
