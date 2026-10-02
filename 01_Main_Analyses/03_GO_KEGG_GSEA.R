# Purpose: Discovery GO/KEGG ORA and GSEA at the retained primary thresholds
# Manuscript/panels: Figure 1d–f
# Reviewer comments: not a separate reviewer analysis
# Execution requires the inputs and working-directory layout described in README.md; "01_Main_Analyses/03_GO_KEGG_GSEA.R"
# ==============================================================================
# Step 4: 功能富集分析 (GO / KEGG / GSEA)
# ==============================================================================

rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

source(Sys.getenv("VS_CONFIG_R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(GSVA)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(pheatmap)
  library(ggpubr)
})

setwd(PROJECT_DIR)

# ==============================================================================
# 1. 加载数据
# ==============================================================================
load(file.path(RDATA_DIR, "Step3_DEA_Results.RData"))
load(file.path(RDATA_DIR, "Step2_Clustering.RData"))

message("✅ 数据加载成功")

# ==============================================================================
# 2. 辅助函数定义
# ==============================================================================

# --- Top N 基因提取 (ORA 用) ---
prepare_topn_genes <- function(deg_table, group_name,
                               top_n = PARAMS$enrich_top_n) {
  sig <- deg_table %>%
    filter(adj.P.Val < PARAMS$dea_fdr_cutoff & logFC > 0) %>%
    arrange(desc(logFC)) %>%
    head(top_n)
  
  ids <- bitr(sig$Gene, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  message(group_name, ": ", nrow(ids), " 个基因成功转换 ID")
  return(ids$ENTREZID)
}

# --- GSEA 排序列表构建 ---
build_gsea_vector <- function(deg_table) {
  df <- deg_table %>%
    arrange(desc(logFC)) %>%
    dplyr::select(Gene, logFC)
  ids <- bitr(df$Gene, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  df_m <- merge(df, ids, by.x = "Gene", by.y = "SYMBOL")
  df_m <- df_m[!duplicated(df_m$ENTREZID), ]
  vec  <- setNames(df_m$logFC, df_m$ENTREZID)
  sort(vec, decreasing = TRUE)
}

# --- 安全的 GO 富集 ---
run_go_safe <- function(gene_ids) {
  if (length(gene_ids) < 10) return(NULL)
  enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff = 0.2, readable = TRUE)
}

# --- 安全的 KEGG 富集 ---
run_kegg_safe <- function(gene_ids) {
  if (length(gene_ids) < 10) return(NULL)
  tryCatch(
    enrichKEGG(gene = gene_ids, organism = "hsa", pvalueCutoff = 0.05),
    error = function(e) {
      warning("KEGG 查询失败 (网络问题?): ", e$message)
      return(NULL)
    }
  )
}

# --- 安全的 GSEA ---
run_gsea_safe <- function(gsea_vec, label) {
  message("🚀 GSEA: ", label)
  gseGO(geneList = gsea_vec, OrgDb = org.Hs.eg.db, ont = "BP",
        minGSSize = 10, maxGSSize = 500,
        pvalueCutoff = PARAMS$gsea_pval_cutoff,
        verbose = FALSE, seed = PARAMS$global_seed)
}

# ==============================================================================
# 3. ORA: GO + KEGG (Top N 策略)
# ==============================================================================
genes_c1 <- prepare_topn_genes(deg_c1, "C1")
genes_c2 <- prepare_topn_genes(deg_c2, "C2")
genes_c3 <- prepare_topn_genes(deg_c3, "C3")

go_c1 <- run_go_safe(genes_c1)
go_c2 <- run_go_safe(genes_c2)
go_c3 <- run_go_safe(genes_c3)

kegg_c1 <- run_kegg_safe(genes_c1)
kegg_c2 <- run_kegg_safe(genes_c2)
kegg_c3 <- run_kegg_safe(genes_c3)

# 绘制并保存 dotplot
enrich_fig_dir <- file.path(FIG_DIR, "Step4_Enrichment")

plot_save_dotplot <- function(enrich_obj, filename, title_txt, colour) {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) return(invisible(NULL))
  p <- dotplot(enrich_obj, showCategory = 10, title = title_txt) +
    scale_colour_gradient(low = colour, high = "grey20") +
    theme_publication(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
  save_pdf(p, file.path(enrich_fig_dir, filename), width = 8, height = 6)
}

plot_save_dotplot(go_c1, "GO_C1.pdf", "C1: GO Biological Process", SUBTYPE_COLORS["C1"])
plot_save_dotplot(go_c2, "GO_C2.pdf", "C2: GO Biological Process", SUBTYPE_COLORS["C2"])
plot_save_dotplot(go_c3, "GO_C3.pdf", "C3: GO Biological Process", SUBTYPE_COLORS["C3"])

plot_save_dotplot(kegg_c1, "KEGG_C1.pdf", "C1: KEGG Pathway", SUBTYPE_COLORS["C1"])
plot_save_dotplot(kegg_c2, "KEGG_C2.pdf", "C2: KEGG Pathway", SUBTYPE_COLORS["C2"])
plot_save_dotplot(kegg_c3, "KEGG_C3.pdf", "C3: KEGG Pathway", SUBTYPE_COLORS["C3"])

# ==============================================================================
# 4. GSEA 分析 (全部亚型)
# ==============================================================================
vec_c1 <- build_gsea_vector(deg_c1)
vec_c2 <- build_gsea_vector(deg_c2)
vec_c3 <- build_gsea_vector(deg_c3)

gse_c1 <- run_gsea_safe(vec_c1, "C1")
gse_c2 <- run_gsea_safe(vec_c2, "C2")
gse_c3 <- run_gsea_safe(vec_c3, "C3")

# --- GSEA 山峰图 ---
plot_gsea_mountain <- function(gse_obj, subtype_key, top_k = 3) {
  if (is.null(gse_obj) || nrow(gse_obj@result) == 0) return(invisible(NULL))
  
  pos_df <- gse_obj@result %>% filter(NES > 0) %>% arrange(p.adjust, desc(NES))
  if (nrow(pos_df) == 0) {
    message("⚠️ ", subtype_key, " 无 NES>0 的正向通路")
    return(invisible(NULL))
  }
  
  ids <- head(pos_df$ID, top_k)
  base_col <- SUBTYPE_COLORS[subtype_key]
  cols <- scales::alpha(base_col, c(1, 0.7, 0.45))[seq_along(ids)]
  
  p <- gseaplot2(gse_obj, geneSetID = ids,
                 title = paste0(subtype_key, ": Top Activated Pathways"),
                 color = cols, pvalue_table = TRUE,
                 ES_geom = "line", subplots = 1:2)
  
  # 调整字体
  p[[1]] <- p[[1]] + theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14),
    axis.text    = element_text(size = 11),
    legend.text  = element_text(size = 11)
  )
  p[[2]] <- p[[2]] + theme(
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 11)
  )
  
  outfile <- file.path(enrich_fig_dir, paste0("GSEA_Mountain_", subtype_key, ".pdf"))
  ggsave(outfile, p, width = 10, height = 6, device = "pdf")
  message("✅ ", outfile)
}

plot_gsea_mountain(gse_c1, "C1")
plot_gsea_mountain(gse_c2, "C2")
plot_gsea_mountain(gse_c3, "C3")

# --- GSEA 总览气泡图 ---
merge_gsea_results <- function(gse_list, names_list) {
  dfs <- lapply(seq_along(gse_list), function(i) {
    obj <- gse_list[[i]]
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    df <- obj@result
    df$Subtype <- names_list[i]
    df
  })
  bind_rows(dfs)
}

all_gsea <- merge_gsea_results(
  list(gse_c1, gse_c2, gse_c3),
  c("C1: Proliferative", "C2: Mesenchymal", "C3: Immune")
)

if (nrow(all_gsea) > 0) {
  top_paths <- all_gsea %>%
    filter(NES > 0) %>%
    group_by(Subtype) %>%
    slice_max(order_by = NES, n = 5) %>%
    pull(Description) %>%
    unique()
  
  plot_gsea_df <- all_gsea %>%
    filter(Description %in% top_paths) %>%
    mutate(
      Subtype = factor(Subtype, levels = c("C1: Proliferative",
                                           "C2: Mesenchymal",
                                           "C3: Immune")),
      Description = str_wrap(Description, width = 45)
    )
  
  # 锁定 Y 轴顺序
  y_order <- plot_gsea_df %>%
    arrange(Subtype, NES) %>%
    pull(Description) %>%
    unique()
  plot_gsea_df$Description <- factor(plot_gsea_df$Description, levels = y_order)
  
  p_summary <- ggplot(plot_gsea_df, aes(x = Subtype, y = Description)) +
    geom_point(aes(size = NES, colour = p.adjust)) +
    scale_colour_gradientn(
      colours = c("#E64B35", "#F39B7F", "grey80"),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size(range = c(3, 8)) +
    labs(x = NULL, y = NULL, colour = "FDR", size = "NES",
         title = "Functional Landscape of VS Subtypes") +
    theme_publication(base_size = 12) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  
  save_pdf(p_summary,
           file.path(enrich_fig_dir, "GSEA_Summary_Dotplot.pdf"),
           width = 9, height = 8)
}


# 
# # 诊断各组有多少条通路能通过各种筛选
# cat("=== C1 ===\n")
# cat("总通路数:", nrow(gse_c1@result), "\n")
# cat("p.adjust < 0.05:", sum(gse_c1@result$p.adjust < 0.05, na.rm=TRUE), "\n")
# cat("NES > 0 & p.adjust < 0.05:", 
#     sum(gse_c1@result$NES > 0 & gse_c1@result$p.adjust < 0.05, na.rm=TRUE), "\n")
# cat("NES > 0 (任意p):", sum(gse_c1@result$NES > 0, na.rm=TRUE), "\n")
# 
# 
# # 诊断 get_density 的实际输出
# test_c1 <- get_density(gse_c1, "C1")
# cat("C1 density 总行数:", nrow(test_c1), "\n")
# cat("C1 density 涉及通路数:", length(unique(test_c1$Description)), "\n")
# print(unique(test_c1$Description))
# 
# # 同时检查 pick_top 选出的3条通路，看它们的成员基因能不能在 geneList 里找到
# sel <- top_all %>% filter(Group == "C1")
# for (i in seq_len(nrow(sel))) {
#   members <- gse_c1@geneSets[[sel$ID[i]]]
#   overlap  <- sum(names(gse_c1@geneList) %in% members)
#   cat(sprintf("通路 %s: geneSets成员=%d, 与geneList重叠=%d\n",
#               sel$ID[i], length(members), overlap))
# }
# 
# # 检查 geneSets 的 key 格式
# cat("result$ID 示例:", head(gse_c1@result$ID, 3), "\n")
# cat("geneSets 名称示例:", head(names(gse_c1@geneSets), 3), "\n")
# 
# ============================================================
# GSEA Ridgeplot（最终版）
# 样式：左侧C123色条 | NES点 | 橙黄密度岭
# ============================================================
# GSEA Ridgeplot（最终版 v4）
# 使用 ggnewscale 双填色：ridge=紫白橙渐变，NES点=蓝白红
# 左侧添加 C1/C2/C3 分组色条
# ============================================================
library(ggridges)
library(patchwork)
library(stringr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggnewscale)

TOP_N  <- 3     # 每组展示通路数
WRAP_W <- 30    # 通路名换行宽度

# ── 1. 选 Top N 通路（NES>0，不限阈值）──────────────────────
pick_top <- function(gse, grp, n = TOP_N) {
  gse@result %>%
    filter(!is.na(p.adjust), !is.na(NES), NES > 0) %>%
    arrange(p.adjust, desc(NES)) %>%
    slice_head(n = n) %>%
    mutate(
      Group     = grp,
      neglog10p = -log10(p.adjust + 1e-10),
      Desc_wrap = str_wrap(Description, WRAP_W)
    ) %>%
    dplyr::select(ID, Description, Desc_wrap, p.adjust, neglog10p, NES, Group)
}

top_all <- bind_rows(
  pick_top(gse_c1, "C1"),
  pick_top(gse_c2, "C2"),
  pick_top(gse_c3, "C3")
)
cat("各组通路数:\n"); print(table(top_all$Group))

# ── 2. 提取基因密度数据 ────────────────────────────────────
get_density <- function(gse, grp) {
  sel <- top_all %>% filter(Group == grp)
  gl  <- gse@geneList
  lapply(seq_len(nrow(sel)), function(i) {
    members <- gse@geneSets[[sel$ID[i]]]
    if (is.null(members)) return(NULL)
    vals <- gl[names(gl) %in% members]
    if (!length(vals)) return(NULL)
    data.frame(Description = sel$Description[i], Group = grp,
               x = as.numeric(vals))
  }) %>% bind_rows()
}

density_df <- bind_rows(
  get_density(gse_c1, "C1"),
  get_density(gse_c2, "C2"),
  get_density(gse_c3, "C3")
)

# ── 3. Y 轴因子顺序（C1顶 → C3底）─────────────────────────
path_order <- top_all %>%
  arrange(factor(Group, levels = c("C3","C2","C1")), p.adjust) %>%
  pull(Desc_wrap) %>% unique()

# 用 join 保证 Desc_wrap 完全一致（避免二次 str_wrap 细微差异）
density_df <- density_df %>%
  left_join(
    top_all %>%
      dplyr::select(Description, Desc_wrap, neglog10p, NES, Group) %>%
      distinct(),
    by = c("Description", "Group")
  ) %>%
  mutate(
    Desc_wrap = factor(Desc_wrap, levels = path_order),
    Group     = factor(Group,     levels = c("C1","C2","C3"))
  ) %>%
  filter(!is.na(Desc_wrap))

cat("密度数据各组通路数:\n")
print(table(density_df$Group, density_df$Desc_wrap))

# ── 4. 计算 NES 点位置和 x 轴范围 ───────────────────────────
# 截断 x 轴（去掉两端极端值，让图更紧凑）
x_lo  <- quantile(density_df$x, 0.005, na.rm = TRUE)
x_hi  <- quantile(density_df$x, 0.998, na.rm = TRUE)
NES_X <- x_lo - (x_hi - x_lo) * 0.18   # NES点固定在密度分布左侧

padj_max <- max(density_df$neglog10p, na.rm = TRUE)
nes_lim  <- max(abs(top_all$NES), na.rm = TRUE)

meta_df <- density_df %>%
  distinct(Desc_wrap, Group, NES, neglog10p) %>%
  mutate(x_nes = NES_X)

# ── 5. 分组颜色 ──────────────────────────────────────────────
grp_col <- c(
  C1 = unname(SUBTYPE_COLORS["C1"]),
  C2 = unname(SUBTYPE_COLORS["C2"]),
  C3 = unname(SUBTYPE_COLORS["C3"])
)

# ── 6. 主图（双填色：ridges + NES点）───────────────────────
p_main <- ggplot(density_df, aes(x = x, y = Desc_wrap)) +
  
  # --- Ridge 密度：紫白橙渐变（颜色差异最大化）---
  geom_density_ridges(
    aes(fill = neglog10p),
    scale = 1.5, alpha = 0.90,
    rel_min_height = 0.01,
    color = "white", linewidth = 0.25
  ) +
  scale_fill_gradient2(
    low      = "#998ec3",       # 紫色（低显著性）
    mid      = "#f7f7f7",       # 白色（中间）
    high     = "#DE712F",       # 橙色（高显著性）
    midpoint = padj_max * 0.4,  # 偏低使橙色更突出
    limits   = c(0, padj_max),
    oob      = squish,
    name     = expression(-log[10](p.adjust))
  ) +
  
  # --- 切换 fill scale ---
  new_scale_fill() +
  
  # --- NES 点：蓝白红渐变 + 大小编码显著性 ---
  geom_point(
    data  = meta_df,
    aes(x = x_nes, y = Desc_wrap, fill = NES, size = neglog10p),
    shape = 21, color = "grey20", stroke = 0.4
  ) +
  scale_fill_gradient2(
    low      = "#4575b4",   # 蓝（NES<0）
    mid      = "#f7f7f7",   # 白（NES≈0）
    high     = "#d73027",   # 红（NES>0）
    midpoint = 0,
    limits   = c(-nes_lim, nes_lim),
    name     = "NES"
  ) +
  scale_size(range = c(3, 9), guide = "none") +
  
  # --- 坐标轴（截短 x 轴）---
  coord_cartesian(xlim = c(NES_X - 0.3, x_hi)) +
  
  labs(
    x     = "Ranked Gene Score (logFC)",
    y     = NULL,
    title = "GSEA Ridgeplot"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.border     = element_rect(color = "#000000", linewidth = 1),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.y      = element_text(size = 8.5, color = "#000000", hjust = 1),
    axis.text.x      = element_text(size = 9,   color = "#000000"),
    axis.title.x     = element_text(size = 11),
    legend.position  = "right",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(4, 6, 4, 0)
  )

# ── 7. 左侧 C1/C2/C3 分组色条 ───────────────────────────────
strip_df <- meta_df %>%
  distinct(Desc_wrap, Group) %>%
  mutate(x = 1)

p_strip <- ggplot(strip_df, aes(x = x, y = Desc_wrap, fill = Group)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = grp_col, name = "Subtype") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin     = margin(4, 2, 4, 4)
  )

# ── 8. 单独提取 Subtype 图例 ─────────────────────────────────
p_leg_src <- ggplot(strip_df, aes(x = x, y = Desc_wrap, fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = grp_col, name = "Subtype") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  )
grp_legend <- cowplot::get_legend(p_leg_src)

# ── 9. 拼图 ────────────────────────────────────────────────
p_combined <- p_strip + p_main +
  plot_layout(widths = c(1, 20), guides = "collect")

p_final <- cowplot::plot_grid(
  p_combined,
  grp_legend,
  nrow       = 1,
  rel_widths = c(1, 0.1)
)

# ── 10. 保存 ────────────────────────────────────────────────
n_paths <- nrow(top_all)
ggsave(
  file.path(enrich_fig_dir, "GSEA_Ridgeplot_final.pdf"),
  p_final,
  width  = 10,
  height = 1.8 + n_paths * 0.72,
  device = "pdf"
)
message("✅ GSEA_Ridgeplot_final.pdf 已保存")



# ==============================================================================
# 5. ssGSEA 通路活性热图 (Fig 3C)
# ==============================================================================
# 精选目标通路
target_go <- list(
  "Mitosis"                  = "GO:0007067",
  "Chromosome Segregation"   = "GO:0007059",
  "ECM Organization"         = "GO:0030198",
  "Ossification"             = "GO:0001503",
  "Adaptive Immune Response" = "GO:0002250",
  "T Cell Activation"        = "GO:0042110"
)

# 提取基因集
get_go_genes <- function(go_id) {
  res <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = go_id,
                          columns = "SYMBOL", keytype = "GO"),
    error = function(e) NULL
  )
  if (!is.null(res)) unique(res$SYMBOL) else NULL
}

pathway_gene_sets <- lapply(target_go, get_go_genes)
pathway_gene_sets <- pathway_gene_sets[!sapply(pathway_gene_sets, is.null)]

# 运行 ssGSEA (兼容新旧版本 GSVA)
expr_mat <- as.matrix(clean_tpm)

if (packageVersion("GSVA") >= "1.50.0") {
  ssgsea_par   <- ssgseaParam(exprData = expr_mat, geneSets = pathway_gene_sets)
  gsva_result  <- gsva(ssgsea_par, verbose = FALSE)
} else {
  gsva_result  <- gsva(expr_mat, pathway_gene_sets,
                       method = "ssgsea", kcdf = "Gaussian")
}

gsva_matrix <- as.matrix(gsva_result)

# 绘制热图
meta_ordered <- meta_step2[order(meta_step2$Subtype), ]
plot_mat_gsva <- gsva_matrix[, meta_ordered$SampleID]

ann_col_gsva <- data.frame(Subtype = meta_ordered$Subtype)
rownames(ann_col_gsva) <- meta_ordered$SampleID

ann_colors_gsva <- list(
  Subtype = SUBTYPE_COLORS[sort(unique(meta_ordered$Subtype))]
)

gap_gsva <- cumsum(table(meta_ordered$Subtype))
gap_gsva <- gap_gsva[-length(gap_gsva)]

pdf(file.path(enrich_fig_dir, "ssGSEA_Pathway_Heatmap.pdf"),
    width = 8, height = 5)
pheatmap(plot_mat_gsva,
         scale        = "row",
         color        = HEATMAP_COLORS,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_col    = ann_col_gsva,
         annotation_colors = ann_colors_gsva,
         gaps_col     = gap_gsva,
         gaps_row     = c(2, 4),
         show_colnames = FALSE,
         main         = "Pathway Activity (ssGSEA)",
         border_color = NA)
dev.off()

# ==============================================================================
# 6. ssGSEA 小提琴图 (Fig 3C 备选)
# ==============================================================================
gsva_df <- as.data.frame(t(gsva_matrix))
gsva_df$SampleID <- rownames(gsva_df)

vio_data <- left_join(gsva_df, meta_step2[, c("SampleID", "Subtype")],
                      by = "SampleID") %>%
  pivot_longer(cols = -c(SampleID, Subtype),
               names_to = "Pathway", values_to = "Score") %>%
  filter(Pathway %in% names(target_go)) %>%
  mutate(Pathway = factor(Pathway, levels = names(target_go)))

my_comparisons <- list(c("C1", "C2"), c("C2", "C3"), c("C1", "C3"))

p_violin <- ggplot(vio_data, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 0.7, alpha = 0.4, colour = "black") +
  scale_fill_manual(values = SUBTYPE_COLORS) +
  stat_compare_means(method = "kruskal.test",
                     label.y.npc = "bottom", size = 3, label = "p.format") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE,
                     tip.length = 0.01, size = 4) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 2) +
  labs(y = "ssGSEA Score", x = NULL) +
  theme_publication(base_size = 12) +
  theme(legend.position = "none")

save_pdf(p_violin,
         file.path(enrich_fig_dir, "ssGSEA_Violin.pdf"),
         width = 8, height = 10)

# ==============================================================================
# 7. 保存
# ==============================================================================
save(go_c1, go_c2, go_c3, kegg_c1, kegg_c2, kegg_c3,
     gse_c1, gse_c2, gse_c3, gsva_matrix,
     file = file.path(RDATA_DIR, "Step4_Enrichment.RData"))

message("🎉 Step 4 完成!")
sessionInfo()
