# Purpose: Discovery thresholds, subtype palette and plotting functions
# Manuscript/panels: shared dependency
# Reviewer comments: not a separate reviewer analysis
# Execution requires the inputs and working-directory layout described in README.md; "01_Main_Analyses/_dependencies/discovery_config.R"
# ==============================================================================
# 全局配置: 听神经瘤 Bulk RNA-seq 分析流程
# ==============================================================================
# 说明: 所有脚本在开头 source 此文件, 统一管理路径/配色/阈值

# --- 项目路径 ---
PROJECT_DIR <- Sys.getenv("VS_LEGACY_WORK_ROOT")
OLD_DIR <- Sys.getenv("VS_SERVER_ROOT")
BULK_CLEAN_DIR <- file.path(Sys.getenv("VS_INPUT_ROOT"), "bulk_cleaned")

# 子目录
RDATA_DIR   <- file.path(PROJECT_DIR, "RData")
FIG_DIR     <- file.path(PROJECT_DIR, "Figures")
TABLE_DIR   <- file.path(PROJECT_DIR, "Tables")

# --- 创建目录 ---
invisible(lapply(
  c(RDATA_DIR, TABLE_DIR,
    file.path(FIG_DIR, "Step1_QC"),
    file.path(FIG_DIR, "Step2_Clustering"),
    file.path(FIG_DIR, "Step3_DEA"),
    file.path(FIG_DIR, "Step4_Enrichment"),
    file.path(FIG_DIR, "Step2_Ext_Validation")),
  function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
))

# # --- 全局配色方案 (统一全文, 基于 ggsci NPG) ---
# SUBTYPE_COLORS <- c(
#   "C1" = "#00A087",   # 红色 — Proliferative
#   "C2" = "#4DBBD5",   # 蓝色 — Mesenchymal
#   "C3" = "#E64B35 "   # 绿色 — Immune
# )
# 
# # 带全称的配色 (用于外部验证等场景)
# SUBTYPE_COLORS_FULL <- c(
#   "C1_Proliferative" = "#00A087",
#   "C2_Mesenchymal"   = "#4DBBD5",
#   "C3_Immune"        = "#E64B35"
# )
# --- 全局配色方案 (统一全文, 基于 ggsci NPG) ---
SUBTYPE_COLORS <- c(
  "C1" = "#00A087",   # 红色 — Proliferative
  "C2" = "#4DBBD5",   # 蓝色 — Mesenchymal
  "C3" =  "#E64B35"   # 绿色 — Immune
)
SUBTYPE_COLORS <- trimws(SUBTYPE_COLORS)  # 防御性去空格

# 带全称的配色 (用于外部验证等场景)
SUBTYPE_COLORS_FULL <- c(
  "C1_Proliferative" = "#00A087",
  "C2_Mesenchymal"   = "#4DBBD5",
  "C3_Immune"        = "#E64B35"
)
SUBTYPE_COLORS_FULL <- trimws(SUBTYPE_COLORS_FULL)

# 热图配色
HEATMAP_COLORS <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 相关性热图配色
COR_COLORS <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(100)

# --- 分析参数 ---
PARAMS <- list(
  # Step 1: QC
  blacklist_samples = c("S42", "S37", "S27", "S38", "S26", "S20"),
  min_tpm           = 1,       # TPM 表达阈值
  min_sample_frac   = 0.2,     # 至少在 20% 样本中表达
  
  # Step 2: Clustering
  n_top_mad_genes   = 5000,    # MAD 筛选基因数
  ccp_maxK          = 6,
  ccp_reps          = 1000,
  ccp_pItem         = 0.8,
  ccp_seed          = 123456,
  k_final           = 3,       # 最终选择的 k 值 (需 CDF/PAC 支持)
  
  # Step 3: DEA
  dea_fdr_cutoff    = 0.05,
  dea_logfc_cutoff  = 1,
  dea_top_n_markers = 20,
  
  # Step 4: Enrichment
  enrich_top_n      = 300,     # Top N 基因用于 ORA
  gsea_pval_cutoff  = 0.05,
  
  # 随机种子
  global_seed       = 42
)

# --- 顶刊绘图主题 ---
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # 面板
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "black", linewidth = 0.8),
      # 坐标轴
      axis.text         = element_text(colour = "black", size = base_size * 0.9),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      # 图例
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.title      = element_text(face = "bold", size = base_size * 0.9),
      legend.text       = element_text(size = base_size * 0.85),
      # 标题
      plot.title        = element_text(face = "bold", hjust = 0.5, size = base_size * 1.1),
      plot.subtitle     = element_text(hjust = 0.5, size = base_size * 0.9),
      # 分面
      strip.background  = element_rect(fill = "grey95", colour = "black"),
      strip.text        = element_text(face = "bold", size = base_size * 0.9)
    )
}

# --- 半矢量图保存函数 (点多时使用) ---
save_half_vector <- function(plot, filename, width = 7, height = 6, dpi = 300) {
  # 使用 cairo_pdf 设备, 对大量散点自动栅格化
  cairo_pdf(filename, width = width, height = height)
  print(plot)
  dev.off()
  message("✅ 已保存: ", filename)
}

# --- 标准 PDF 保存函数 ---
save_pdf <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename, plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  message("✅ 已保存: ", filename)
}

message("✅ 全局配置加载完毕 | 项目路径: ", PROJECT_DIR)
