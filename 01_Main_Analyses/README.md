# 01_Main_Analyses

| 文件 | 用途 | 正文/图表 | 审稿问题 |
|---|---|---|---|
| [rerun_config.R](_dependencies/rerun_config.R) | Shared input/output configuration for the retained numerical analyses |  |  |
| [discovery_config.R](_dependencies/discovery_config.R) | Discovery thresholds, subtype palette and plotting functions |  |  |
| [01_Consensus_clustering_PCA.R](01_Consensus_clustering_PCA.R) | Primary 5,000-MAD-gene clustering and PCA | Figure 1b; discovery C1/C2/C3 |  |
| [02_Differential_expression.R](02_Differential_expression.R) | Discovery one-versus-rest limma-trend contrasts and 100-gene signatures | Figure 1c; Supplementary Table 2 |  |
| [03_GO_KEGG_GSEA.R](03_GO_KEGG_GSEA.R) | Discovery GO/KEGG ORA and GSEA at the retained primary thresholds | Figure 1d–f |  |
| [04_Official_ESTIMATE.R](04_Official_ESTIMATE.R) | Official ESTIMATE stromal/immune scores and estimated tumour purity | Figure 2a; Supplementary Figure 3b |  |
| [05_Immune_programs_C3_correlation.R](05_Immune_programs_C3_correlation.R) | Recruitment/checkpoint programme scores and C3 purity-adjusted correlation | Figure 2d,e; Supplementary Table 6 | R1 Comment 5.6; R2 Comment 3 |
| [06_External_frozen_signatures.R](06_External_frozen_signatures.R) | Frozen discovery signatures applied to external expression data | Figure 1h; external transcriptional concordance | R2 Comment 4 |
| [07_Singlecell_operational_scores.R](07_Singlecell_operational_scores.R) | Reproduce the manuscript operational module scores and GMM split | Figure 3a–f |  |
| [08_Myeloid_complete_GO_GSEA.R](08_Myeloid_complete_GO_GSEA.R) | Original myeloid contrast and complete GO BP GSEA from the stored input RDS | Figure 3g; full supplementary GO table | R2 Comment 5 |
| [09_MGS_immune_anchoring.R](09_MGS_immune_anchoring.R) | MGS construction with immune/Schwann anchoring independent of subtype labels | Figure 4a–d,f,g |  |
| [10_Clinical_adjusted_models.R](10_Clinical_adjusted_models.R) | Ordinal/Firth models for texture and brainstem adhesion | Figure 5b–e; Supplementary Table 4 | R1 Comment 5.7 |
| [11_MGS_bootstrap_composition.R](11_MGS_bootstrap_composition.R) | 500-replicate MGS bootstrap and official ESTIMATE sequential variance partitioning | Figure 4e,h | R2 Comment 3 |

同一脚本若兼有分析与绘图，只保存一次。`_dependencies` 是所需配置，不是额外分析。
