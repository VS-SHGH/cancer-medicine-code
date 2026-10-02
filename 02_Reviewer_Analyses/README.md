# 02_Reviewer_Analyses

| 文件 | 用途 | 正文/图表 | 审稿问题 |
|---|---|---|---|
| [revision_paths.R](_dependencies/revision_paths.R) | Named input slots for reviewer analyses |  |  |
| [singlecell_config.R](_dependencies/singlecell_config.R) | Configuration used by the paired single-cell sensitivity analyses |  |  |
| [01_Clustering_parameter_sensitivity.R](01_Clustering_parameter_sensitivity.R) | 36-setting clustering sensitivity and consensus margins | Supplementary Figure 1 | R2 Comment 2 |
| [02_Clustering_2000_3000_5000.R](02_Clustering_2000_3000_5000.R) | Focused 2,000/3,000/5,000-gene consensus comparison | Supplementary Figure 1 | R2 Comment 2 |
| [03_Composition_MGS_sensitivity.R](03_Composition_MGS_sensitivity.R) | Composition-associated variation and residualized MGS | Supplementary Figure 2a | R2 Comment 3 |
| [04_Composition_residual_clustering.R](04_Composition_residual_clustering.R) | Clustering after removing composition-associated expression | Supplementary Figure 2b | R2 Comment 3 |
| [05_External_assignment_sensitivity.R](05_External_assignment_sensitivity.R) | Signature coverage, ambiguous assignments and external strata | Supplementary Figure 2c | R2 Comment 4 |
| [06_External_MGS_sample_alignment.R](06_External_MGS_sample_alignment.R) | External MGS matching by sample identity | Supplementary Figure 2d | R1 Comment 5.2; R2 Comment 4 |
| [07_Deconvolution_and_clinical_sensitivity.R](07_Deconvolution_and_clinical_sensitivity.R) | CIBERSORT-pass comparisons, cross-method correlations and adjusted clinical analyses | Supplementary Figure 3a; Supplementary Tables 3,4 | R1 Comments 5.3–5.7 |
| [08_Flow_percentage_zscores.R](08_Flow_percentage_zscores.R) | 26 flow variables: percentage units and mean patient-level z scores | Supplementary Figure 5b | R2 Comment 6 |
| [09_Complete_blood_screen.R](09_Complete_blood_screen.R) | 83 group tests, 249 Dunn comparisons and 168 operative-phenotype correlations | Supplementary Figure 5c,d; Figure 5g–i interpretation | R1 Comment 5.6; R2 Comment 7 |
| [10_Score_purity_sensitivity.R](10_Score_purity_sensitivity.R) | Four score–purity correlations and partial-correlation sensitivity | Supplementary Table 6 | R1 Comment 5.6; R2 Comment 3 |
| [11_External_GO_KEGG_GSEA.R](11_External_GO_KEGG_GSEA.R) | External GO/KEGG ORA and preranked GSEA with frozen assignments | Supplementary Figure 7; Supplementary Table 7 and linked data | R2 Comment 4 |
| [12_External_top300_ORA.R](12_External_top300_ORA.R) | Top-300 upregulated-gene ORA sensitivity | Supplementary Table 7 linked data | R2 Comment 4 |
| [13_External_functional_summary.py](13_External_functional_summary.py) | Retained external functional result summaries | Supplementary Table 7 linked data | R2 Comment 4 |
| [14_Normalized_singlecell_sensitivity.R](14_Normalized_singlecell_sensitivity.R) | Normalized module scores and patient-level pseudobulk preparation | Single-cell sensitivity results | R2 Comment 5 |
| [15_Patient_paired_CAMERA.R](15_Patient_paired_CAMERA.R) | Within-patient paired CAMERA and depth-adjusted pathway sensitivity | Supplementary Figure 3d and linked data | R2 Comment 5 |
| [16_CellChat_normalization_sensitivity.R](16_CellChat_normalization_sensitivity.R) | Fresh normalized CellChat sensitivity and patient-level ligand/receptor support | CellChat supplementary sensitivity | R1 Comment 5.8; R2 Comment 5 |

同一脚本若兼有分析与绘图，只保存一次。`_dependencies` 是所需配置，不是额外分析。
