# Purpose: Original myeloid contrast and complete GO BP GSEA from the stored input RDS
# Manuscript/panels: Figure 3g; full supplementary GO table
# Reviewer comments: R2 Comment 5
# Execution requires the inputs and working-directory layout described in README.md; "01_Main_Analyses/08_Myeloid_complete_GO_GSEA.R"
# Reproduce original VS.4/code/06_single cell.R C3-high vs C3-low DEG + GO-BP GSEA.
# Source project supplied separately via VS_SERVER_ROOT
# Original settings: FindMarkers(logfc.threshold=.25, min.pct=.1), then gseGO
# (SYMBOL, BP, pvalueCutoff=.05, minGSSize=10, maxGSSize=500); seed=42.
options(stringsAsFactors=FALSE, scipen=999)
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(Rcpp);library(clusterProfiler);library(org.Hs.eg.db)})
root <- Sys.getenv('VS_SERVER_ROOT', unset='inputs/legacy_project')
setwd(root)
out <- file.path(root,'Reproduction_Fig3g_20261001')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
obj <- readRDS(file.path(root,'RData/Step07_all_myeloid.rds'))
meta <- obj@meta.data
stopifnot(all(c('C3_status') %in% colnames(meta)))
status <- as.character(meta$C3_status)
cells_hi <- rownames(meta)[status=='C3-high']
cells_lo <- rownames(meta)[status=='C3-low']
stopifnot(length(cells_hi)==24668L,length(cells_lo)==51853L)
mat <- GetAssayData(obj,assay='RNA',layer='data')
stopifnot(all(c(cells_hi,cells_lo) %in% colnames(mat)))
# Match Seurat's historical average-log2-fold-change convention used for this project.
mat_hi <- mat[,cells_hi,drop=FALSE]; mat_lo <- mat[,cells_lo,drop=FALSE]
legacy_mean <- function(m) log2(Matrix::rowSums(expm1(m))/ncol(m)+1)
fc <- legacy_mean(mat_hi)-legacy_mean(mat_lo)
pct_hi <- Matrix::rowSums(mat_hi>0)/ncol(mat_hi)
pct_lo <- Matrix::rowSums(mat_lo>0)/ncol(mat_lo)
keep <- pmax(pct_hi,pct_lo)>=0.1 & abs(fc)>=0.25
keep[is.na(keep)] <- FALSE
features <- rownames(mat)[keep]
n1 <- length(cells_hi); n2 <- length(cells_lo); nr <- sum(keep); n <- n1+n2
# Seurat's default Wilcoxon test uses asymptotic rank-sum inference at this sample size.
# Rank in small row blocks to avoid a multi-GB dense matrix; matrixStats uses compiled code.
gene_idx <- which(keep)
pval <- numeric(length(gene_idx))
for (start in seq.int(1L,length(gene_idx),by=32L)) {
  ii <- start:min(start+31L,length(gene_idx))
  xb <- as.matrix(mat[gene_idx[ii],c(cells_hi,cells_lo),drop=FALSE])
  rk <- matrixStats::rowRanks(xb,ties.method='average')
  sr1 <- rowSums(rk[,seq_len(n1),drop=FALSE])
  U <- sr1 - n1*(n1+1)/2
  variance <- n1*n2/(n*(n-1))*rowSums((rk-(n+1)/2)^2)
  delta <- abs(U-n1*n2/2)
  z <- pmax(0,delta-0.5)/sqrt(variance)
  pval[ii] <- ifelse(variance>0,2*pnorm(z,lower.tail=FALSE),1)
  rm(xb,rk)
  if (start %% 512L == 1L) message('Rank-sum rows completed: ',max(ii),' / ',length(gene_idx))
}
deg <- data.frame(gene=features,avg_log2FC=fc[keep],pct.1=pct_hi[keep],pct.2=pct_lo[keep],p_val=pval,p_val_adj=p.adjust(pval,'bonferroni',n=nrow(mat)))
write.csv(deg,file.path(out,'DEG_all_genes_original_method_rerun.csv'),row.names=FALSE)
# Compare significant genes to the source project's saved DEG output, if present.
oldfile<-file.path(root,'Figures/Step07_C3Myeloid/DEG_C3high_vs_C3low.csv')
if(file.exists(oldfile)){old<-read.csv(oldfile);write.csv(data.frame(old_sig_n=nrow(old),rerun_sig_n=sum(deg$p_val_adj<.05 & abs(deg$avg_log2FC)>.5),overlap=sum(old$gene %in% deg$gene[deg$p_val_adj<.05 & abs(deg$avg_log2FC)>.5])),file.path(out,'DEG_cache_comparison.csv'),row.names=FALSE)}
# Original source sets seed globally in config.
set.seed(42)
gene_list<-deg$avg_log2FC; names(gene_list)<-deg$gene; gene_list<-sort(gene_list,decreasing=TRUE); gene_list<-gene_list[is.finite(gene_list)&!is.na(names(gene_list))&nzchar(names(gene_list))];gene_list<-gene_list[!duplicated(names(gene_list))]
run_go<-function(cutoff){set.seed(42);as.data.frame(gseGO(geneList=gene_list,OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType='SYMBOL',ont='BP',pvalueCutoff=cutoff,minGSSize=10,maxGSSize=500,verbose=FALSE))}
orig<-run_go(.05); full<-run_go(1)
write.csv(orig,file.path(out,'GSEA_original_code_pvalueCutoff_0.05.csv'),row.names=FALSE)
write.csv(full,file.path(out,'GSEA_complete_all_GO_BP_terms.csv'),row.names=FALSE)
cat('myeloid cells:',ncol(mat),' C3-high:',n1,' C3-low:',n2,'\n')
cat('genes ranked:',length(gene_list),'\n')
cat('DEGs passing Seurat filters:',nrow(deg),'\n')
cat('GSEA rows returned at original pvalueCutoff=.05:',nrow(orig),'\n')
cat('Complete GO-BP terms at cutoff=1:',nrow(full),'\n')
cat('Complete table FDR<.05:',sum(full$p.adjust<.05,na.rm=TRUE),' positive:',sum(full$NES>0&full$p.adjust<.05,na.rm=TRUE),' negative:',sum(full$NES<0&full$p.adjust<.05,na.rm=TRUE),'\n')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
