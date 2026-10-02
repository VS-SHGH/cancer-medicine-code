# Purpose: Focused 2,000/3,000/5,000-gene consensus comparison
# Manuscript/panels: Supplementary Figure 1
# Reviewer comments: R2 Comment 2
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/02_Clustering_2000_3000_5000.R"
#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(ConsensusClusterPlus); library(cluster); library(ggplot2); library(pheatmap)})
script_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script_path <- normalizePath(sub('^--file=','',script_arg[1]))
source(file.path(dirname(dirname(script_path)), '00_paths.R'))
OUT <- file.path(RESULT_DIR,'R2_1_clustering','focused_2000_3000_5000')
FIG <- file.path(FIGURE_DIR,'R2_1_clustering','focused_2000_3000_5000')
LOG <- file.path(LOG_DIR,'R2_1_clustering','focused_2000_3000_5000')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE); dir.create(FIG,recursive=TRUE,showWarnings=FALSE); dir.create(LOG,recursive=TRUE,showWarnings=FALSE)

e <- new.env(parent=emptyenv()); load(BULK_STEP1_RDATA,envir=e)
e2 <- new.env(parent=emptyenv()); load(BULK_CLUSTER_RDATA,envir=e2)
tpm <- as.matrix(e$clean_tpm); ref <- setNames(as.character(e2$meta_step2$Subtype),e2$meta_step2$SampleID); tpm <- tpm[,names(ref),drop=FALSE]
logt <- log2(tpm+1); madv <- apply(logt,1,mad,na.rm=TRUE)

adjusted_rand <- function(x,y){tab<-table(x,y);c2<-function(z)z*(z-1)/2;n<-sum(tab);a<-sum(c2(tab));b<-sum(c2(rowSums(tab)));c<-sum(c2(colSums(tab)));ex<-b*c/c2(n);(a-ex)/((b+c)/2-ex)}
consensus_metrics <- function(cc,x,k,g){cm<-cc[[k]]$consensusMatrix;v<-cm[lower.tri(cm)];si<-silhouette(as.integer(cc[[k]]$consensusClass),as.dist(1-cor(x)));data.frame(genes=g,k=k,PAC=mean(v>.1&v<.9),CDF_area=1-mean(v),mean_silhouette=mean(si[,'sil_width']),min_silhouette=min(si[,'sil_width']),cluster_sizes=paste(as.integer(table(cc[[k]]$consensusClass)),collapse='/'))}
map_labels <- function(lab){ids<-names(lab);r<-as.integer(sub('^C','',ref[ids]));p<-as.matrix(expand.grid(rep(list(1:3),3)));p<-p[apply(p,1,function(z)length(unique(z))==3),,drop=FALSE];a<-apply(p,1,function(z)mean(z[as.integer(lab)]==r));setNames(paste0('C',p[which.max(a),as.integer(lab)]),ids)}

allm<-list(); all_lab<-list(); ccs<-list()
for (g in c(2000L,3000L,5000L)) {
  genes<-names(sort(madv,decreasing=TRUE))[seq_len(g)]; x<-logt[genes,,drop=FALSE]; x<-sweep(x,1,apply(x,1,median),'-')
  set.seed(123456)
  cc<-ConsensusClusterPlus(x,maxK=6,reps=1000,pItem=.8,pFeature=1,clusterAlg='hc',distance='pearson',innerLinkage='average',finalLinkage='average',seed=123456,title=file.path(LOG,paste0('G',g)),plot=NULL,writeTable=FALSE,verbose=FALSE)
  ccs[[as.character(g)]]<-cc
  allm[[as.character(g)]]<-do.call(rbind,lapply(2:6,function(k)consensus_metrics(cc,x,k,g)))
  mapped<-map_labels(cc[[3]]$consensusClass)
  all_lab[[as.character(g)]]<-data.frame(genes=g,SampleID=names(mapped),reference_Subtype=ref[names(mapped)],focused_Subtype=mapped,exact_match=mapped==ref[names(mapped)])
}
metrics<-do.call(rbind,allm); metrics$delta_CDF_area<-ave(metrics$CDF_area,metrics$genes,FUN=function(z)c(NA,diff(z)))
labs<-do.call(rbind,all_lab)
ari<-do.call(rbind,lapply(split(labs,labs$genes),function(z)data.frame(genes=z$genes[1],ARI_vs_primary=adjusted_rand(z$focused_Subtype,z$reference_Subtype),exact_agreement=mean(z$exact_match),changed_samples=paste(z$SampleID[!z$exact_match],collapse=';'))))
write.csv(metrics,file.path(OUT,'focused_2000_3000_5000_k2_to_k6_metrics.csv'),row.names=FALSE)
write.csv(labs,file.path(OUT,'focused_2000_3000_5000_k3_patient_labels.csv'),row.names=FALSE)
write.csv(ari,file.path(OUT,'focused_2000_3000_5000_k3_concordance.csv'),row.names=FALSE)

subtype_cols<-c(C1='#B05A7A',C2='#D99B3D',C3='#4C78A8')
for(g in c(2000L,3000L,5000L)){
  cc<-ccs[[as.character(g)]]; cm<-cc[[3]]$consensusMatrix; ids<-names(cc[[3]]$consensusClass); if(is.null(rownames(cm)))dimnames(cm)<-list(ids,ids); ord<-order(ref[ids],ids); ann<-data.frame(Subtype=ref[ids][ord]);rownames(ann)<-ids[ord]; gaps<-cumsum(as.integer(table(ref[ids][ord])));gaps<-gaps[-length(gaps)]
  f<-file.path(FIG,paste0('consensus_matrix_k3_',g,'_MAD'))
  pdf(paste0(f,'.pdf'),width=7.2,height=6.8);pheatmap(cm[ids[ord],ids[ord]],cluster_rows=FALSE,cluster_cols=FALSE,annotation_col=ann,annotation_colors=list(Subtype=subtype_cols),gaps_row=gaps,gaps_col=gaps,show_rownames=FALSE,show_colnames=FALSE,border_color=NA,color=colorRampPalette(c('white','#CFE8F3','#2166AC'))(100),main=paste0('Consensus matrix (k=3; ',format(g,big.mark=','),' MAD genes)'));dev.off()
}

p1<-ggplot(metrics,aes(k,mean_silhouette,group=factor(genes),colour=factor(genes)))+geom_line(linewidth=.8)+geom_point(size=2)+scale_x_continuous(breaks=2:6)+scale_colour_manual(values=c('2000'='#8C6BB1','3000'='#E08214','5000'='#2C7FB8'),name='MAD genes')+labs(x='Number of clusters (k)',y='Mean silhouette width')+theme_classic(base_size=9,base_family='sans')+theme(legend.position='top')
p2<-ggplot(metrics,aes(k,PAC,group=factor(genes),colour=factor(genes)))+geom_line(linewidth=.8)+geom_point(size=2)+scale_x_continuous(breaks=2:6)+scale_colour_manual(values=c('2000'='#8C6BB1','3000'='#E08214','5000'='#2C7FB8'),name='MAD genes')+labs(x='Number of clusters (k)',y='PAC (0.1–0.9)')+theme_classic(base_size=9,base_family='sans')+theme(legend.position='top')
savep<-function(p,n,w=5.8,h=4.2){ggsave(file.path(FIG,paste0(n,'.pdf')),p,width=w,height=h,device=pdf);ggsave(file.path(FIG,paste0(n,'.svg')),p,width=w,height=h,device=svglite::svglite);ggsave(file.path(FIG,paste0(n,'.tiff')),p,width=w,height=h,device=ragg::agg_tiff,dpi=600,compression='lzw');ggsave(file.path(FIG,paste0(n,'_preview.png')),p,width=w,height=h,device=ragg::agg_png,dpi=300)}
savep(p1,'focused_mean_silhouette_by_k');savep(p2,'focused_PAC_by_k')

sel<-subset(metrics,k==3);write.csv(sel,file.path(OUT,'focused_k3_comparison.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(LOG,'sessionInfo.txt'))
writeLines(c('Focused sensitivity analysis: top 2,000, 3,000, and 5,000 MAD genes.','Primary transformation: log2(TPM+1), per-gene median centring.','ConsensusClusterPlus: Pearson distance, average linkage, 80% item subsampling, 1,000 resamples, k=2–6.','k=3 rationale: highest mean silhouette in the candidate range, balanced cluster sizes, and no two-sample cluster; PAC alone was not used as the selection rule.'),file.path(LOG,'analysis_note.txt'))
print(sel);print(ari)
