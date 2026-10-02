# Purpose: External frozen signature and nine-GO score heatmap
# Manuscript/panels: Figure 1h
# Reviewer comments: R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "03_Final_Figure_Code/01_Figure1h.R"
suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(grid);library(GSVA);library(AnnotationDbi);library(org.Hs.eg.db)})
out <- '27_Figure1h_OriginalStyle_20260927'
s <- read.csv('03_Results/R2_3_external_validation/external_signature_scores_recomputed.csv',check.names=FALSE)
stopifnot(nrow(s)==88,!anyDuplicated(s$SampleID))
e<-new.env();load('01_Data/Source_Links/BULK_ENRICH_RDATA',envir=e)
ids<-c('GO:0006260','GO:0098813','GO:0000819','GO:0045229','GO:0043062','GO:0030198','GO:0002443','GO:0002449','GO:0002250')
sets<-list();defs<-list();maps<-list()
for(i in seq_along(ids)){
 g<-get(paste0('gse_c',ceiling(i/3)),e);id<-ids[i];ent<-as.character(g@geneSets[[id]])
 mp<-AnnotationDbi::select(org.Hs.eg.db,keys=ent,keytype='ENTREZID',columns='SYMBOL');mp$GO_ID<-id
 maps[[i]]<-mp;sets[[id]]<-unique(na.omit(mp$SYMBOL))
 rr<-g@result[g@result$ID==id,];defs[[i]]<-data.frame(GO_ID=id,Description=rr$Description,Discovery_group=paste0('C',ceiling(i/3)),Discovery_NES=rr$NES,Discovery_BH=rr$p.adjust)
}
d<-do.call(rbind,defs);write.csv(d,file.path(out,'pathway_definitions.csv'),row.names=FALSE)
write.csv(do.call(rbind,maps),file.path(out,'GO_Entrez_SYMBOL_mapping.csv'),row.names=FALSE)
raw<-matrix(NA_real_,9,88,dimnames=list(ids,s$SampleID));clin<-list();cov<-list()
for(ds in unique(s$Dataset)){
 dr<-file.path('01_Data/Source_Links',paste0('EXTERNAL_',ds,'_DIR'))
 x<-read.csv(file.path(dr,paste0(ds,'_Expression_Log2.csv')),check.names=FALSE)
 m<-as.matrix(x[,-1]);rownames(m)<-as.character(x[[1]]);storage.mode(m)<-'numeric'
 if(anyDuplicated(rownames(m)))m<-limma::avereps(m)
 samples<-s$SampleID[s$Dataset==ds];stopifnot(all(samples%in%colnames(m)))
 m<-m[,samples,drop=FALSE];stopifnot(all(is.finite(m)))
 cov[[ds]]<-data.frame(Dataset=ds,GO_ID=ids,Mapped_genes=lengths(sets),Observed_genes=sapply(sets,function(g)sum(g%in%rownames(m))))
 z<-GSVA::gsva(m,sets,method='ssgsea',ssgsea.norm=TRUE,verbose=FALSE,parallel.sz=1)
 raw[,samples]<-z[ids,samples];clin[[ds]]<-read.csv(file.path(dr,'Clinical_Data.csv'))
}
write.csv(do.call(rbind,cov),file.path(out,'gene_coverage.csv'),row.names=FALSE)
write.csv(data.frame(GO_ID=rownames(raw),raw,check.names=FALSE),file.path(out,'functional_ssGSEA_raw.csv'),row.names=FALSE)
sig<-t(as.matrix(s[,paste0('Signature_C',1:3)]));colnames(sig)<-s$SampleID
rownames(sig)<-c('C1 fixed signature','C2 fixed signature','C3 fixed signature')
mat<-rbind(sig,raw);rownames(mat)[4:12]<-d$Description
for(ds in unique(s$Dataset)){
 ix<-which(s$Dataset==ds);mat[,ix]<-t(scale(t(mat[,ix,drop=FALSE])))
}
stopifnot(all(is.finite(mat)))
group<-ifelse(s$Ambiguous,'Ambiguous',substr(s$Frozen_Assignment,1,2))
s$Display_group<-group
ord<-order(factor(group,levels=c('C1','C2','C3','Ambiguous')),s$Dataset,s$SampleID)
mat<-mat[,ord];s<-s[ord,];cl<-do.call(rbind,clin);cl<-cl[match(s$SampleID,cl$Sample_geo_accession),]
write.csv(data.frame(Feature=rownames(mat),mat,check.names=FALSE),file.path(out,'heatmap_within_cohort_Z.csv'),row.names=FALSE)
write.csv(cbind(s,Sex=cl$Sex,Age=cl$Age,Size=cl$Size),file.path(out,'sample_order_and_annotations.csv'),row.names=FALSE)
cols<-c(C1='#00A087',C2='#4DBBD5',C3='#E64B35',Ambiguous='#BDBDBD')
sex<-as.character(cl$Sex);size<-as.character(cl$Size);nf<-s$NF2_disease_status
mkcols<-function(v,pal){lv<-sort(unique(na.omit(v)));setNames(rep(pal,length.out=length(lv)),lv)}
ha<-HeatmapAnnotation(Study=s$Dataset,Assignment=s$Display_group,Sex=sex,Age=as.numeric(cl$Age),`NF2 disease`=nf,Size=size,
 col=list(Study=c(GSE141801='#A6CEE3',GSE39645='#B2DF8A'),Assignment=cols,Sex=mkcols(sex,c('#CAB2D6','#FFFF99')),
 Age=colorRamp2(c(20,50,80),c('#F7FBFF','#6BAED6','#08306B')),`NF2 disease`=mkcols(nf,c('#FDBF6F','#E5E5E5','#999999')),Size=mkcols(size,c('#FB9A99','#FDDDB5','#C7E9C0'))),
 annotation_name_side='left',na_col='#EEEEEE',simple_anno_size=unit(2,'mm'),annotation_name_gp=gpar(fontsize=6),annotation_legend_param=list(labels_gp=gpar(fontsize=6),title_gp=gpar(fontsize=6)))
make<-function(){Heatmap(mat,name='Within-cohort Z',col=colorRamp2(c(-2,0,2),c('#008837','#F7F7F7','#7B3294')),
 top_annotation=ha,cluster_rows=FALSE,cluster_columns=FALSE,show_column_names=FALSE,
 column_split=factor(s$Display_group,levels=names(cols)),cluster_column_slices=FALSE,column_gap=unit(1,'mm'),
 row_split=factor(c(rep('Fixed signatures',3),rep('C1 functions',3),rep('C2 functions',3),rep('C3 functions',3)),levels=c('Fixed signatures','C1 functions','C2 functions','C3 functions')),
 cluster_row_slices=FALSE,row_gap=unit(1.2,'mm'),row_title=NULL,row_names_side='left',row_names_gp=gpar(fontsize=6),row_names_max_width=unit(58,'mm'),
 column_title_gp=gpar(fontsize=7,fontface='bold'),rect_gp=gpar(col='white',lwd=.15),border=TRUE,
 heatmap_legend_param=list(at=c(-2,0,2),labels=c('≤−2','0','≥2'),title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=6)))}
render<-function(){pushViewport(viewport(gp=gpar(fontfamily='Arial')));draw(make(),newpage=FALSE,merge_legends=TRUE,heatmap_legend_side='right',annotation_legend_side='right',padding=unit(c(3,3,3,3),'mm'));popViewport()}
w<-210/25.4;h<-115/25.4
svglite::svglite(file.path(out,'Fig1h_original_style_functions.svg'),width=w,height=h);render();dev.off()
quartz(type='pdf',file=file.path(out,'Fig1h_original_style_functions.pdf'),width=w,height=h,family='Arial');render();dev.off()
ragg::agg_png(file.path(out,'Fig1h_original_style_functions_preview.png'),width=w,height=h,units='in',res=300);render();dev.off()
ragg::agg_tiff(file.path(out,'Fig1h_original_style_functions.tiff'),width=w,height=h,units='in',res=600,compression='lzw');render();dev.off()
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
