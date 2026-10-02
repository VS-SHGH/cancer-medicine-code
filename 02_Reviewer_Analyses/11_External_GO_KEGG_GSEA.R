# Purpose: External GO/KEGG ORA and preranked GSEA with frozen assignments
# Manuscript/panels: Supplementary Figure 7; Supplementary Table 7 and linked data
# Reviewer comments: R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/11_External_GO_KEGG_GSEA.R"
# Run from the active Cancer Medicine workspace with Rscript.
# Frozen external assignments; no manuscript, historical input or existing result is modified.
suppressPackageStartupMessages({library(limma);library(fgsea);library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db);library(GSVA)})
set.seed(20260929)
out <- '33_External_Functional_Analysis_20260929'
wr <- function(x,n) write.csv(x,file.path(out,'results',paste0(n,'.csv')),row.names=FALSE)
ann <- read.csv('27_Figure1h_OriginalStyle_20260927/sample_order_and_annotations.csv',check.names=FALSE)
stopifnot(nrow(ann)==88,!anyDuplicated(ann$SampleID))
e<-new.env();load('01_Data/Source_Links/BULK_DEA_RDATA',e)
classifier<-unique(unlist(e$gene_sets));writeLines(classifier,file.path(out,'reference/classifier_symbols.txt'))
# GO BP annotations are frozen by saving this local database-derived collection.
mp<-AnnotationDbi::select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype='ENTREZID'),keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL'))
mp<-unique(mp[!is.na(mp$GOALL)&mp$ONTOLOGYALL=='BP',c('GOALL','ENTREZID')])
go<-lapply(split(mp$ENTREZID,mp$GOALL),unique)
gd<-AnnotationDbi::select(GO.db,keys=names(go),keytype='GOID',columns='TERM')
go_names<-setNames(gd$TERM,gd$GOID)
kg<-read.delim(file.path(out,'reference/KEGG_pathway_genes.tsv'),header=FALSE)
kn<-read.delim(file.path(out,'reference/KEGG_pathway_names.tsv'),header=FALSE)
kg[[1]]<-sub('^path:','',kg[[1]]);kg[[2]]<-sub('^hsa:','',kg[[2]])
kegg<-lapply(split(kg[[2]],kg[[1]]),unique);kn[[1]]<-sub('^path:','',kn[[1]])
kn[[2]]<-sub(' - Homo sapiens \\(human\\)$','',kn[[2]])
collections<-list(GO_BP=go,KEGG=kegg);descriptions<-list(GO_BP=go_names,KEGG=setNames(kn[[2]],kn[[1]]))
saveRDS(list(sets=collections,descriptions=descriptions),file.path(out,'reference/gene_sets.rds'))
sm<-AnnotationDbi::select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype='SYMBOL'),keytype='SYMBOL',columns='ENTREZID')
sm<-sm[!is.na(sm$ENTREZID),];sm<-sm[order(sm$SYMBOL,as.numeric(sm$ENTREZID)),]
# Deterministic single Entrez ID per symbol; duplicate Entrez IDs later retain smallest gene-level P.
sm<-sm[!duplicated(sm$SYMBOL),];wr(sm,'symbol_entrez_mapping')
readmat<-function(ds){x<-read.csv(file.path('01_Data/Source_Links',paste0('EXTERNAL_',ds,'_DIR'),paste0(ds,'_Expression_Log2.csv')),check.names=FALSE);m<-as.matrix(x[,-1]);rownames(m)<-x[[1]];storage.mode(m)<-'numeric';if(anyDuplicated(rownames(m)))m<-avereps(m);stopifnot(all(is.finite(m)));m}
ora<-function(selected,universe,sets){
 universe<-intersect(universe,unique(unlist(sets)));selected<-intersect(selected,universe)
 ss<-lapply(sets,intersect,y=universe);ss<-ss[lengths(ss)>=10 & lengths(ss)<=500]
 h<-lapply(ss,intersect,y=selected);k<-lengths(h);K<-lengths(ss);N<-length(universe);n<-length(selected)
 p<-if(n>0)phyper(k-1,K,N-K,n,lower.tail=FALSE) else rep(1,length(k))
 data.frame(ID=names(ss),Count=k,Set_size=K,Selected_n=n,Universe_n=N,P=p,FDR=p.adjust(p,'BH'),Genes=vapply(h,paste,collapse=';',FUN.VALUE=''),row.names=NULL)
}
summ<-list();qs<-list();gseout<-list();oraout<-list()
for(ds in unique(ann$Dataset)){
 m0<-readmat(ds)
 for(scope in c('all_unambiguous','nonirradiated_sporadic')){
 a<-ann[ann$Dataset==ds & !ann$Ambiguous,]
 if(scope=='nonirradiated_sporadic')a<-a[a$Primary_Nonirradiated_Sporadic,]
 stopifnot(all(a$SampleID%in%colnames(m0)))
 m<-m0[,a$SampleID,drop=FALSE];m<-m[apply(m,1,sd)>0,,drop=FALSE]
 g<-factor(a$Display_group,levels=c('C1','C2','C3'));design<-model.matrix(~0+g);colnames(design)<-levels(g)
 stopifnot(all(table(g)>=3),qr(design)$rank==3)
 fit<-eBayes(contrasts.fit(lmFit(m,design),makeContrasts(C1=C1-(C2+C3)/2,C2=C2-(C1+C3)/2,C3=C3-(C1+C2)/2,levels=design)),trend=TRUE)
 for(j in 1:3){
 label<-paste(ds,scope,paste0('C',j),sep='__');message(label)
 tt<-topTable(fit,coef=j,number=Inf,sort.by='none');tt$SYMBOL<-rownames(tt);tt$ENTREZID<-sm$ENTREZID[match(tt$SYMBOL,sm$SYMBOL)]
 wr(tt,paste0(label,'__DE_all'))
 tt<-tt[!is.na(tt$ENTREZID),];tt<-tt[order(tt$P.Value,tt$SYMBOL),];tt<-tt[!duplicated(tt$ENTREZID),]
 qs[[label]]<-data.frame(Dataset=ds,Scope=scope,Group=paste0('C',j),N_C1=sum(g=='C1'),N_C2=sum(g=='C2'),N_C3=sum(g=='C3'),Tested_symbols=nrow(m),Mapped_unique_Entrez=nrow(tt),Up=sum(tt$adj.P.Val<.05&tt$logFC>0),Down=sum(tt$adj.P.Val<.05&tt$logFC<0))
 for(coll in names(collections)){
  for(direction in c('Up','Down')){
   sel<-tt$ENTREZID[tt$adj.P.Val<.05 & if(direction=='Up')tt$logFC>0 else tt$logFC<0]
   r<-ora(sel,tt$ENTREZID,collections[[coll]]);r$Description<-unname(descriptions[[coll]][r$ID]);r$Direction<-direction
   r$Dataset<-ds;r$Scope<-scope;r$Group<-paste0('C',j);r$Collection<-coll
   oraout[[paste(label,coll,direction)]]<-r
  }
  # Match discovery GSEA ranking (logFC); no DEG threshold is used for GSEA.
  ranks<-sort(setNames(tt$logFC,tt$ENTREZID),decreasing=TRUE)
  set.seed(20260929+j)
  rr<-as.data.frame(fgseaMultilevel(collections[[coll]],ranks,minSize=10,maxSize=500,eps=1e-8,nPermSimple=1000,nproc=1))
  rr$leadingEdge<-vapply(rr$leadingEdge,paste,collapse=';',FUN.VALUE='')
  rr$Description<-unname(descriptions[[coll]][rr$pathway]);rr$Dataset<-ds;rr$Scope<-scope;rr$Group<-paste0('C',j);rr$Collection<-coll
  gseout[[paste(label,coll)]]<-rr
  wr(rr,paste0(label,'__',coll,'_GSEA'))
 }
 wr(do.call(rbind,oraout[grepl(label,names(oraout),fixed=TRUE)]),paste0(label,'__ORA'))
 }
}
}
wr(do.call(rbind,qs),'analysis_sample_and_DE_counts')
wr(do.call(rbind,gseout),'GSEA_all_complete')
wr(do.call(rbind,oraout),'ORA_all_complete')
# Direct like-for-like comparison: the same nine GO ssGSEA scores in discovery and external cohorts.
defs<-read.csv('27_Figure1h_OriginalStyle_20260927/pathway_definitions.csv')
map9<-read.csv('27_Figure1h_OriginalStyle_20260927/GO_Entrez_SYMBOL_mapping.csv')
sets9<-lapply(split(map9$SYMBOL,map9$GO_ID),function(x)unique(na.omit(x)))
e<-new.env();load('01_Data/Source_Links/BULK_STEP1_RDATA',e)
meta<-read.csv('03_Results/R2_2_composition/composition_scores_by_sample.csv')
mat<-e$clean_tpm[,meta$SampleID,drop=FALSE]
z<-gsva(log2(mat+1),sets9,method='ssgsea',ssgsea.norm=TRUE,verbose=FALSE,parallel.sz=1)
wr(data.frame(GO_ID=rownames(z),z,check.names=FALSE),'discovery_nine_GO_ssGSEA_raw')
raw<-read.csv('27_Figure1h_OriginalStyle_20260927/functional_ssGSEA_raw.csv',check.names=FALSE);ext<-as.matrix(raw[,-1]);rownames(ext)<-raw[[1]]
tabs<-list()
for(ds in c('Discovery',unique(ann$Dataset))){
 if(ds=='Discovery'){v<-z;a<-data.frame(SampleID=meta$SampleID,Display_group=meta$Subtype)}else{a<-ann[ann$Dataset==ds,];v<-ext[,a$SampleID,drop=FALSE]}
 zz<-t(scale(t(v)))
 for(id in rownames(v))for(gr in c('C1','C2','C3')){
  ids<-a$SampleID[a$Display_group==gr];tabs[[paste(ds,id,gr)]]<-data.frame(Dataset=ds,GO_ID=id,Description=defs$Description[match(id,defs$GO_ID)],Group=gr,N=length(ids),Median_raw=median(v[id,ids]),Median_within_cohort_Z=median(zz[id,ids]))
 }
}
wr(do.call(rbind,tabs),'nine_GO_discovery_external_medians')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('COMPLETE')
writeLines(capture.output(warnings()),file.path(out,'warnings.txt'))
