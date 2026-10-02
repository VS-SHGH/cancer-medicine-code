# Purpose: Top-300 upregulated-gene ORA sensitivity
# Manuscript/panels: Supplementary Table 7 linked data
# Reviewer comments: R2 Comment 4
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/12_External_top300_ORA.R"
# Run after analyse.R. Matches discovery's top-300 upregulated-gene ORA selection.
# Uses an assayed, annotation-matched background rather than all annotated human genes.
out<-'33_External_Functional_Analysis_20260929'
ref<-readRDS(file.path(out,'reference/gene_sets.rds'))
files<-list.files(file.path(out,'results'),pattern='__DE_all.csv$',full.names=TRUE)
stopifnot(length(files)==12)
all<-list()
for(f in files){
 tt<-read.csv(f);lab<-strsplit(sub('__DE_all.csv$','',basename(f)),'__',fixed=TRUE)[[1]]
 # Select on symbol-level results before mapping, matching the historical workflow.
 up<-tt[tt$adj.P.Val<.05 & tt$logFC>0,];up<-up[order(-up$logFC,up$SYMBOL),];up<-head(up,300)
 for(coll in names(ref$sets)){
  universe<-intersect(unique(na.omit(tt$ENTREZID)),unique(unlist(ref$sets[[coll]])))
  selected<-intersect(unique(na.omit(up$ENTREZID)),universe)
  ss<-lapply(ref$sets[[coll]],intersect,y=universe);ss<-ss[lengths(ss)>=10 & lengths(ss)<=500]
  h<-lapply(ss,intersect,y=selected);k<-lengths(h);K<-lengths(ss);N<-length(universe);n<-length(selected)
  p<-if(n)phyper(k-1,K,N-K,n,lower.tail=FALSE)else rep(1,length(k))
  r<-data.frame(Dataset=lab[1],Scope=lab[2],Group=lab[3],Collection=coll,ID=names(ss),Description=unname(ref$descriptions[[coll]][names(ss)]),Count=k,Set_size=K,Selected_n=n,Universe_n=N,P=p,FDR=p.adjust(p,'BH'),Genes=vapply(h,paste,collapse=';',FUN.VALUE=''),row.names=NULL)
  all[[paste(basename(f),coll)]]<-r
 }
}
write.csv(do.call(rbind,all),file.path(out,'results/ORA_top300_up_complete.csv'),row.names=FALSE)
