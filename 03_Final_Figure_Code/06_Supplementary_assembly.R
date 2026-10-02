# Purpose: Final multi-panel assembly; extracted required helpers only
# Manuscript/panels: Supplementary Figures 1–3,5,6
# Reviewer comments: not a separate reviewer analysis
# Execution requires the inputs and working-directory layout described in README.md; "03_Final_Figure_Code/06_Supplementary_assembly.R"
# Reassemble archived quantitative outputs; no statistical models are rerun.
# All rows are retained. Colours and each image aspect ratio are preserved.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr);library(grid);library(png);library(ragg);library(svglite)})
old <- '09_Revision_Package/Figures'
out <- '35_Supplement_Reassembly_20260930/Figures';dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Reuse native plotting objects without running old export commands.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(readr); library(grid); library(png); library(ragg); library(svglite)})
root <- normalizePath('.')
package_dir <- file.path(root,'09_Revision_Package')
figure_dir <- file.path(package_dir,'Figures')
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
baseline <- file.path(dirname(root),'01_Baseline_Cancer_Medicine_Submission')
replacement <- file.path(root,'04_Figures/Replacement_Panels_Original_Style')
group_colors <- c(C1='#00A087',C2='#4DBBD5',C3='#E64B35')
theme_set(theme_classic(base_family='Arial',base_size=8)+theme(axis.text=element_text(colour='black'),legend.position='bottom',plot.title=element_text(size=9,face='bold')))
read_table <- function(path) read_csv(path,show_col_types=FALSE)
crop_image <- function(number,left,top,right,bottom) {
  pixels <- readPNG(file.path(baseline,'02_Figures',paste0('Figure',number,'_00.png')))
  height <- dim(pixels)[1]; width <- dim(pixels)[2]
  pixels[seq(max(1,round(top*height)),min(height,round(bottom*height))),seq(max(1,round(left*width)),min(width,round(right*width))),,drop=FALSE]
}
image_file <- function(path) readPNG(path)
place <- function(object,x,y,width,height,label=NULL) {
  pushViewport(viewport(x=x,y=1-y,width=width,height=height,just=c('left','top')))
  if (inherits(object,'ggplot')) print(object,vp=viewport(width=.96,height=.94)) else {
    if(is.character(object)) object <- image_file(object)
    ratio <- dim(object)[2]/dim(object)[1]
    panel_ratio <- convertWidth(unit(1,'npc'),'mm',valueOnly=TRUE)/convertHeight(unit(1,'npc'),'mm',valueOnly=TRUE)
    if(ratio>panel_ratio) grid.raster(object,width=unit(.96,'npc'),height=unit(.96*panel_ratio/ratio,'npc')) else grid.raster(object,height=unit(.94,'npc'),width=unit(.94*ratio/panel_ratio,'npc'))
  }
  if(!is.null(label)) grid.text(label,x=unit(.2,'mm'),y=unit(1,'npc')-unit(.8,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=10,fontface='bold'))
  popViewport()
}
place_bottom_global <- function(object,x,y,width,height,label=NULL) {
  pushViewport(viewport(x=x,y=y,width=width,height=height,just=c('left','bottom')))
  if (inherits(object,'ggplot')) print(object,vp=viewport(width=.97,height=.94)) else {
    if (is.character(object)) object <- image_file(object)
    ratio <- dim(object)[2]/dim(object)[1]
    panel_ratio <- convertWidth(unit(1,'npc'),'mm',valueOnly=TRUE)/convertHeight(unit(1,'npc'),'mm',valueOnly=TRUE)
    if (ratio > panel_ratio) grid.raster(object,width=unit(.97,'npc'),height=unit(.97*panel_ratio/ratio,'npc'))
    else grid.raster(object,height=unit(.94,'npc'),width=unit(.94*ratio/panel_ratio,'npc'))
  }
  if(!is.null(label)) grid.text(label,x=unit(.2,'mm'),y=unit(1,'npc')-unit(.8,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=10,fontface='bold'))
  popViewport()
}
export <- function(stem,draw_fun,width_mm=183,height_mm=225) {
  width_in <- width_mm/25.4; height_in <- height_mm/25.4
  svglite(file.path(figure_dir,paste0(stem,'.svg')),width=width_in,height=height_in); grid.newpage(); draw_fun(); dev.off()
  quartz(type='pdf',file=file.path(figure_dir,paste0(stem,'.pdf')),width=width_in,height=height_in,family='Arial'); grid.newpage(); draw_fun(); dev.off()
  agg_tiff(file.path(figure_dir,paste0(stem,'.tiff')),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw',background='white'); grid.newpage(); draw_fun(); dev.off()
  agg_png(file.path(figure_dir,paste0(stem,'_preview.png')),width=width_mm,height=height_mm,units='mm',res=220,background='white'); grid.newpage(); draw_fun(); dev.off()
}
panel <- function(stem) file.path(replacement,paste0(stem,'_preview.png'))

# Reuse the documented graphics helpers without rerunning all main-figure exports.
margins<-read_table('03_Results/R2_1_clustering/baseline_k3_patient_consensus_stability.csv') %>% arrange(Subtype,consensus_margin) %>% mutate(SampleID=factor(SampleID,levels=SampleID))
margin_plot<-ggplot(margins,aes(consensus_margin,SampleID,fill=Subtype))+geom_col(width=.8)+scale_fill_manual(values=group_colors)+labs(x='Within-group consensus margin',y='Patient',fill='Group')+theme(axis.text.y=element_text(size=5),legend.key.size=unit(2,'mm'))
labels<-read_table('03_Results/R2_2_composition/residualized_k3_patient_labels.csv') %>% count(original_Subtype,residualized_Subtype) %>% complete(original_Subtype,residualized_Subtype,fill=list(n=0))
label_plot<-ggplot(labels,aes(residualized_Subtype,original_Subtype,fill=n))+geom_tile(colour='white')+geom_text(aes(label=n),size=3)+scale_fill_gradient(low='white',high='#2166AC',name='Patients')+labs(x='Residualized group (mapped label)',y='Original group',title='Composition-adjusted label agreement',subtitle='ARI = 0.078; label agreement = 50%')+theme(legend.position='right',plot.title=element_text(size=7),plot.subtitle=element_text(size=6))
counts<-read_table('03_Results/R2_3_external_validation/external_assignment_summary.csv') %>% mutate(Assignment=sub('_.*','',Assignment))
count_plot<-ggplot(counts,aes(Stratum,n,fill=Assignment))+geom_col(width=.7)+facet_wrap(~Dataset)+scale_fill_manual(values=c(group_colors,Ambiguous='grey60'))+labs(x=NULL,y='Samples',fill=NULL)+theme(axis.text.x=element_text(angle=20,hjust=1,size=6),strip.background=element_blank(),strip.text=element_text(size=7),legend.key.size=unit(2,'mm'),legend.text=element_text(size=6))
composition<-read_table('03_Results/R2_2_composition/composition_scores_by_sample.csv') %>% pivot_longer(c(Stromal_Score,Immune_Score,ESTIMATE_Score,Tumour_Purity),names_to='Score',values_to='Value')
cp<-ggplot(composition,aes(Subtype,Value,fill=Subtype))+geom_boxplot(outlier.shape=NA,linewidth=.25)+geom_point(position=position_jitter(width=.08,height=0,seed=123),size=.6)+facet_wrap(~Score,scales='free_y')+scale_fill_manual(values=group_colors,guide='none')+labs(x=NULL,y=NULL)+theme(strip.background=element_blank(),strip.text=element_text(size=6))
flow<-read_table('03_Results/R2_5_flow_methods/R2_5_flow_zscore_group_matrix.csv') %>% mutate(label=factor(label,levels=rev(label))) %>% pivot_longer(ends_with('_z_mean'),names_to='Group',values_to='Z') %>% mutate(Group=sub('_z_mean','',Group))
fp<-ggplot(flow,aes(Group,label,fill=Z))+geom_tile(colour='white',linewidth=.15)+scale_fill_gradient2(low='#809CC0',mid='white',high='#D87570',name='Mean z')+labs(x=NULL,y=NULL)+theme(axis.text.y=element_text(size=5.8),legend.key.height=unit(2,'mm'),legend.key.width=unit(5,'mm'))
trim<-function(im){ mask<-pmin(im[,,1],im[,,2],im[,,3])<.97;if(dim(im)[3]==4)mask<-mask & im[,,4]>.05; rr<-which(rowSums(mask)>0);cc<-which(colSums(mask)>0);if(!length(rr))return(im);im[max(1,min(rr)-8):min(dim(im)[1],max(rr)+8),max(1,min(cc)-8):min(dim(im)[2],max(cc)+8),,drop=FALSE]}
img_cache<-new.env();pngim<-function(f){if(!exists(f,envir=img_cache,inherits=FALSE))assign(f,trim(readPNG(f)),envir=img_cache);get(f,envir=img_cache)}
cut<-function(f,x1,y1,x2,y2){im<-readPNG(f);trim(im[max(1,round(y1*dim(im)[1])):min(dim(im)[1],round(y2*dim(im)[1])),max(1,round(x1*dim(im)[2])):min(dim(im)[2],round(x2*dim(im)[2])),,drop=FALSE])}
place2<-function(obj,x,y,w,h,lab=NULL){
 pushViewport(viewport(x=x,y=1-y,width=w,height=h,just=c('left','top')))
 if(inherits(obj,'ggplot')) print(obj,vp=viewport(y=.47,height=.92,width=.97)) else {
  if(is.character(obj))obj<-pngim(obj)
  r<-dim(obj)[2]/dim(obj)[1];pr<-convertWidth(unit(1,'npc'),'mm',TRUE)/convertHeight(unit(1,'npc'),'mm',TRUE)
  grid.raster(obj,x=.5,y=.46,width=unit(min(.97,.89*r/pr),'npc'),height=unit(min(.89,.97*pr/r),'npc'))
 }
 if(!is.null(lab))grid.text(lab,x=.005,y=.997,just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=10,fontface='bold'));popViewport()
}
expo<-function(n,h,fun,w=183){stem<-file.path(out,paste0('Supplementary_Figure_S',n))
 svglite(paste0(stem,'.svg'),w/25.4,h/25.4);grid.newpage();fun();dev.off()
 if(Sys.info()[['sysname']]=='Darwin') quartz(type='pdf',file=paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='Arial') else grDevices::cairo_pdf(filename=paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='Arial');grid.newpage();fun();dev.off()
 agg_png(paste0(stem,'.png'),width=w,height=h,units='mm',res=300,background='white');grid.newpage();fun();dev.off()
 agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw',background='white');grid.newpage();fun();dev.off()
}
expo(1,220,function(){
 place2('04_Figures/R2_1_clustering/focused_2000_3000_5000/focused_PAC_by_k_preview.png',.01,.01,.48,.25,'a')
 place2('04_Figures/R2_1_clustering/focused_2000_3000_5000/focused_mean_silhouette_by_k_preview.png',.51,.01,.48,.25,'b')
 place2(margin_plot,.01,.28,.48,.39,'c');place2('04_Figures/R2_1_clustering/R2_1_k3_sensitivity_ARI_preview.png',.51,.28,.48,.39,'d')
 for(i in 1:3)place2(file.path(old,paste0('S1_consensus_',c(2000,3000,5000)[i],'.png')),.005+(i-1)/3,.69,.326,.30,letters[i+4])
})
# The updated S2 panel d matches the verified 0.193244335685485 external test.
updated<-file.path(old,'Supplementary_Figure_S2_updated.png')
em<-read_csv('09_Revision_Package/Supplementary_Data/R2_3_external_validation/external_mgs_figure_input.csv',show_col_types=FALSE)%>%mutate(Class=factor(Class,levels=c('Injury-like','nmSC Core')))
et<-read_csv('09_Revision_Package/Supplementary_Data/R2_3_external_validation/external_mgs_tests.csv',show_col_types=FALSE)
stopifnot(nrow(em)==88);for(ds in et$Dataset)stopifnot(abs(wilcox.test(MGS~Class,data=filter(em,Dataset==ds),exact=FALSE)$p.value-et$p[et$Dataset==ds])<1e-12)
et$label<-c('Wilcoxon P = 3.02 × 10⁻⁷','Wilcoxon P = 0.19')
emplot<-ggplot(em,aes(Class,MGS,fill=Class))+geom_boxplot(width=.58,outlier.shape=NA,linewidth=.3)+geom_point(position=position_jitter(width=.09,height=0,seed=123),size=.65,alpha=.82)+facet_wrap(~Dataset,nrow=1)+geom_text(data=et,aes(x=1.5,y=1.07,label=label),inherit.aes=FALSE,size=2.1)+scale_fill_manual(values=c('Injury-like'='#CD5C5C','nmSC Core'='#4682B4'),guide='none')+scale_y_continuous(limits=c(-.03,1.14),breaks=c(0,.25,.5,.75,1))+labs(x=NULL,y='Independently fitted MGS (0–1)',title='External MGS sensitivity')+theme(axis.text.x=element_text(size=6),plot.title=element_text(size=8),strip.text=element_text(size=7))
expo(2,225,function(){
 place2(panel('Fig2_composition_residualized_MGS'),.01,.01,.49,.26,'a');place2(label_plot,.51,.01,.48,.26,'b')
 place2(count_plot,.01,.29,.49,.26,'c');place2(emplot,.51,.29,.48,.26,'d')
 place2(file.path(old,'Supplementary_Figure_S8_external_functional_concordance.png'),.01,.57,.98,.42,'e')
})
expo(3,168,function(){
 place2(panel('Fig2_CIBERSORT_quality_pass_rate'),.01,.01,.38,.48,'a');place2(cp,.41,.01,.58,.48,'b')
 place2(panel('Fig3g_audited_GO_programs'),.01,.52,.58,.47,'c');place2('04_Figures/R2_4_immune_interpretation/R2_4_paired_CAMERA_summary_preview.png',.61,.52,.38,.47,'d')
})
# S4 is copied byte-for-byte unchanged.
file.copy(file.path(old,'Supplementary_Figure_S4.pdf'),file.path(out,'Supplementary_Figure_S4.pdf'),overwrite=TRUE)
file.copy(file.path(old,'Supplementary_Figure_S4_preview.png'),file.path(out,'Supplementary_Figure_S4.png'),overwrite=TRUE)
kw<-read_csv('03_Results/R2_6_variable_selection/R2_6_complete_subtype_screen_83_features.csv',show_col_types=FALSE)
co<-read_csv('03_Results/R2_6_variable_selection/R2_6_complete_clinical_correlation_screen_168_tests.csv',show_col_types=FALSE)
stopifnot(nrow(kw)==83,nrow(co)==168,length(unique(co$Feature))==84)
# Wrap complete screens into three adjacent blocks instead of narrow 83/84-row strips.
features<-sort(unique(kw$Feature)); cf<-sort(unique(co$Feature))
kwlim<-range(-log10(pmax(kw$FDR_BH_83,1e-300)),na.rm=TRUE);qlim<-range(-log10(pmax(co$Q_global_168_BH,1e-300)),na.rm=TRUE)
kp<-lapply(1:3,function(i){keep<-features[seq((i-1)*28+1,min(i*28,83))];dd<-kw%>%filter(Feature%in%keep)%>%mutate(Feature=factor(Feature,levels=rev(keep)))
 ggplot(dd,aes('Group test',Feature,fill=-log10(pmax(FDR_BH_83,1e-300))))+geom_tile(colour='white',linewidth=.15)+geom_text(aes(label=sprintf('q=%.3f',FDR_BH_83),colour=FDR_BH_83<.1),size=1.9)+scale_colour_manual(values=c('FALSE'='white','TRUE'='black'),guide='none')+scale_fill_viridis_c(option='C',limits=kwlim,name='−log10 BH q')+labs(x=NULL,y=NULL)+theme_classic(base_family='Arial',base_size=7)+theme(axis.text.y=element_text(size=6),axis.ticks=element_blank(),legend.position=if(i==3)'bottom' else 'none',legend.key.width=unit(4,'mm'),legend.key.height=unit(2,'mm'),plot.margin=margin(1,1,1,1))
})
corps<-lapply(1:3,function(i){keep<-cf[((i-1)*28+1):(i*28)];dd<-co%>%filter(Feature%in%keep)%>%mutate(Feature=factor(Feature,levels=rev(keep)),Endpoint=factor(Endpoint,levels=c('Tumor_texture_ordinal','Brainstem_adhesion'),labels=c('Texture','Adhesion')))
 ggplot(dd,aes(Endpoint,Feature,fill=Spearman_rho))+geom_tile(colour='white',linewidth=.15)+geom_point(aes(size=-log10(pmax(Q_global_168_BH,1e-300))),alpha=.65)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),name='Spearman ρ')+scale_size_continuous(name='−log10 global q',range=c(.2,1.4),limits=qlim)+labs(x=NULL,y=NULL)+theme_classic(base_family='Arial',base_size=7)+theme(axis.text.y=element_text(size=6),axis.ticks=element_blank(),legend.position=if(i==3)'bottom' else 'none',legend.key.width=unit(3,'mm'),legend.key.height=unit(2,'mm'),legend.text=element_text(size=5),legend.title=element_text(size=5),plot.margin=margin(1,1,1,1))+guides(fill=guide_colourbar(order=1),size=guide_legend(order=2,nrow=1))
})
expo(5,270,function(){
 place2(panel('Fig5g_adjusted_texture_adhesion_forest'),.01,.01,.57,.29,'a');place2(fp,.61,.01,.38,.29,'b')
 for(i in 1:3){place2(kp[[i]],.005+(i-1)/3,.32,.326,.31,if(i==1)'c' else NULL);place2(corps[[i]],.005+(i-1)/3,.66,.326,.33,if(i==1)'d' else NULL)}
},210)
# Reposition complete gating blocks; no gates, counts, axes, or events are edited.
p1<-file.path(old,'Supplementary_Figure_S6_Flow_Cytometry_Gating_page1.png');p2<-file.path(old,'Supplementary_Figure_S6_Flow_Cytometry_Gating_page2.png')
# Crop only surrounding blank page margins, maintaining original aspect ratio.
gates<-list(); im1<-readPNG(p1);im2<-readPNG(p2)
xb1<-c(32,326,655,972,1283);yb1<-c(55,407,756,1116,1472)
for(r in 1:4)for(c in 1:4){if(r==4&&c==1)next;gates[[length(gates)+1]]<-trim(im1[yb1[r]:(yb1[r+1]-1),xb1[c]:(xb1[c+1]-1),,drop=FALSE])}
xb2<-list(c(1,318,650,982,1290),c(15,334,666,1000),c(35,368,690));yb2<-c(85,490,939,1350)
for(r in 1:3)for(c in 1:c(4,3,2)[r])gates[[length(gates)+1]]<-trim(im2[yb2[r]:(yb2[r+1]-1),xb2[[r]][c]:(xb2[[r]][c+1]-1),,drop=FALSE])
stopifnot(length(gates)==24)
expo(6,275,function(){for(i in 1:24)place2(gates[[i]],.005+((i-1)%%4)*.25,.003+floor((i-1)/4)/6,.245,.164,letters[i])},210)
writeLines(c('S1: old S1 a-d + old S1 continued e-g','S2: old S2 a-d + old S8 as e','S3: old S3 a-b + old S7 a-b as c-d','S4: unchanged','S5: old a,d top; full b and c each wrapped into three blocks (83 /168 tests retained)','S6: both original gating sheets combined, isotropic scaling only'),file.path(out,'panel_mapping.txt'))
