# Purpose: Four score–purity correlations and partial-correlation sensitivity
# Manuscript/panels: Supplementary Table 6
# Reviewer comments: R1 Comment 5.6; R2 Comment 3
# Execution requires the inputs and working-directory layout described in README.md; "02_Reviewer_Analyses/10_Score_purity_sensitivity.R"
args <- commandArgs(trailingOnly=TRUE)
d <- read.csv('03_Results/R1_reviewer1_audit/R1_immune_program_scores_all38.csv')
stopifnot(nrow(d)==38,!anyDuplicated(d$SampleID), !anyNA(d))
out <- list()
for(scope in c('All38','C3')) {
 x <- if(scope=='C3') d[d$Subtype=='C3',] else d
 for(v in c('Immune_Recruitment_Score','Checkpoint_Score')) {
  z <- cor.test(x[[v]],x$Tumour_Purity,method='spearman',exact=FALSE)
  out[[length(out)+1]] <- data.frame(scope=scope,n=nrow(x),score=v,rho=unname(z$estimate),p=z$p.value)
 }
}
tab <- do.call(rbind,out);tab$q_BH_four_tests <- p.adjust(tab$p,'BH');write.csv(tab,'21_Score_Purity_Audit_20260925/score_purity_correlations.csv',row.names=FALSE)
x <- d[d$Subtype=='C3',]; a <- x$Immune_Recruitment_Score;b <- x$Checkpoint_Score;z <- x$Tumour_Purity
legacy <- cor.test(resid(lm(a~z)),resid(lm(b~z)),method='spearman',exact=FALSE)
ra <- rank(a);rb<-rank(b);rz<-rank(z);r <- cor(resid(lm(ra~rz)),resid(lm(rb~rz)));df<-nrow(x)-3;p<-2*pt(-abs(r*sqrt(df/(1-r*r))),df)
s <- data.frame(method=c('Historical_Spearman_of_raw_linear_residuals','Rank_partial_correlation_t_approximation'),rho=c(unname(legacy$estimate),r),p=c(legacy$p.value,p),n=nrow(x))
write.csv(s,'21_Score_Purity_Audit_20260925/purity_adjustment_sensitivity.csv',row.names=FALSE)
print(tab);print(s);cat('Suppression/checkpoint Spearman:',cor(d$Immune_Suppression_Score,d$Checkpoint_Score,method='spearman'),'\n');capture.output(sessionInfo(),file='21_Score_Purity_Audit_20260925/sessionInfo.txt')
