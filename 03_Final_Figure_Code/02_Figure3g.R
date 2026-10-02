# Purpose: NES lollipop plot of all 12 terms passing the original cutoff
# Manuscript/panels: Figure 3g
# Reviewer comments: R2 Comment 5
# Execution requires the inputs and working-directory layout described in README.md; "03_Final_Figure_Code/02_Figure3g.R"
# Replot the reproduced C3-high vs C3-low GO-BP GSEA results in the historical NES lollipop style.
# Input: GSEA_original_code_pvalueCutoff_0.05.csv (12 terms from the reproduced run)
# Output: separate revision files; no original figure is overwritten.

suppressPackageStartupMessages({library(ggplot2); library(stringr)})
root <- file.path(Sys.getenv('VS_SERVER_ROOT', unset='inputs/legacy_project'), 'Reproduction_Fig3g_20261001')
res <- read.csv(file.path(root, 'GSEA_original_code_pvalueCutoff_0.05.csv'), check.names=FALSE)
stopifnot(nrow(res) == 12L, all(c('Description','NES','p.adjust') %in% names(res)))
stopifnot(all(is.finite(res$NES)), all(is.finite(res$p.adjust)), all(res$p.adjust > 0))

# Preserve all 12 terms returned by the original gseGO cutoff; sort by NES as in the prior panel.
plot_df <- res
plot_df$logFDR <- -log10(plot_df$p.adjust)
plot_df$Description <- stringr::str_wrap(plot_df$Description, width=37)
plot_df <- plot_df[order(plot_df$NES), , drop=FALSE]
plot_df$Description <- factor(plot_df$Description, levels=plot_df$Description)
max_nes <- max(abs(plot_df$NES), na.rm=TRUE) * 1.3

p <- ggplot(plot_df, aes(x=NES, y=Description)) +
  geom_vline(xintercept=0, linetype='dashed', color='grey30', linewidth=0.6) +
  geom_segment(aes(x=0, xend=NES, y=Description, yend=Description),
               color='grey70', linewidth=1.2) +
  geom_point(aes(size=logFDR, fill=NES), shape=21, color='black', stroke=0.6) +
  scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B', midpoint=0, name='NES') +
  scale_size_continuous(range=c(3,8), name=expression(-log[10](FDR))) +
  scale_x_continuous(limits=c(-max_nes,max_nes), expand=expansion(mult=0.02)) +
  theme_classic(base_size=12, base_family='Helvetica') +
  labs(x='Normalized Enrichment Score (NES)', y=NULL,
       title='Key GO Pathways in C3-high Myeloid') +
  theme(plot.title=element_text(face='bold', hjust=0.5, size=13),
        axis.text.y=element_text(color='black', face='bold', size=10, lineheight=0.8),
        axis.text.x=element_text(color='black'),
        axis.title.x=element_text(face='bold', margin=ggplot2::margin(t=10)),
        legend.position='right', legend.background=element_blank(),
        legend.box.margin=ggplot2::margin(0,0,0,-10),
        panel.grid.major.y=element_line(color='grey95', linetype='dotted'))

# Preserve a tidy plotting table as source data for the panel.
write.csv(plot_df, file.path(root, 'Fig3g_panel_source_data.csv'), row.names=FALSE)
grDevices::pdf(file.path(root, 'Fig3g_GSEA_rerun_panel.pdf'), width=10, height=6.5, family='Helvetica', useDingbats=FALSE)
print(p)
grDevices::dev.off()
ggsave(file.path(root, 'Fig3g_GSEA_rerun_panel.png'), p, width=10, height=6.5, dpi=600, bg='white')
print(p)
cat('Panel plotted from', nrow(plot_df), 'terms; NES-positive:', sum(plot_df$NES>0),
    'NES-negative:', sum(plot_df$NES<0), '\n')
