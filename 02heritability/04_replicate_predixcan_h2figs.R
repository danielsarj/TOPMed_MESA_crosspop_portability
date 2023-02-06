library(tidyverse)
library(ggrepel)
library(data.table)
library(viridis)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep='')
setwd('/home/daniel/MESA_heritability')

parser <- ArgumentParser()
parser$add_argument("--r2", help="LD r2 to use")
parser$add_argument("--mode", help="GCTA REML mode used")
args <- parser$parse_args()

for (t in c('PBMC','Mono','Tcell')){
  for (p in c('AFA','EUR','HIS','CHN')){
    if (p == 'CHN' & t != 'PBMC'){
      next
    } else {
      correct_genes <- fread(t %&% '_' %&% p %&% '_hsqestimates_r' %&% as.character(r2) %&% '_' %&% mode %&% '.txt.gz') %>%
        mutate(tissue=t, pop=p)
      if (exists('correct_genes_df')){
        correct_genes_df <- rbind(correct_genes_df, correct_genes)
      } else {correct_genes_df <- correct_genes}
    }
  }
}

correct_genes_df_sig <- correct_genes_df %>% filter(pval <= 0.05)
fwrite(correct_genes_df_sig, 'plots/significant_h2estimates_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'.txt', col.names=T, sep='\t', quote=F)

for (m_pop in c('AFA','EUR','HIS','CHN')){
  for (method in c('ENunf', 'mashr', 'matrixeqtl')){
    for (g_pop in c('GBR', 'YRI', 'ALL', 'CEU', 'FIN', 'TSI')){
      
      spearman_f <- fread('/home/daniel/Geuvadis/WGS_predixcan/Geuvadis.'%&% g_pop %&%'_PBMC.'%&% m_pop %&%'_'%&% method %&%'_expression_spearman_correlation.txt') %>% 
        filter(!is.na(estimate)) %>% select(gene_id, estimate) %>% mutate(model=method, mesa_pop=m_pop, geu_pop=g_pop)
      
      if (exists('spearmans.df')){
        spearmans.df <- rbind(spearmans.df, spearman_f)
      } else {spearmans.df <- spearman_f}
    }
  }
}

mas <- spearmans.df %>% filter(model=='mashr') %>% rename(mashr_estimate=estimate)
mat <- spearmans.df %>% filter(model=='matrixeqtl') %>% rename(matrixeqtl_estimate=estimate)
en <- spearmans.df %>% filter(model=='ENunf') %>% rename(EN_estimate=estimate)

### making h2 plots per population

genes_df_sig_pbmc_perpop <- correct_genes_df_sig %>% filter(tissue=='PBMC') %>% group_by(pop) %>% arrange(h2) %>% 
  mutate(ymin = pmax(h2-2*se), ymax = pmin(h2+2*se), position=1:n(), lower_bound = NA)

for (i in 1:nrow(genes_df_sig_pbmc_perpop)){
  if (genes_df_sig_pbmc_perpop$ymin[i] <= 0.01){
    genes_df_sig_pbmc_perpop$lower_bound[i] <- '<= 0.01'
  } else {genes_df_sig_pbmc_perpop$lower_bound[i] <- '> 0.01'}
}

h2plot <- ggplot(genes_df_sig_pbmc_perpop, aes(x=position, y=h2, ymin=ymin, ymax=ymax, color=lower_bound)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~pop) + scale_color_viridis_d(direction=-1) +
  labs(x='Genes ordered by h²', y='Gene expression cis-heritability (h²)', col='h² lower bound')
h2plot
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% mode %&%'_r' %&% as.character(r2) %&%'.png', width=8, height=5)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% mode %&%'_r' %&% as.character(r2) %&%'.pdf', width=8, height=5)
h2_summ <- genes_df_sig_pbmc_perpop %>% group_by(pop, lower_bound) %>% summarise(n=n(), avg=mean(h2), median=median(h2)) 
fwrite(h2_summ, '/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_summary.txt', col.names=T, sep='\t', quote=F)

summary_df <- genes_df_sig_pbmc_perpop %>% filter(lower_bound=='> 0.01') %>% group_by(pop) %>% summarise(n_genes = n(), classification='sig and her')
summary_df$sample_size <- c(334, 104, 528, 312)
correlation <- cor.test(summary_df$n_genes, summary_df$sample_size, method = 'pearson')
ggplot(summary_df, aes(x=sample_size, y=n_genes, color=classification)) + geom_point(size=4) + 
  geom_label_repel(label=summary_df$pop) + scale_color_viridis_d() + 
  labs(x='Sample size', y='Number of genes with positive h²')
ggsave('/home/daniel/MESA_heritability/plots/PBMC_popsample_vs_ngenes_' %&% mode %&%'_r' %&% as.character(r2) %&%'.png', width=8, height=5)
ggsave('/home/daniel/MESA_heritability/plots/PBMC_popsample_vs_ngenes_' %&% mode %&%'_r' %&% as.character(r2) %&%'.pdf', width=8, height=5)

genes_df_sig_pbmc_perpop_matrix_wgeuvadis <- inner_join(mat, genes_df_sig_pbmc_perpop,  by=c('gene_id'='gene', 'mesa_pop'='pop'))
genes_df_sig_pbmc_perpop_mashr_wgeuvadis <- inner_join(mas, genes_df_sig_pbmc_perpop,  by=c('gene_id'='gene', 'mesa_pop'='pop'))
genes_df_sig_pbmc_perpop_en_wgeuvadis <- inner_join(en, genes_df_sig_pbmc_perpop,  by=c('gene_id'='gene', 'mesa_pop'='pop'))

h2plot_mashr <- ggplot(genes_df_sig_pbmc_perpop_mashr_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop+lower_bound, ncol=6)
h2plot_mashr + geom_point(data=genes_df_sig_pbmc_perpop_mashr_wgeuvadis, aes(x=position, y=mashr_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop+lower_bound, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_splitbylowerbound_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMASHRestimatesinGeuvadis.png', width=12, height=20)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_splitbylowerbound_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMASHRestimatesinGeuvadis.pdf', width=12, height=20)
h2plot_mashr <- ggplot(genes_df_sig_pbmc_perpop_mashr_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6)
h2plot_mashr + geom_point(data=genes_df_sig_pbmc_perpop_mashr_wgeuvadis, aes(x=position, y=mashr_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMASHRestimatesinGeuvadis.png', width=12, height=10)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMASHRestimatesinGeuvadis.pdf', width=12, height=10)

h2plot_matrix <- ggplot(genes_df_sig_pbmc_perpop_matrix_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop+lower_bound, ncol=6)
h2plot_matrix + geom_point(data=genes_df_sig_pbmc_perpop_matrix_wgeuvadis, aes(x=position, y=matrixeqtl_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop+lower_bound, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_splitbylowerbound_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMATRIXestimatesinGeuvadis.png', width=12, height=20)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_splitbylowerbound_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMATRIXestimatesinGeuvadis.pdf', width=12, height=20)
h2plot_matrix <- ggplot(genes_df_sig_pbmc_perpop_matrix_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6)
h2plot_matrix + geom_point(data=genes_df_sig_pbmc_perpop_matrix_wgeuvadis, aes(x=position, y=matrixeqtl_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMATRIXestimatesinGeuvadis.png', width=12, height=10)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMATRIXestimatesinGeuvadis.pdf', width=12, height=10)

h2plot_en <- ggplot(genes_df_sig_pbmc_perpop_en_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop+lower_bound, ncol=6)
h2plot_en + geom_point(data=genes_df_sig_pbmc_perpop_en_wgeuvadis, aes(x=position, y=EN_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop+lower_bound, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_splitbylowerbound_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wENestimatesinGeuvadis.png', width=12, height=20)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_splitbylowerbound_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wENestimatesinGeuvadis.pdf', width=12, height=20)
h2plot_en <- ggplot(genes_df_sig_pbmc_perpop_en_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6)
h2plot_en + geom_point(data=genes_df_sig_pbmc_perpop_en_wgeuvadis, aes(x=position, y=EN_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wENestimatesinGeuvadis.png', width=12, height=10)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_perpop_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wENestimatesinGeuvadis.pdf', width=12, height=10)

### making plots taking the highest h2 across populations

for (g in unique(correct_genes_df_sig$gene)){
  tmp <- correct_genes_df_sig %>% filter(gene==g) %>% group_by(tissue) %>% slice(which.max(h2))
  if (exists('genes_df_sig_maxh2')){
    genes_df_sig_maxh2 <- rbind(genes_df_sig_maxh2, tmp)
  } else {genes_df_sig_maxh2 <- tmp}
}
genes_df_sig_maxh2_pbmc <- genes_df_sig_maxh2 %>% filter(tissue=='PBMC') %>% arrange(h2) %>% 
  mutate(ymin = pmax(h2 - 2 * se), ymax = pmin(h2 + 2 * se))
genes_df_sig_maxh2_pbmc <- genes_df_sig_maxh2_pbmc %>% mutate(position=1:nrow(genes_df_sig_maxh2_pbmc))

h2plot <- ggplot(genes_df_sig_maxh2_pbmc, aes(x=position, y=h2, ymin=ymin, ymax=ymax) ) + 
  geom_pointrange(col='gray') + geom_point()
h2plot
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'.png', width=6, height=4)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'.pdf', width=6, height=4)

genes_df_sig_maxh2_pbmc_matrix_wgeuvadis <- inner_join(mat, genes_df_sig_maxh2_pbmc,  by=c('gene_id'='gene'))
genes_df_sig_maxh2_pbmc_mashr_wgeuvadis <- inner_join(mas, genes_df_sig_maxh2_pbmc,  by=c('gene_id'='gene'))
genes_df_sig_maxh2_pbmc_en_wgeuvadis <- inner_join(en, genes_df_sig_maxh2_pbmc,  by=c('gene_id'='gene'))

h2plot_mashr <- ggplot(genes_df_sig_maxh2_pbmc_mashr_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6)
h2plot_mashr + geom_point(data=genes_df_sig_maxh2_pbmc_mashr_wgeuvadis, aes(x=position, y=mashr_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMASHRestimatesinGeuvadis.png', width=12, height=10)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMASHRestimatesinGeuvadis.pdf', width=12, height=10)

h2plot_matrix <- ggplot(genes_df_sig_maxh2_pbmc_matrix_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6)
h2plot_matrix + geom_point(data=genes_df_sig_maxh2_pbmc_matrix_wgeuvadis, aes(x=position, y=matrixeqtl_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMATRIXestimatesinGeuvadis.png', width=12, height=10)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wMATRIXestimatesinGeuvadis.pdf', width=12, height=10)

h2plot_en <- ggplot(genes_df_sig_maxh2_pbmc_en_wgeuvadis, aes(x=position, y=h2, ymin=ymin, ymax=ymax)) + 
  geom_pointrange(col='gray') + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6)
h2plot_en + geom_point(data=genes_df_sig_maxh2_pbmc_en_wgeuvadis, aes(x=position, y=EN_estimate), color='red', size=0.4, alpha=0.1) + facet_wrap(~mesa_pop+geu_pop, ncol=6)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wENestimatesinGeuvadis.png', width=12, height=10)
ggsave('/home/daniel/MESA_heritability/plots/h2estimates_PBMC_' %&% args$mode %&%'_r' %&% as.character(args$r2) %&%'_wENestimatesinGeuvadis.pdf', width=12, height=10)
