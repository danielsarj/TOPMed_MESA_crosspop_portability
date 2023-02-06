library(tidyverse)
library(data.table)
library(viridis)
"%&%" <- function(a,b) paste(a,b, sep='')
setwd('/home/daniel/MESA_heritability')

for (t in c('PBMC','Mono','Tcell')){
  for (p in c('AFA','EUR','HIS','CHN')){
    if (p == 'CHN' & t != 'PBMC'){
      next
    } else {
      if (file.exists(t %&% '_' %&% p %&% '_genes_that_went_wrong_constrained.txt.gz')==T){
        error_genes <- fread(t %&% '_' %&% p %&% '_genes_that_went_wrong_constrained.txt.gz') %>%
          mutate(tissue=t, pop=p)
        if (exists('error_genes_df')){
          error_genes_df <- rbind(error_genes_df, error_genes)
        } else {error_genes_df <- error_genes}
      }
      
      correct_genes <- fread(t %&% '_' %&% p %&% '_hsqestimates_constrained.txt.gz') %>%
        mutate(tissue=t, pop=p)
      if (exists('correct_genes_df')){
        correct_genes_df <- rbind(correct_genes_df, correct_genes)
      } else {correct_genes_df <- correct_genes}
    }
  }
}

# error vs correct
correct_n_nonc <- correct_genes_df_nonc %>% group_by(tissue, pop) %>% summarise(n = n()) %>% mutate(status = 'estimated')
error_n_nonc <- error_genes_df_nonc %>% group_by(tissue, pop) %>% summarise(n = n()) %>% mutate(status = 'error')
final_n_nonc <- rbind(correct_n_nonc, error_n_nonc)
correct_n <- correct_genes_df %>% group_by(tissue, pop) %>% summarise(n = n()) %>% mutate(status = 'estimated')
error_n <- error_genes_df %>% group_by(tissue, pop) %>% summarise(n = n()) %>% mutate(status = 'error')
final_n <- rbind(correct_n, error_n)
ggplot(final_n_nonc, aes(x=pop, y=n, fill=status)) + geom_col() + facet_wrap(~tissue) + scale_fill_viridis_d()
ggplot(final_n, aes(x=pop, y=n, fill=status)) + geom_col() + facet_wrap(~tissue) + scale_fill_viridis_d()

# all genes
ggplot(correct_genes_df_nonc, aes(x=h2, fill=pop)) + geom_density() + facet_wrap(~tissue) + scale_fill_viridis_d()
ggplot(correct_genes_df, aes(x=h2, fill=pop)) + geom_density() + facet_wrap(~tissue) + scale_fill_viridis_d()

ggplot(correct_genes_df_nonc, aes(x=pop, y=h2)) + geom_boxplot() + facet_wrap(~tissue) + scale_fill_viridis_d()

# only significant genes
correct_genes_sig_nonc <- correct_genes_df_nonc %>% filter(pval <= 0.05)
correct_genes_sig <- correct_genes_df %>% filter(pval <= 0.05)

correct_genes_sig_nonc %>% filter(tissue=='PBMC') %>% select(gene) %>% unique() %>% nrow()
correct_genes_sig %>% filter(tissue=='PBMC', pop=='AFA') %>% select(gene) %>% unique() %>% nrow()

ggplot(correct_genes_sig_nonc, aes(x=h2, fill=pop)) + geom_density(alpha=.5) + facet_wrap(~tissue) + scale_fill_viridis_d()
ggplot(correct_genes_sig, aes(x=h2, fill=pop)) + geom_density(alpha=.5) + facet_wrap(~tissue) + scale_fill_viridis_d()


ggplot(correct_genes_sig_nonc, aes(x=pop, y=h2)) + geom_violin() + geom_boxplot(width=0.1) + facet_wrap(~tissue) + scale_fill_viridis_d()


# comparing to GEUVADIS
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

# filter genes based on h2 and pval
correct_genes_df_pos_sig <- correct_genes_df %>% filter(h2 > 0 & pval <= 0.05)
correct_genes_df_sig <- correct_genes_df %>% filter(pval <= 0.05)

# get max h2 per tissue across populations
for (g in unique(correct_genes_df_pos_sig$gene)){
  tmp <- correct_genes_df_pos_sig %>% filter(gene==g) %>% group_by(tissue) %>% slice(which.max(h2))
  if (exists('genes_df_pos_sig_maxh2')){
    genes_df_pos_sig_maxh2 <- rbind(genes_df_pos_sig_maxh2, tmp)
  } else {genes_df_pos_sig_maxh2 <- tmp}
}
 
for (g in unique(correct_genes_df_sig$gene)){
  tmp <- correct_genes_df_sig %>% filter(gene==g) %>% group_by(tissue) %>% slice(which.max(h2))
  if (exists('genes_df_sig_maxh2')){
    genes_df_sig_maxh2 <- rbind(genes_df_sig_maxh2, tmp)
  } else {genes_df_sig_maxh2 <- tmp}
}

genes_df_pos_sig_maxh2_pbmc <- genes_df_pos_sig_maxh2 %>% filter(tissue=='PBMC')
mas_h2_pbmc <- inner_join(mas, genes_df_pos_sig_maxh2_pbmc, by=c('gene_id'='gene'))
mat_h2_pbmc <- inner_join(mat, genes_df_pos_sig_maxh2_pbmc, by=c('gene_id'='gene'))
en_h2_pbmc <- inner_join(en, genes_df_pos_sig_maxh2_pbmc, by=c('gene_id'='gene'))
ggplot(mat_h2_pbmc, aes(x=h2, y=matrixeqtl_estimate)) + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6) +
  geom_smooth(method='lm', color='red')
ggplot(mas_h2_pbmc, aes(x=h2, y=mashr_estimate)) + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6) +
  geom_smooth(method='lm', color='red')
ggplot(en_h2_pbmc, aes(x=h2, y=EN_estimate)) + geom_point() + facet_wrap(~mesa_pop+geu_pop, ncol=6) +
  geom_smooth(method='lm', color='red')

spearmans.scatter.df <- full_join(mas, mat, by = c("gene_id", "mesa_pop", "geu_pop")) %>% full_join(en, by = c("gene_id", "mesa_pop", "geu_pop"))

spearmans.scatter.df_pos_sig_maxh2 <- spearmans.scatter.df %>% filter(gene_id %in% unique(genes_df_pos_sig_maxh2_pbmc$gene))
ggplot(spearmans.scatter.df_pos_sig_maxh2, aes(x=mashr_estimate, y=EN_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")
ggplot(spearmans.scatter.df_pos_sig_maxh2, aes(x=mashr_estimate, y=matrixeqtl_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")

spearmans.scatter.df_pos_sig_maxh2_mashr_OR_EN <- spearmans.scatter.df_pos_sig_maxh2 %>% filter(mashr_estimate > 0 | EN_estimate > 0)
ggplot(spearmans.scatter.df_pos_sig_maxh2_mashr_OR_EN, aes(x=mashr_estimate, y=EN_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")

spearmans.scatter.df_pos_sig_maxh2_mashr_AND_EN <- spearmans.scatter.df_pos_sig_maxh2 %>% filter(mashr_estimate > 0 & EN_estimate > 0)
ggplot(spearmans.scatter.df_pos_sig_maxh2_mashr_AND_EN, aes(x=mashr_estimate, y=EN_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")

### 
for (i in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
genes_df <- correct_genes_df %>% filter(h2 > i & pval <= 0.05) %>% filter(tissue=='PBMC')
spearmans.scatter.df_mod <- spearmans.scatter.df %>% filter(gene_id %in% unique(genes_df$gene))

ggplot(spearmans.scatter.df_mod, aes(x=mashr_estimate, y=EN_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")
#ggsave('/home/daniel/MESA_heritability/plots/Geuvadis_scatterplot_MASHRvsEN_allSpearman_sigH2_H2' %&% i %&%'.png', height = 8, width = 10)
#ggsave('/home/daniel/MESA_heritability/plots/Geuvadis_scatterplot_MASHRvsEN_allSpearman_sigH2_H2' %&% i %&%'.pdf', height = 8, width = 10)
#ggplot(spearmans.scatter.df_mod, aes(x=mashr_estimate, y=matrixeqtl_estimate)) + geom_point() + 
#  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")

spearmans.scatter.df_mod_mashr_OR_EN <- spearmans.scatter.df_mod  %>% filter(mashr_estimate > 0 | EN_estimate > 0)
ggplot(spearmans.scatter.df_mod_mashr_OR_EN, aes(x=mashr_estimate, y=EN_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")
#ggsave('/home/daniel/MESA_heritability/plots/Geuvadis_scatterplot_MASHRvsEN_ORSpearman_sigH2_H2' %&% i %&%'.png', height = 8, width = 10)
#ggsave('/home/daniel/MESA_heritability/plots/Geuvadis_scatterplot_MASHRvsEN_ORSpearman_sigH2_H2' %&% i %&%'.pdf', height = 8, width = 10)

spearmans.scatter.df_mod_mashr_AND_EN <- spearmans.scatter.df_mod %>% filter(mashr_estimate > 0 & EN_estimate > 0)
ggplot(spearmans.scatter.df_mod_mashr_AND_EN, aes(x=mashr_estimate, y=EN_estimate)) + geom_point() + 
  facet_wrap(~mesa_pop+geu_pop, ncol=6) + geom_smooth(method='lm', color='red') + geom_abline(color="blue")
#ggsave('/home/daniel/MESA_heritability/plots/Geuvadis_scatterplot_MASHRvsEN_ANDSpearman_sigH2_H2' %&% i %&%'.png', height = 8, width = 10)
#ggsave('/home/daniel/MESA_heritability/plots/Geuvadis_scatterplot_MASHRvsEN_ANDSpearman_sigH2_H2' %&% i %&%'.pdf', height = 8, width = 10)
}
