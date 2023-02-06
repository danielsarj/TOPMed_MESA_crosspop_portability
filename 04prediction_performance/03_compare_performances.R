library(data.table)
library(tidyverse)
library(ggpubr)
library(viridis)
'%&%' = function(a,b) paste(a,b,sep='')
setwd('/home/daniel/Geuvadis/WGS_predixcan')

# reading significant non-zero heritability estimates
h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01)

# reading geuvadis outputs
for (tis in c('PBMC','Mono','Tcell')){
  for (m_pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & m_pop=='CHN'){
      next
    } else {
      for (method in c('ENunf', 'mashr', 'matrixeqtl')){
        for (g_pop in c('GBR', 'YRI', 'ALL', 'CEU', 'FIN', 'TSI')){
          
          g_output <- fread('Geuvadis.'%&% g_pop %&%'_'%&% tis %&%'.'%&% m_pop %&%'_'%&% method %&%'_expression_spearman_correlation.txt') %>% 
            filter(!is.na(estimate)) %>% select(gene_id, estimate) %>% mutate(tissue=tis, model=method, mesa_pop=m_pop, geu_pop=g_pop)
          
          if (exists('geuvadis.spearmans.df')){
            geuvadis.spearmans.df <- rbind(geuvadis.spearmans.df, g_output)
          } else {geuvadis.spearmans.df <- g_output}
        }
      }
    }
  }
}

# filtering geuvadis outputs so they only contain genes in the h2 dataframe
for (tis in c('PBMC','Mono','Tcell')){
  genes.in.tissue <- h2estimates %>% filter(tissue==tis) %>% select(gene) %>% unique()
  geuvadis.filt <- geuvadis.spearmans.df %>% filter(tissue==tis, gene_id %in% genes.in.tissue$gene)
  
  if (exists('geuvadis.spearmans.df.filtered')){
    geuvadis.spearmans.df.filtered <- rbind(geuvadis.spearmans.df.filtered, geuvadis.filt)
  } else {geuvadis.spearmans.df.filtered <- geuvadis.filt}
}

# computing median spearman estimate by tissue/model/MESA populatio/Geuvadis population
geuvadis.median.spearmans.df <- geuvadis.spearmans.df.filtered %>% group_by(tissue, model, mesa_pop, geu_pop) %>%
  select(estimate) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm=T))
geuvadis.median.spearmans.df$model <- gsub('ENunf', 'EN', geuvadis.median.spearmans.df$model)
geuvadis.median.spearmans.df$model <- gsub('mashr', 'MASHR', geuvadis.median.spearmans.df$model)
geuvadis.median.spearmans.df$model <- gsub('matrixeqtl', 'MatrixeQTL', geuvadis.median.spearmans.df$model)
fwrite(geuvadis.median.spearmans.df, 'Geuvadis_median_spearman_heatmap_filteredbyh2.txt', col.names=T, sep='\t', quote=F)

# making heatmap
ggplot(geuvadis.median.spearmans.df, aes(x=model, y=geu_pop, fill=median_estimate)) + facet_wrap(~tissue+mesa_pop, ncol=4) +
  geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis populations') + xlab('Models') + ylab('Geuvadis populations') +
  theme(axis.text.x = element_text(angle = 40, hjust=1))
ggsave('Geuvadis_median_spearman_heatmap_filteredbyh2.pdf', height=8, width=10)
ggsave('Geuvadis_median_spearman_heatmap_filteredbyh2.tiff', height=8, width=10)  
ggsave('Geuvadis_median_spearman_heatmap_filteredbyh2.png', height=8, width=10)  

# reformating geuvadis dataframe
for (method in c('ENunf', 'mashr', 'matrixeqtl')){
  tmp <- geuvadis.spearmans.df.filtered %>% filter(model==method)
  colnames(tmp)[2] <- c(method%&%'_estimate')
  
  if (exists('geuvadis.spearmans.reformatted.df')){
    geuvadis.spearmans.reformatted.df <- full_join(geuvadis.spearmans.reformatted.df, tmp, by = c('gene_id', 'tissue', 'mesa_pop', 'geu_pop'))
  } else {geuvadis.spearmans.reformatted.df<-tmp}
}
geuvadis.spearmans.reformatted.df <- geuvadis.spearmans.reformatted.df %>% select(gene_id, tissue, mesa_pop, geu_pop, ENunf_estimate, mashr_estimate, matrixeqtl_estimate)

# wilcox test to see if models are significantly different
for (tis in c('PBMC','Mono','Tcell')){
  for (m_pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & m_pop=='CHN'){
      next
    } else {
      for (g_pop in c('GBR', 'YRI', 'ALL', 'CEU', 'FIN', 'TSI')){
        
        mashr_matrix_all <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, matrixeqtl_estimate) %>% melt() %>% drop_na()
        test <- wilcox.test(value~variable, data=mashr_matrix_all)
        mashr_matrix_all <- test$p.value
        
        mashr_matrix_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, matrixeqtl_estimate) %>% drop_na() %>% melt() 
        test <- wilcox.test(value~variable, data=mashr_matrix_intersection)
        mashr_matrix_intersection <- test$p.value
        
        mashr_EN_all <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, ENunf_estimate) %>% melt() %>% drop_na()
        test <- wilcox.test(value~variable, data=mashr_EN_all)
        mashr_EN_all <- test$p.value
        
        mashr_EN_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, ENunf_estimate) %>% drop_na() %>% melt() 
        test <- wilcox.test(value~variable, data=mashr_EN_intersection)
        mashr_EN_intersection <- test$p.value
        
        matrix_EN_all <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(matrixeqtl_estimate, ENunf_estimate) %>% melt() %>% drop_na()
        test <- wilcox.test(value~variable, data=matrix_EN_all)
        matrix_EN_all <- test$p.value
        
        matrix_EN_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(matrixeqtl_estimate, ENunf_estimate) %>% drop_na() %>% melt() 
        test <- wilcox.test(value~variable, data=matrix_EN_intersection)
        matrix_EN_intersection <- test$p.value
        
        new_line <- c(tis, m_pop, g_pop, mashr_matrix_all, mashr_EN_all, matrix_EN_all, mashr_matrix_intersection, mashr_EN_intersection, matrix_EN_intersection)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
  }
}
colnames(wilcox.tests.pvalues) <- c('tissue','mesa_pop','geu_pop', 'mashr_matrix_all', 'mashr_EN_all', 'matrix_EN_all',  'mashr_matrix_intersection', 'mashr_EN_intersection', 'matrix_EN_intersection')
wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
wilcox.tests.pvalues$mashr_matrix_all <- as.numeric(wilcox.tests.pvalues$mashr_matrix_all)
wilcox.tests.pvalues$mashr_EN_all <- as.numeric(wilcox.tests.pvalues$mashr_EN_all)
wilcox.tests.pvalues$matrix_EN_all <- as.numeric(wilcox.tests.pvalues$matrix_EN_all)
wilcox.tests.pvalues$mashr_matrix_intersection <- as.numeric(wilcox.tests.pvalues$mashr_matrix_intersection)
wilcox.tests.pvalues$mashr_EN_intersection <- as.numeric(wilcox.tests.pvalues$mashr_EN_intersection)
wilcox.tests.pvalues$matrix_EN_intersection <- as.numeric(wilcox.tests.pvalues$matrix_EN_intersection)
fwrite(wilcox.tests.pvalues, 'Geuvadis_wilcox_spearman_tests_pvalues_filteredbyh2.txt', col.names=T, sep='\t', quote=F)

### taking a look at only GBR and YRI

# filtering by geuvadis population
gbr.yri.spermans.df <- geuvadis.spearmans.df.filtered %>% filter(geu_pop=='GBR' | geu_pop=='YRI')
gbr.yri.spermans.df$model <- gsub('ENunf', 'EN', gbr.yri.spermans.df$model)
gbr.yri.spermans.df$model <- gsub('mashr', 'MASHR', gbr.yri.spermans.df$model)
gbr.yri.spermans.df$model <- gsub('matrixeqtl', 'MatrixeQTL', gbr.yri.spermans.df$model)

# violin plots per tissue
for (tis in unique(gbr.yri.spermans.df$tissue)){
  # all genes and models
  gbr.yri.spermans.df %>% filter(tissue==tis) %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') 
  
  # all genes and EN and MASHR
  gbr.yri.spermans.df %>% filter(tissue==tis, model=='EN' | model =='MASHR') %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') + stat_compare_means(method='wilcox.test', label='p.signif')
  if (tis!='PBMC'){
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMASHR_violinplot_allgenes_filteredbyh2.pdf', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMASHR_violinplot_allgenes_filteredbyh2.tiff', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMASHR_violinplot_allgenes_filteredbyh2.png', height=4, width=8)  
  } else {
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMASHR_violinplot_allgenes_filteredbyh2.pdf', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMASHR_violinplot_allgenes_filteredbyh2.tiff', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMASHR_violinplot_allgenes_filteredbyh2.png', height=6, width=8)  
  }
  
  # all genes and EN and MatrixeQTL
  gbr.yri.spermans.df %>% filter(tissue==tis, model=='EN' | model =='MatrixeQTL') %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') + stat_compare_means(method='wilcox.test', label='p.signif')
  if (tis!='PBMC'){
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMatrixeQTL_violinplot_allgenes_filteredbyh2.pdf', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMatrixeQTL_violinplot_allgenes_filteredbyh2.tiff', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMatrixeQTL_violinplot_allgenes_filteredbyh2.png', height=4, width=8)  
  } else {
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMatrixeQTL_violinplot_allgenes_filteredbyh2.pdf', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMatrixeQTL_violinplot_allgenes_filteredbyh2.tiff', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMatrixeQTL_violinplot_allgenes_filteredbyh2.png', height=6, width=8)  
  }
  
  # all genes and MASHR and MatrixeQTL
  gbr.yri.spermans.df %>% filter(tissue==tis, model=='MASHR' | model =='MatrixeQTL') %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') + stat_compare_means(method='wilcox.test', label='p.signif')
  if (tis!='PBMC'){
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_MASHRandMatrixeQTL_violinplot_allgenes_filteredbyh2.pdf', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_MASHRandMatrixeQTL_violinplot_allgenes_filteredbyh2.tiff', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_MASHRandMatrixeQTL_violinplot_allgenes_filteredbyh2.png', height=4, width=8)  
  } else {
    ggsave('Geuvadis_PBMC_GBR_YRI_MASHRandMatrixeQTL_violinplot_allgenes_filteredbyh2.pdf', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_MASHRandMatrixeQTL_violinplot_allgenes_filteredbyh2.tiff', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_MASHRandMatrixeQTL_violinplot_allgenes_filteredbyh2.png', height=6, width=8)  
  }
  
  # assessing genes that are shared between models and populations
  for (p in c('GBR', 'YRI')){
    for (m in c('AFA','EUR','HIS','CHN')){
      if (tis!='PBMC' & m=='CHN'){
        next
      } else {
        en.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='EN', geu_pop==p, mesa_pop==m)
        mashr.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='MASHR', geu_pop==p, mesa_pop==m)
        matrixeqtl.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='MatrixeQTL', geu_pop==p, mesa_pop==m)
        
        en.mashr.shared.genes <- intersect(en.genes$gene_id, mashr.genes$gene_id) 
        en.mashr.joint.df <- rbind(en.genes, mashr.genes) %>% filter(gene_id %in% en.mashr.shared.genes)
        
        en.matrixeqtl.shared.genes <- intersect(en.genes$gene_id, matrixeqtl.genes$gene_id) 
        en.matrixeqtl.joint.df <- rbind(en.genes, matrixeqtl.genes) %>% filter(gene_id %in% en.matrixeqtl.shared.genes)
        
        mashr.matrixeqtl.shared.genes <- intersect(mashr.genes$gene_id, matrixeqtl.genes$gene_id) 
        mashr.matrixeqtl.joint.df <- rbind(mashr.genes, matrixeqtl.genes) %>% filter(gene_id %in% mashr.matrixeqtl.shared.genes)
        
        if (exists('shared.genes.en.mashr')){
          shared.genes.en.mashr <- rbind(shared.genes.en.mashr, en.mashr.joint.df)
        } else {shared.genes.en.mashr <- en.mashr.joint.df}
        
        if (exists('shared.genes.en.matrixeqtl')){
          shared.genes.en.matrixeqtl <- rbind(shared.genes.en.matrixeqtl, en.matrixeqtl.joint.df)
        } else {shared.genes.en.matrixeqtl <- en.matrixeqtl.joint.df}
        
        if (exists('shared.genes.mashr.matrixeqtl')){
          shared.genes.mashr.matrixeqtl <- rbind(shared.genes.mashr.matrixeqtl, mashr.matrixeqtl.joint.df)
        } else {shared.genes.mashr.matrixeqtl <- mashr.matrixeqtl.joint.df}
      }
    }
  }
  
  # shared genes and EN and MASHR
  shared.genes.en.mashr %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') + stat_compare_means(method='wilcox.test', label='p.signif')
  if (tis!='PBMC'){
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMASHR_violinplot_sharedgenes_filteredbyh2.pdf', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMASHR_violinplot_sharedgenes_filteredbyh2.tiff', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMASHR_violinplot_sharedgenes_filteredbyh2.png', height=4, width=8)  
  } else {
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMASHR_violinplot_sharedgenes_filteredbyh2.pdf', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMASHR_violinplot_sharedgenes_filteredbyh2.tiff', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMASHR_violinplot_sharedgenes_filteredbyh2.png', height=6, width=8)  
  }
  
  # shared genes and EN and MatrixeQTL
  shared.genes.en.matrixeqtl %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') + stat_compare_means(method='wilcox.test', label='p.signif')
  if (tis!='PBMC'){
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.pdf', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.tiff', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_ENandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.png', height=4, width=8)  
  } else {
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.pdf', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.tiff', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_ENandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.png', height=6, width=8)  
  }
  
  # all genes and MASHR and MatrixeQTL
  shared.genes.mashr.matrixeqtl %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
    geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
    labs(title='Spearman correlation in Geuvadis ('%&% tis %&%')', fill='Model') + ylab('Spearman correlation') + xlab('Geuvadis populations') + stat_compare_means(method='wilcox.test', label='p.signif')
  if (tis!='PBMC'){
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_MASHRandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.pdf', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_MASHRandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.tiff', height=4, width=8)  
    ggsave('Geuvadis_'%&% tis %&%'_GBR_YRI_MASHRandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.png', height=4, width=8)  
  } else {
    ggsave('Geuvadis_PBMC_GBR_YRI_MASHRandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.pdf', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_MASHRandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.tiff', height=6, width=8)  
    ggsave('Geuvadis_PBMC_GBR_YRI_MASHRandMatrixeQTL_violinplot_sharedgenes_filteredbyh2.png', height=6, width=8)  
  }
  
  rm(shared.genes.en.mashr)
  rm(shared.genes.en.matrixeqtl)
  rm(shared.genes.mashr.matrixeqtl)
}
