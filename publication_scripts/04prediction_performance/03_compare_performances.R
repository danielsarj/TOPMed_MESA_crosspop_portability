library(data.table)
library(tidyverse)
library(ggpubr)
library(viridis)
'%&%' = function(a,b) paste(a,b,sep='')
setwd('/home/daniel/Geuvadis/WGS_predixcan')

# reading significant non-zero heritability estimates
h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01)

# reading geuvadis outputs
for (tis in c('PBMC', 'Mono', 'Tcell')){
  for (m_pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & m_pop=='CHN'){
      next
    } else {
      for (method in c('ENunf', 'mashr', 'matrixeqtl', 'TIGAR_1e-04', 'JTI')){
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
geuvadis.median.spearmans.df$model <- gsub('TIGAR_1e-04', 'TIGAR', geuvadis.median.spearmans.df$model)
#fwrite(geuvadis.median.spearmans.df, 'Geuvadis_median_spearman_heatmap_filteredbyh2.txt', col.names=T, sep='\t', quote=F)

# making heatmap
ggplot(geuvadis.median.spearmans.df, aes(x=model, y=geu_pop, fill=median_estimate)) + facet_wrap(~tissue+mesa_pop, ncol=4) +
  geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis populations', fill='Value') + xlab('Models') + ylab('Geuvadis populations') +
  theme(axis.text.x = element_text(angle = 40, hjust=1))
#ggsave('Geuvadis_median_spearman_heatmap_filteredbyh2.pdf', height=8, width=10)
#ggsave('Geuvadis_median_spearman_heatmap_filteredbyh2.tiff', height=8, width=10)  
#ggsave('Geuvadis_median_spearman_heatmap_filteredbyh2.png', height=8, width=10)  

# reformating geuvadis dataframe
for (method in c('ENunf', 'mashr', 'matrixeqtl', 'TIGAR_1e-04', 'JTI')){
  tmp <- geuvadis.spearmans.df.filtered %>% filter(model==method)
  colnames(tmp)[2] <- c(method%&%'_estimate')
  
  if (exists('geuvadis.spearmans.reformatted.df')){
    geuvadis.spearmans.reformatted.df <- full_join(geuvadis.spearmans.reformatted.df, tmp, by = c('gene_id', 'tissue', 'mesa_pop', 'geu_pop'))
  } else {geuvadis.spearmans.reformatted.df<-tmp}
}
geuvadis.spearmans.reformatted.df <- geuvadis.spearmans.reformatted.df %>% select(gene_id, tissue, mesa_pop, geu_pop, ENunf_estimate, mashr_estimate, matrixeqtl_estimate, JTI_estimate, 'TIGAR_1e-04_estimate')
colnames(geuvadis.spearmans.reformatted.df) <- c('gene_id', 'tissue', 'mesa_pop', 'geu_pop', 'EN_estimate', 'mashr_estimate', 'matrixeqtl_estimate', 'JTI_estimate', 'TIGAR_estimate')

# pairwise wilcox test to see if models are significantly different
for (tis in c('PBMC','Mono','Tcell')){
  for (m_pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & m_pop=='CHN'){
      next
    } else {
      for (g_pop in c('GBR', 'YRI', 'ALL', 'CEU', 'FIN', 'TSI')){
        
        mashr_matrix_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, matrixeqtl_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(mashr_matrix_intersection)/2
        medians <- mashr_matrix_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=mashr_matrix_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'MASHR', 'MatrixeQTL', n_genes, medians$medians[1], medians$medians[2], pval)
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
        
        mashr_EN_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, EN_estimate) %>% drop_na() %>% melt() 
        n_genes <- nrow(mashr_EN_intersection)/2
        medians <- mashr_EN_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=mashr_EN_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'MASHR', 'EN', n_genes, medians$medians[1], medians$medians[2], pval)
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        mashr_TIGAR_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, TIGAR_estimate) %>% drop_na() %>% melt() 
        n_genes <- nrow(mashr_TIGAR_intersection)/2
        medians <- mashr_TIGAR_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=mashr_TIGAR_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'MASHR', 'TIGAR', n_genes, medians$medians[1], medians$medians[2], pval)
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        mashr_JTI_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(mashr_estimate, JTI_estimate) %>% drop_na() %>% melt() 
        n_genes <- nrow(mashr_JTI_intersection)/2
        medians <- mashr_JTI_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=mashr_JTI_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'MASHR', 'JTI', n_genes, medians$medians[1], medians$medians[2], pval)
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)

        EN_matrix_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(EN_estimate, matrixeqtl_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(EN_matrix_intersection)/2
        medians <- EN_matrix_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=EN_matrix_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'EN', 'MatrixeQTL', n_genes, medians$medians[1], medians$medians[2], pval)            
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        EN_TIGAR_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(EN_estimate, TIGAR_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(EN_TIGAR_intersection)/2
        medians <- EN_TIGAR_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=EN_TIGAR_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'EN', 'TIGAR', n_genes, medians$medians[1], medians$medians[2], pval)  
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        EN_JTI_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(EN_estimate, JTI_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(EN_JTI_intersection)/2
        medians <- EN_JTI_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=EN_JTI_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'EN', 'JTI', n_genes, medians$medians[1], medians$medians[2], pval)   
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        matrix_JTI_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(matrixeqtl_estimate, JTI_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(matrix_JTI_intersection)/2
        medians <- matrix_JTI_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=matrix_JTI_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'MatrixeQTL', 'JTI', n_genes, medians$medians[1], medians$medians[2], pval)          
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        matrix_TIGAR_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(matrixeqtl_estimate, TIGAR_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(matrix_TIGAR_intersection)/2
        medians <- matrix_TIGAR_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=matrix_TIGAR_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'MatrixeQTL', 'TIGAR', n_genes, medians$medians[1], medians$medians[2], pval)        
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        
        TIGAR_JTI_intersection <- geuvadis.spearmans.reformatted.df %>% filter(tissue==tis, mesa_pop==m_pop, geu_pop==g_pop) %>% 
          select(TIGAR_estimate, JTI_estimate) %>% drop_na() %>% melt()
        n_genes <- nrow(TIGAR_JTI_intersection)/2
        medians <- TIGAR_JTI_intersection %>% group_by(variable) %>% summarise(medians=median(value))
        test <- wilcox.test(value~variable, data=TIGAR_JTI_intersection)
        pval <- test$p.value
        new_line <- c(tis, m_pop, g_pop, 'TIGAR', 'JTI', n_genes, medians$medians[1], medians$medians[2], pval)  
        wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
      }
    }
  }
}
colnames(wilcox.tests.pvalues) <- c('Tissue','MESA_pop','GEU_pop', 'Method1', 'Method2', 'N_genes', 'Median_Spearman_Method1', 'Median_Spearman_Method2', 'Pvalue')
wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
wilcox.tests.pvalues$N_genes <- as.numeric(wilcox.tests.pvalues$N_genes)
wilcox.tests.pvalues$Median_Spearman_Method1 <- as.numeric(wilcox.tests.pvalues$Median_Spearman_Method1)
wilcox.tests.pvalues$Median_Spearman_Method2 <- as.numeric(wilcox.tests.pvalues$Median_Spearman_Method2)
wilcox.tests.pvalues$Pvalue <- as.numeric(wilcox.tests.pvalues$Pvalue)
#fwrite(wilcox.tests.pvalues, 'Geuvadis_wilcox_spearman_tests_pvalues_filteredbyh2.txt', col.names=T, sep='\t', quote=F)

### taking a look at only GBR and YRI
# filtering by geuvadis population
gbr.yri.spermans.df <- geuvadis.spearmans.df.filtered %>% filter(geu_pop=='GBR' | geu_pop=='YRI')
gbr.yri.spermans.df$model <- gsub('ENunf', 'EN', gbr.yri.spermans.df$model)
gbr.yri.spermans.df$model <- gsub('mashr', 'MASHR', gbr.yri.spermans.df$model)
gbr.yri.spermans.df$model <- gsub('matrixeqtl', 'MatrixeQTL', gbr.yri.spermans.df$model)
gbr.yri.spermans.df$model <- gsub('TIGAR_1e-04', 'TIGAR', gbr.yri.spermans.df$model)

# violin plots per tissue
for (tis in unique(gbr.yri.spermans.df$tissue)){
  # assessing genes that are shared between models and populations
  for (p in c('GBR', 'YRI')){
    for (m in c('AFA','EUR')){
      if (tis!='PBMC' & m=='CHN'){
        next
      } else {
        
        en.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='EN', geu_pop==p, mesa_pop==m)
        mashr.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='MASHR', geu_pop==p, mesa_pop==m)
        matrixeqtl.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='MatrixeQTL', geu_pop==p, mesa_pop==m)
        tigar.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='TIGAR', geu_pop==p, mesa_pop==m)
        jti.genes <- gbr.yri.spermans.df %>% filter(tissue==tis, model=='JTI', geu_pop==p, mesa_pop==m)
        
        full.shared.genes <- intersect(mashr.genes$gene_id, en.genes$gene_id) %>% intersect(jti.genes$gene_id) %>% 
          intersect(tigar.genes$gene_id) %>% intersect(matrixeqtl.genes$gene_id)
        en.mashr.jti.joint.df <- rbind(en.genes, mashr.genes, matrixeqtl.genes, tigar.genes, jti.genes) %>% filter(gene_id %in% full.shared.genes)
        
        if (exists('shared.genes.final.df')){
          shared.genes.final.df <- rbind(shared.genes.final.df, en.mashr.jti.joint.df)
        } else {shared.genes.final.df <- en.mashr.jti.joint.df}
      }
    }
  }
}

shared.genes.final.df %>% filter(tissue=='PBMC') %>% ggplot(., aes(x=geu_pop, y=estimate, fill=model)) + geom_violin() + facet_wrap(~mesa_pop) +
  geom_boxplot(width=0.4, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) +
  labs(title='Spearman correlation in Geuvadis', fill='Method') + ylab('Spearman correlation') + xlab('Geuvadis populations')
#ggsave('Geuvadis_spearman_correlation_intersection_all_methods_PBMC.pdf', height=4, width=8)
#ggsave('Geuvadis_spearman_correlation_intersection_all_methods_PBMC.tiff', height=4, width=8)  
#ggsave('Geuvadis_spearman_correlation_intersection_all_methods_PBMC.png', height=4, width=8)  
for (p in c('GBR', 'YRI')){
  for (m in c('AFA','EUR')){
    print(p %&% ' ' %&% m)
    shared.genes.final.df %>% filter(geu_pop==p, mesa_pop==m) %>% compare_means(estimate~model, ., method="wilcox.test", p.adjust.method='bonferroni') %>% print()
  }
}
