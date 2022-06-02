library(data.table)
library(tidyverse)
library(viridis)
'%&%' = function(a,b) paste(a,b,sep='')
setwd('/home/daniel/Geuvadis/WGS_predixcan')

for (tissue in c('PBMC','Mono','Tcell')){
  if (tissue=='PBMC'){
    #loading geuvadis predictions
    GeuGBR_AFA_mashr<-fread('Geuvadis.GBR_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_mashr<-fread('Geuvadis.GBR_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_CHN_mashr<-fread('Geuvadis.GBR_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_mashr<-fread('Geuvadis.GBR_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_matrixeqtl<-fread('Geuvadis.GBR_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_matrixeqtl<-fread('Geuvadis.GBR_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_CHN_matrixeqtl<-fread('Geuvadis.GBR_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_matrixeqtl<-fread('Geuvadis.GBR_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_mashr<-fread('Geuvadis.YRI_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_mashr<-fread('Geuvadis.YRI_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_CHN_mashr<-fread('Geuvadis.YRI_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_mashr<-fread('Geuvadis.YRI_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_matrixeqtl<-fread('Geuvadis.YRI_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_matrixeqtl<-fread('Geuvadis.YRI_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_CHN_matrixeqtl<-fread('Geuvadis.YRI_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_matrixeqtl<-fread('Geuvadis.YRI_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_ENunf<-fread('Geuvadis.GBR_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_ENunf<-fread('Geuvadis.GBR_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_CHN_ENunf<-fread('Geuvadis.GBR_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_ENunf<-fread('Geuvadis.GBR_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_ENunf<-fread('Geuvadis.YRI_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_ENunf<-fread('Geuvadis.YRI_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_CHN_ENunf<-fread('Geuvadis.YRI_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_ENunf<-fread('Geuvadis.YRI_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    
    #binding rows by method
    mashr_spearmans<-bind_rows(GeuGBR_AFA_mashr, GeuGBR_EUR_mashr, GeuGBR_CHN_mashr, GeuGBR_HIS_mashr,
                               GeuYRI_AFA_mashr, GeuYRI_EUR_mashr, GeuYRI_CHN_mashr, GeuYRI_HIS_mashr)
    matrixeqtl_spearmans<-bind_rows(GeuGBR_AFA_matrixeqtl, GeuGBR_EUR_matrixeqtl, GeuGBR_CHN_matrixeqtl, GeuGBR_HIS_matrixeqtl,
                                    GeuYRI_AFA_matrixeqtl, GeuYRI_EUR_matrixeqtl, GeuYRI_CHN_matrixeqtl, GeuYRI_HIS_matrixeqtl)
    ENunf_spearmans<-bind_rows(GeuGBR_AFA_ENunf, GeuGBR_EUR_ENunf, GeuGBR_CHN_ENunf, GeuGBR_HIS_ENunf,
                               GeuYRI_AFA_ENunf, GeuYRI_EUR_ENunf, GeuYRI_CHN_ENunf, GeuYRI_HIS_ENunf)
    
    #joining everything in a single df
    joint_df <- full_join(matrixeqtl_spearmans, mashr_spearmans, by=c("gene_id","model_pop","Pop", "estimate", "model")) %>% full_join(ENunf_spearmans, by=c("gene_id","model_pop","Pop", "estimate", "model"))
    
    #violin plot w/o intersection of genes
    ggplot(data = joint_df, aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (all genes)')
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot.pdf', height=5, width=8)
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot.png', height=5, width=8)
    
    #getting the intersection across models per MESA pop
    for (pop in c('AFA','EUR','HIS','CHN')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      matrixeqtl_genes <- matrixeqtl_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% intersect(matrixeqtl_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else if (pop=='HIS'){
        his_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else {
        chn_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      }
    }
    
    #violinplot w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions, chn_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, 'PBMC_median_spearmans_GBR_YRI_intersection_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_joint_spearmans, aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes)')
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection.pdf', height=5, width=8)
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection.png', height=5, width=8)
    
    #wilcox rank sum test after intersection
    for (m in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (p in c('GBR', 'YRI')){
        matrix <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='ENunf') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = matrix)
        matrix <- test$p.value
        
        ENunf <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='MatrixeQTL') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        ENunf_matrix <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='Mashr') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf_matrix)
        ENunf_matrix <- test$p.value
        
        new_line <- c(m, p, matrix, ENunf, ENunf_matrix)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_matrix', 'mashr_ENunf', 'matrix_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_matrix <- as.numeric(wilcox.tests.pvalues$mashr_matrix)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    wilcox.tests.pvalues$matrix_ENunf <- as.numeric(wilcox.tests.pvalues$matrix_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_PBMC_GBR_YRI_pvalues_intersection.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    #getting the intersection across models per MESA pop - this time, excluding MatrixeQTL
    for (pop in c('AFA','EUR','HIS','CHN')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else if (pop=='HIS'){
        his_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else {
        chn_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      }
    }
    
    #violinplot w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions, chn_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, 'PBMC_median_spearmans_GBR_YRI_intersection_noMatrix_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_joint_spearmans,  aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes)') + stat_compare_means(method='wilcox.test', label='p.format') 
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection_noMatrix.pdf', height=7, width=8)
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection_noMatrix.png', height=7, width=8)
    
    #wilcox rank sum test after intersection w/o MatrixeQTL
    for (m in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (p in c('GBR', 'YRI')){
        ENunf <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='MatrixeQTL') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        new_line <- c(m, p, ENunf)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_PBMC_GBR_YRI_pvalues_intersection_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    #violinplot w/ intersection of genes with rho>0.1 in GBR or YRI in either EN and/or MASHR
    en_g <- fread('PBMC.ENunf.GBRandYRI.genes.spearman01.txt', header=F)
    mashr_g <- fread('PBMC.mashr.GBRandYRI.genes.spearman01.txt', header=F)
    en_or_mashr <- union(en_g, mashr_g) %>% pull()
    en_and_mashr <- intersect(en_g, mashr_g) %>% pull()
    
    intersect_joint_spearmans_and <- intersect_joint_spearmans %>% filter(gene_id %in% en_and_mashr)
    intersect_joint_spearmans_or <- intersect_joint_spearmans %>% filter(gene_id %in% en_or_mashr)
    
    intersect_median_df <- intersect_joint_spearmans_and %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, 'PBMC_median_spearmans_GBR_YRI_intersection_ENunf.and.MASHR_noMatrix_df.txt', col.names = T, sep = ' ')
    intersect_median_df <- intersect_joint_spearmans_or %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, 'PBMC_median_spearmans_GBR_YRI_intersection_ENunf.or.MASHR_noMatrix_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_joint_spearmans_and,  aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes, rho>0.1 in ENunf and MASHR)') + stat_compare_means(method='wilcox.test', label='p.format') 
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection_ENunf.and.MASHR_noMatrix.pdf', height=7, width=8)
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection_ENunf.and.MASHR_noMatrix.png', height=7, width=8)
    ggplot(data = intersect_joint_spearmans_or,  aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes, rho>0.1 in ENunf or MASHR)') + stat_compare_means(method='wilcox.test', label='p.format') 
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection_ENunf.or.MASHR_noMatrix.pdf', height=7, width=8)
    ggsave('Geuvadis_PBMC_GBR_YRI_violinplot_intersection_ENunf.or.MASHR_noMatrix.png', height=7, width=8)
    
    #wilcox rank sum test after intersection w/o MatrixeQTL
    for (m in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (p in c('GBR', 'YRI')){
        ENunf <- intersect_joint_spearmans_and %>% filter(grepl(m, model_pop), Pop == p) %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        new_line <- c(m, p, ENunf)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_PBMC_GBR_YRI_pvalues_intersection_ENunf.and.MASHR_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    for (m in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (p in c('GBR', 'YRI')){
        ENunf <- intersect_joint_spearmans_or %>% filter(grepl(m, model_pop), Pop == p) %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        new_line <- c(m, p, ENunf)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_PBMC_GBR_YRI_pvalues_intersection_ENunf.or.MASHR_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
  
  } else {
    GeuGBR_AFA_mashr<-fread('Geuvadis.GBR_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_mashr<-fread('Geuvadis.GBR_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_mashr<-fread('Geuvadis.GBR_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_matrixeqtl<-fread('Geuvadis.GBR_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_matrixeqtl<-fread('Geuvadis.GBR_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_matrixeqtl<-fread('Geuvadis.GBR_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_mashr<-fread('Geuvadis.YRI_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_mashr<-fread('Geuvadis.YRI_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_mashr<-fread('Geuvadis.YRI_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_matrixeqtl<-fread('Geuvadis.YRI_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_matrixeqtl<-fread('Geuvadis.YRI_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_matrixeqtl<-fread('Geuvadis.YRI_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_ENunf<-fread('Geuvadis.GBR_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_ENunf<-fread('Geuvadis.GBR_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_ENunf<-fread('Geuvadis.GBR_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_ENunf<-fread('Geuvadis.YRI_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_ENunf<-fread('Geuvadis.YRI_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_ENunf<-fread('Geuvadis.YRI_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    
    #binding rows by method
    mashr_spearmans<-bind_rows(GeuGBR_AFA_mashr, GeuGBR_EUR_mashr, GeuGBR_HIS_mashr,
                               GeuYRI_AFA_mashr, GeuYRI_EUR_mashr, GeuYRI_HIS_mashr)
    matrixeqtl_spearmans<-bind_rows(GeuGBR_AFA_matrixeqtl, GeuGBR_EUR_matrixeqtl, GeuGBR_HIS_matrixeqtl,
                                    GeuYRI_AFA_matrixeqtl, GeuYRI_EUR_matrixeqtl, GeuYRI_HIS_matrixeqtl)
    ENunf_spearmans<-bind_rows(GeuGBR_AFA_ENunf, GeuGBR_EUR_ENunf, GeuGBR_HIS_ENunf,
                               GeuYRI_AFA_ENunf, GeuYRI_EUR_ENunf, GeuYRI_HIS_ENunf)
    
    #joining everything in a single df
    joint_df <- full_join(matrixeqtl_spearmans, mashr_spearmans, by=c("gene_id","model_pop","Pop", "estimate", "model")) %>% full_join(ENunf_spearmans, by=c("gene_id","model_pop","Pop", "estimate", "model"))
    
    #violin plot w/o intersection of genes
    ggplot(data = joint_df, aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (all genes)')
    ggsave('Geuvadis_'%&% tissue %&%'_GBR_YRI_violinplot.pdf', height=4, width=8)
    ggsave('Geuvadis_'%&% tissue %&%'_GBR_YRI_violinplot.png', height=4, width=8)
    
    #getting the intersection across models per MESA pop
    for (pop in c('AFA','EUR','HIS')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      matrixeqtl_genes <- matrixeqtl_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% intersect(matrixeqtl_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else {
        his_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      }
    }
    
    #violinplot w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, tissue %&%'_median_spearmans_intersection_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_joint_spearmans, aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes)')
    ggsave('Geuvadis_'%&% tissue %&%'_GBR_YRI_violinplot_intersection.pdf', height=4, width=8)
    ggsave('Geuvadis_'%&% tissue %&%'_GBR_YRI_violinplot_intersection.png', height=4, width=8)
    
    #wilcox rank sum test after intersection
    for (m in c('AFA', 'EUR', 'HIS')){
      for (p in c('GBR', 'YRI')){
        matrix <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='ENunf') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = matrix)
        matrix <- test$p.value
        
        ENunf <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='MatrixeQTL') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        ENunf_matrix <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='Mashr') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf_matrix)
        ENunf_matrix <- test$p.value
        
        new_line <- c(m, p, matrix, ENunf, ENunf_matrix)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_matrix', 'mashr_ENunf', 'matrix_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_matrix <- as.numeric(wilcox.tests.pvalues$mashr_matrix)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    wilcox.tests.pvalues$matrix_ENunf <- as.numeric(wilcox.tests.pvalues$matrix_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_'%&% tissue %&%'_GBR_YRI_pvalues_intersection.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    #getting the intersection across models per MESA pop - this time, excluding MatrixeQTL
    for (pop in c('AFA','EUR','HIS')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else {
        his_predictions <- joint_df %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      }
    }
    
    #violinplot w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, tissue %&%'_median_spearmans_GBR_YRI_intersection_noMatrix_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_joint_spearmans,  aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes)') + stat_compare_means(method='wilcox.test', label='p.format') 
    ggsave('Geuvadis_'%&% tissue %&%'_GBR_YRI_violinplot_intersection_noMatrix.pdf', height=4, width=8)
    ggsave('Geuvadis_'%&% tissue %&%'_GBR_YRI_violinplot_intersection_noMatrix.png', height=4, width=8)
    
    #wilcox rank sum test after intersection w/o MatrixeQTL
    for (m in c('AFA', 'EUR', 'HIS')){
      for (p in c('GBR', 'YRI')){
        ENunf <- intersect_joint_spearmans %>% filter(grepl(m, model_pop), Pop == p, !model=='MatrixeQTL') %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        new_line <- c(m, p, ENunf)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_'%&% tissue %&%'_GBR_YRI_pvalues_intersection_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    #violinplot w/ intersection of genes with rho>0.1 in GBR or YRI in either EN and/or MASHR
    en_g <- fread(tissue %&% '.ENunf.GBRandYRI.genes.spearman01.txt', header=F)
    mashr_g <- fread(tissue %&% '.mashr.GBRandYRI.genes.spearman01.txt', header=F)
    en_or_mashr <- union(en_g, mashr_g) %>% pull()
    en_and_mashr <- intersect(en_g, mashr_g) %>% pull()
    
    intersect_joint_spearmans_and <- intersect_joint_spearmans %>% filter(gene_id %in% en_and_mashr)
    intersect_joint_spearmans_or <- intersect_joint_spearmans %>% filter(gene_id %in% en_or_mashr)
    
    intersect_median_df <- intersect_joint_spearmans_and %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, tissue %&% '_median_spearmans_GBR_YRI_intersection_ENunf.and.MASHR_noMatrix_df.txt', col.names = T, sep = ' ')
    intersect_median_df <- intersect_joint_spearmans_or %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, tissue %&% '_median_spearmans_GBR_YRI_intersection_ENunf.or.MASHR_noMatrix_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_joint_spearmans_and,  aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes, rho>0.1 in ENunf and MASHR)') + stat_compare_means(method='wilcox.test', label='p.format') 
    ggsave('Geuvadis_' %&% tissue %&% '_GBR_YRI_violinplot_intersection_ENunf.and.MASHR_noMatrix.pdf', height=7, width=8)
    ggsave('Geuvadis_' %&% tissue %&% '_GBR_YRI_violinplot_intersection_ENunf.and.MASHR_noMatrix.png', height=7, width=8)
    ggplot(data = intersect_joint_spearmans_or,  aes(x = Pop, y = estimate, fill = model)) + facet_wrap(~model_pop) +
      geom_violin() + geom_boxplot(width=0.1, position=position_dodge(0.9)) + scale_fill_viridis(discrete=T, alpha=0.5) + 
      labs(title='Spearman correlation in Geuvadis (intersection of genes, rho>0.1 in ENunf or MASHR)') + stat_compare_means(method='wilcox.test', label='p.format') 
    ggsave('Geuvadis_' %&% tissue %&% '_GBR_YRI_violinplot_intersection_ENunf.or.MASHR_noMatrix.pdf', height=7, width=8)
    ggsave('Geuvadis_' %&% tissue %&% '_GBR_YRI_violinplot_intersection_ENunf.or.MASHR_noMatrix.png', height=7, width=8)
    
    #wilcox rank sum test after intersection w/o MatrixeQTL
    for (m in c('AFA', 'EUR', 'HIS')){
      for (p in c('GBR', 'YRI')){
        ENunf <- intersect_joint_spearmans_and %>% filter(grepl(m, model_pop), Pop == p) %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        new_line <- c(m, p, ENunf)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_' %&% tissue %&% '_GBR_YRI_pvalues_intersection_ENunf.and.MASHR_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    for (m in c('AFA', 'EUR', 'HIS')){
      for (p in c('GBR', 'YRI')){
        ENunf <- intersect_joint_spearmans_or %>% filter(grepl(m, model_pop), Pop == p) %>% 
          select(model, estimate)
        test <- wilcox.test(estimate ~ model, data = ENunf)
        ENunf <- test$p.value
        
        new_line <- c(m, p, ENunf)
        
        if (exists('wilcox.tests.pvalues')){
          wilcox.tests.pvalues <- rbind(wilcox.tests.pvalues, new_line)
        } else {wilcox.tests.pvalues <- new_line}
      }
    }
    colnames(wilcox.tests.pvalues) <- c('model_pop','geu_pop', 'mashr_ENunf')
    wilcox.tests.pvalues <- as.data.frame(wilcox.tests.pvalues)
    wilcox.tests.pvalues$mashr_ENunf <- as.numeric(wilcox.tests.pvalues$mashr_ENunf)
    fwrite(wilcox.tests.pvalues, 'wilcox_' %&% tissue %&% '_GBR_YRI_pvalues_intersection_ENunf.or.MASHR_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues) 
    
  }
}