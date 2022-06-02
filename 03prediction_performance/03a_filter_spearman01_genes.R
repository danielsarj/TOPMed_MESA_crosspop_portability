setwd('/home/daniel/Geuvadis/WGS_predixcan')
library(data.table)
library(dplyr)
'%&%' = function(a,b) paste (a,b,sep='')

for (tissue in c('PBMC', 'Mono', 'Tcell')){
  for (model in c('EN', 'ENunf', 'mashr', 'matrixeqtl')){
    for (pop in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (geu in c('GBR', 'YRI')){
        if (tissue != 'PBMC' & pop == 'CHN'){
          next
        } else {
          table <- fread('Geuvadis.'%&% geu %&%'_'%&% tissue %&%'.'%&% pop %&%'_'%&% model %&%'_expression_spearman_correlation.txt') %>% 
            filter(estimate >= 0.1) %>% select(gene_id)
          table <- table[order(table$gene_id),]
          if (exists('final')){
            final <- rbind(final, table)
          } else {final <- table}
        }
      }
    }
    final <- final %>% unique()
    fwrite(final, tissue %&% '.'%&% model %&%'.GBRandYRI.genes.spearman01.txt', col.names = F, quote = F)
    rm(final)
  }
}
