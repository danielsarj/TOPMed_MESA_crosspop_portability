library(data.table)
library(janitor)
library(tidyverse)
'%&%' = function(a,b) paste (a,b,sep='')
setwd('/home/chris/topmed_expression_whole_genome/expression/raw_expression')

exp_tbl <- fread('TOPMed_MESA_RNAseq_Pilot_RSEMv1.3.0.rsem_genes_tpm.txt') %>% select(-'transcript_id(s)') %>% 
  separate('gene_id', into=c('gene_id', 'right'), sep='\\.') %>% select(-right) %>% filter(grepl('ENSG', gene_id))

for (t in c('PBMC', 'Mono', 'Tcell')){
  heritable_genes <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% 
    filter(h2-2*se > 0.01, tissue==t) %>% select(gene) %>% unique() %>% pull()
  
  afa_ids <- fread('/home/daniel/MESA_genotypes_subset/IDs_files/' %&% t %&% '_ALL_fulldata.txt') %>% filter(pop=='AFA') %>% select(torid) %>% pull()
  eur_ids <- fread('/home/daniel/MESA_genotypes_subset/IDs_files/' %&% t %&% '_ALL_fulldata.txt') %>% filter(pop=='EUR') %>% select(torid) %>% pull()
  his_ids <- fread('/home/daniel/MESA_genotypes_subset/IDs_files/' %&% t %&% '_ALL_fulldata.txt') %>% filter(pop=='HIS') %>% select(torid) %>% pull()
  
  afa_exp <- exp_tbl %>% select(gene_id, all_of(afa_ids)) %>% rowwise() %>% mutate(AFA = median(c_across(where(is.numeric)), na.rm=T)) %>% 
    filter(AFA > 0.1) %>% select(gene_id, AFA) 
  eur_exp <- exp_tbl %>% select(gene_id, all_of(eur_ids)) %>% rowwise() %>% mutate(EUR = median(c_across(where(is.numeric)), na.rm=T)) %>% 
    filter(EUR > 0.1) %>% select(gene_id, EUR) 
  his_exp <- exp_tbl %>% select(gene_id, all_of(his_ids)) %>% rowwise() %>% mutate(HIS = median(c_across(where(is.numeric)), na.rm=T)) %>% 
    filter(HIS > 0.1) %>% select(gene_id, HIS) 
    
  if (t == 'PBMC'){
    chn_ids <- fread('/home/daniel/MESA_genotypes_subset/IDs_files/PBMC_ALL_fulldata.txt') %>% filter(pop=='CHN') %>% select(torid) %>% pull()
    
    chn_exp <- exp_tbl %>% select(gene_id, all_of(chn_ids)) %>% rowwise() %>% mutate(CHN = median(c_across(where(is.numeric)), na.rm=T)) %>% 
      filter(CHN > 0.1) %>% select(gene_id, CHN) 
    
    corr_exp <- full_join(afa_exp, eur_exp, by=c('gene_id')) %>% full_join(his_exp, by=c('gene_id')) %>%
      full_join(chn_exp, by=c('gene_id')) %>% drop_na() %>% select(-gene_id) %>% cor(method='pearson') %>% as.data.frame()
    
    fwrite(corr_exp, '/home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/PBMC_transcriptome_wide_median_correlation.txt', col.names=T, row.names=T, quote=F, sep='\t')
    
    corr_exp <- full_join(afa_exp, eur_exp, by=c('gene_id')) %>% full_join(his_exp, by=c('gene_id')) %>%
      full_join(chn_exp, by=c('gene_id')) %>% drop_na() %>% filter(gene_id %in% heritable_genes) %>% select(-gene_id) %>% cor(method='pearson') %>% as.data.frame()
    
    fwrite(corr_exp, '/home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/PBMC_transcriptome_wide_median_correlation_h2filt.txt', col.names=T, row.names=T, quote=F, sep='\t')
    
  } else {
    corr_exp <- full_join(afa_exp, eur_exp, by=c('gene_id')) %>% full_join(his_exp, by=c('gene_id')) %>%
      drop_na() %>% select(-gene_id) %>% cor(method='pearson') %>% as.data.frame()
    
    fwrite(corr_exp, '/home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/' %&% t %&% '_transcriptome_wide_median_correlation.txt', col.names=T, row.names=T, quote=F, sep='\t')
    
    corr_exp <- full_join(afa_exp, eur_exp, by=c('gene_id')) %>% full_join(his_exp, by=c('gene_id')) %>%
      drop_na() %>% filter(gene_id %in% heritable_genes) %>% select(-gene_id) %>% cor(method='pearson') %>% as.data.frame()
    
    fwrite(corr_exp, '/home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/' %&% t %&% '_transcriptome_wide_median_correlation_h2filt.txt', col.names=T, row.names=T, quote=F, sep='\t')
  }
}