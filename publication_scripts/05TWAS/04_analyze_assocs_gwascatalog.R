library(tidyverse)
library(data.table)

main_df <- fread('PAGE_PanUKBB_significant_results_updt_wGWASCatalog.txt')
all_assocs <- main_df %>% select(gene_name, phenotype) %>% unique()
non_reported_assocs <- main_df %>% filter(GWAS_catalog=='none') %>% select(gene_name, phenotype) %>% unique()

non_reported_assocs_permodel <- main_df %>% filter(GWAS_catalog=='none') %>% select(gene_name, population, model, phenotype) %>% 
  unique() %>% group_by(population, model) %>% summarise(n_pairs=n())

reported_assocs <- main_df %>% filter(GWAS_catalog!='none') %>% select(gene_name, phenotype) %>% unique()
reported_assocs_permodel <- main_df %>% filter(GWAS_catalog!='none') %>% select(gene_name, population, model, phenotype) %>% 
  unique() %>% group_by(population, model) %>% summarise(n_pairs=n())

reported_assocs_and_unique <- main_df %>% filter(GWAS_catalog!='none', unique=='Y') %>% select(gene_name, phenotype, model) %>% unique()
table(reported_assocs_and_unique$model)

reported_assocs_and_unique_perpop <- main_df %>% filter(GWAS_catalog!='none', unique=='Y') %>% select(gene_name, phenotype, model, population) %>% unique()
assocs_p_m <- reported_assocs_and_unique_perpop %>% group_by(gene_name, phenotype, model) %>% summarise(n=n())
table(reported_assocs_and_unique_perpop$model)

non_reported_assocs_and_unique <- main_df %>% filter(GWAS_catalog=='none', unique=='Y') %>% select(gene_name, phenotype, model) %>% unique()
table(non_reported_assocs_and_unique$model)


### 

a <- main_df %>% filter(gene_name=='BUD13', phenotype=='HDL') %>% select(GWAS_catalog) %>% unique()
a <- a$GWAS_catalog

