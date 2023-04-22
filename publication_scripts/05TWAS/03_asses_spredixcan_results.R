library(data.table)
library(tidyverse)
library(UpSetR)
library(viridis)
library(ggpubr)
library(janitor)
'%&%' = function(a,b) paste (a,b,sep='')
setwd('/home/daniel/SPrediXcan/WGS_SPredixcan')

page_phenos <- c('Height','QRS_duration','C-reactive_protein','BMI','Chronic_kidney','Smoking','Coffee','Diastolic_blood_pressure','Glomerular_filtration_rate','End-stage_renal_disease','Fasting_blood_glucose','Fasting_blood_insulin','Hemoglobin_A1c','HDL','LDL','Hypertension','Mean_corpuscular_hemoglobin','Platelet_count','PR_interval','QT_interval','Systolic_blood_pressure','Total_cholesterol','Triglycerides','TypeII_diabetes','WBC_count','Waist-hip_ratio-50','Waist-hip_ratio-51','Waist-hip_ratio-52')

# reading PAGE S-PrediXcan outputs
for (tis in c('PBMC','Mono','Tcell')){
  for (pop in c('AFA','EUR','HIS','CHN','ALL')){
    if (tis!='PBMC' & pop=='CHN'){
      next
    } else {
      if (pop=='ALL'){
        for (m in c('elasticnet_unfiltered')){
          for (pheno in page_phenos){
            spredixcan_output <- fread('PAGE/WojcikG_'%&% m %&%'_'%&% pheno %&%'_'%&% tis %&%'_'%&% pop %&%'.csv') %>% mutate(tissue=tis, population=pop, model=m, phenotype=pheno) 
          
            if (exists('compiled_page_spredixcan')){
              compiled_page_spredixcan <- rbind(compiled_page_spredixcan, spredixcan_output)
            } else {compiled_page_spredixcan <- spredixcan_output}
          }
        }
      } else {
        for (m in c('mashr','matrixeQTL', 'elasticnet_unfiltered')){
          for (pheno in page_phenos){
            spredixcan_output <- fread('PAGE/WojcikG_'%&% m %&%'_'%&% pheno %&%'_'%&% tis %&%'_'%&% pop %&%'.csv') %>% mutate(tissue=tis, population=pop, model=m, phenotype=pheno) 
            
            if (exists('compiled_page_spredixcan')){
              compiled_page_spredixcan <- rbind(compiled_page_spredixcan, spredixcan_output)
            } else {compiled_page_spredixcan <- spredixcan_output}
          }
        }
      }
    }
  }
}

panukbb_phenos <- c('BMI_calcuated','BMI_estimated','Coffee','C-reactive_protein','Diastolic_blood_pressure_auto','Diastolic_blood_pressure_manual','Fasting_blood_glucose','Glomerular_filtration_rate_cystain_C','Glomerular_filtration_rate_serum_creatinine_and_cystain_C','Glomerular_filtration_rate_serum_creatinine','HDL','Height_sitting','Height_standing','Hemoglobin_A1c','Hypertension_noncancer','Hypertension','LDL','Mean_corpuscular_hemoglobin','Platelet','PR_interval','QRS_duration','Smoking','Systolic_blood_pressure_auto','Systolic_blood_pressure_manual','Total_cholesterol','Triglycerides','TypeII_diabetes','Waist-hip_ratio_hip_circumference','Waist-hip_ratio_waist_circumference','WBC')

# reading PanUKBB S-PrediXcan outputs
for (tis in c('PBMC','Mono','Tcell')){
  for (pop in c('AFA','EUR','HIS','CHN','ALL')){
    if (tis!='PBMC' & pop=='CHN'){
      next
    } else {
      if (pop=='ALL'){
        for (m in c('elasticnet_unfiltered')){
          for (pheno in panukbb_phenos){
            if (pheno=='WBCHQ'){
              spredixcan_output <- fread('PanUKBB/Pan.UK.Biobank_meta_hq_'%&% m %&%'_'%&% pheno %&%'_'%&% tis %&%'_'%&% pop %&%'.csv') %>% mutate(tissue=tis, population=pop, model=m, phenotype=pheno) 
            } else{
              spredixcan_output <- fread('PanUKBB/Pan.UK.Biobank_meta_'%&% m %&%'_'%&% pheno %&%'_'%&% tis %&%'_'%&% pop %&%'.csv') %>% mutate(tissue=tis, population=pop, model=m, phenotype=pheno) 
            }
          
            if (exists('compiled_panukbb_spredixcan')){
              compiled_panukbb_spredixcan <- rbind(compiled_panukbb_spredixcan, spredixcan_output)
            } else {compiled_panukbb_spredixcan <- spredixcan_output}
          }
        }
      } else {
        for (m in c('mashr','matrixeQTL', 'elasticnet_unfiltered')){
          for (pheno in panukbb_phenos){
            if (pheno=='WBCHQ'){
              spredixcan_output <- fread('PanUKBB/Pan.UK.Biobank_meta_hq_'%&% m %&%'_'%&% pheno %&%'_'%&% tis %&%'_'%&% pop %&%'.csv') %>% mutate(tissue=tis, population=pop, model=m, phenotype=pheno) 
            } else{
              spredixcan_output <- fread('PanUKBB/Pan.UK.Biobank_meta_'%&% m %&%'_'%&% pheno %&%'_'%&% tis %&%'_'%&% pop %&%'.csv') %>% mutate(tissue=tis, population=pop, model=m, phenotype=pheno) 
            }
            
            if (exists('compiled_panukbb_spredixcan')){
              compiled_panukbb_spredixcan <- rbind(compiled_panukbb_spredixcan, spredixcan_output)
            } else {compiled_panukbb_spredixcan <- spredixcan_output}
          }
        }
      }
    }
  }
}

# renaming phenotypes so they match
compiled_panukbb_spredixcan$phenotype <- gsub('Height_sitting', 'Height', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Height_standing', 'Height', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('BMI_calcuated', 'BMI', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('BMI_estimated', 'BMI', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Diastolic_blood_pressure_auto', 'Diastolic_blood_pressure', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Diastolic_blood_pressure_manual', 'Diastolic_blood_pressure', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Glomerular_filtration_rate_cystain_C', 'Glomerular_filtration_rate', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Glomerular_filtration_rate_serum_creatinine_and_cystain_C', 'Glomerular_filtration_rate', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Glomerular_filtration_rate_serum_creatinine', 'Glomerular_filtration_rate', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Hypertension_noncancer', 'Hypertension', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Platelet', 'Platelet_count', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Systolic_blood_pressure_auto', 'Systolic_blood_pressure', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Systolic_blood_pressure_manual', 'Systolic_blood_pressure', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('WBC', 'WBC_count', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Waist-hip_ratio_hip_circumference', 'Waist-hip_ratio', compiled_panukbb_spredixcan$phenotype)
compiled_panukbb_spredixcan$phenotype <- gsub('Waist-hip_ratio_waist_circumference', 'Waist-hip_ratio', compiled_panukbb_spredixcan$phenotype)
compiled_page_spredixcan$phenotype <- gsub('Waist-hip_ratio-50', 'Waist-hip_ratio', compiled_page_spredixcan$phenotype)
compiled_page_spredixcan$phenotype <- gsub('Waist-hip_ratio-51', 'Waist-hip_ratio', compiled_page_spredixcan$phenotype)
compiled_page_spredixcan$phenotype <- gsub('Waist-hip_ratio-52', 'Waist-hip_ratio', compiled_page_spredixcan$phenotype)

# merging both spredixcan results
compiled_page_spredixcan <- compiled_page_spredixcan %>% mutate(study='PAGE')
compiled_panukbb_spredixcan <- compiled_panukbb_spredixcan %>% mutate(study='PanUKBB')
compiled_both_studies <- rbind(compiled_page_spredixcan, compiled_panukbb_spredixcan)

# reading significant non-zero heritability estimates
h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01)

# filtering spredixcan outputs so it only contains genes in the h2 df
for (tis in c('PBMC','Mono','Tcell')){
  genes.in.tissue <- h2estimates %>% filter(tissue==tis) %>% select(gene) %>% unique()
  page.filt <- compiled_page_spredixcan %>% filter(tissue==tis, pvalue < 5e-8, gene %in% genes.in.tissue$gene)
  panukbb.filt <- compiled_panukbb_spredixcan %>% filter(tissue==tis, pvalue < 5e-8, gene %in% genes.in.tissue$gene)
  
  if (exists('page_compiled_filtered')){
    page_compiled_filtered <- rbind(page_compiled_filtered, page.filt)
  } else {page_compiled_filtered <- page.filt}
  
  if (exists('panukbb_compiled_filtered')){
    panukbb_compiled_filtered <- rbind(panukbb_compiled_filtered, panukbb.filt)
  } else {panukbb_compiled_filtered <- panukbb.filt}
}

# joining PAGE and PanUKBB based on gene, phenotype and model
reduced_page <- page_compiled_filtered %>% select(gene, gene_name, tissue, population, model, phenotype, zscore, effect_size, pvalue) %>% rename(page_zscore=zscore, page_pvalue=pvalue, page_effectsize=effect_size)
reduced_panukbb <- panukbb_compiled_filtered %>% select(gene, gene_name, tissue, population, model, phenotype, zscore, effect_size, pvalue) %>% rename(panukbb_zscore=zscore, panukbb_pvalue=pvalue, panukbb_effectsize=effect_size)
joint_df <- inner_join(reduced_page, reduced_panukbb, by=c('gene', 'gene_name', 'tissue', 'population', 'model', 'phenotype'))
joint_df$model <- gsub('mashr', 'MASHR', joint_df$model)
joint_df$model <- gsub('matrixeQTL', 'MatrixeQTL', joint_df$model)
joint_df$model <- gsub('elasticnet_unfiltered', 'EN', joint_df$model)
fwrite(joint_df, 'PAGE_PanUKBB_significant_results.txt', col.names=T, quote=F, sep='\t')

# filter based on same direction of effect
joint_df <- joint_df %>% filter((page_effectsize > 0 & panukbb_effectsize > 0) | (page_effectsize < 0 & panukbb_effectsize < 0))

# correlation of zscores in PAGE and PanUKBB
ggplot(joint_df, aes(x=page_zscore, y=panukbb_zscore)) + geom_point() + geom_density2d() +
  geom_smooth(method='lm', color='red') + stat_cor(method='pearson') + facet_wrap(~model)

# summarizing results, counting gene-trait pairs
summ <- joint_df %>% group_by(tissue, population, model, phenotype) %>% summarise(n=n())
summ$model <- gsub('mashr', 'MASHR', summ$model)
summ$model <- gsub('matrixeQTL', 'MatrixeQTL', summ$model)
summ$model <- gsub('elasticnet_unfiltered', 'EN', summ$model)

summ %>% filter(tissue=='PBMC') %>% ggplot(., aes(x=phenotype, y=n, fill=model)) + geom_col(position='dodge') + 
  coord_flip() +  scale_fill_viridis_d(name='Model') + facet_wrap(~population) + xlab('Phenotypes') + ylab('Number of significant genes')
ggsave('PAGE_PanUKBB_PBMC_colgraph_filth2_sig_samedirection.png', height=6, width=8)
ggsave('PAGE_PanUKBB_PBMC_colgraph_filth2_sig_samedirection.pdf', height=6, width=8)
ggsave('PAGE_PanUKBB_PBMC_colgraph_filth2_sig_samedirection.tiff', height=6, width=8)

summ %>% ggplot(., aes(x=phenotype, y=n, fill=model)) + geom_col(position='dodge') + 
  coord_flip() +  scale_fill_viridis_d(name='Model') + facet_wrap(~population) + xlab('Phenotypes') + ylab('Number of significant genes')
ggsave('PAGE_PanUKBB_colgraph_filth2_sig_samedirection.png', height=6, width=8)
ggsave('PAGE_PanUKBB_colgraph_filth2_sig_samedirection.pdf', height=6, width=8)
ggsave('PAGE_PanUKBB_colgraph_filth2_sig_samedirection.tiff', height=6, width=8)

# assessing gene-trait pairs found by only one method
mashr_pairs <- joint_df %>% filter(model=='MASHR') %>% select(gene, phenotype) %>% unique()
matrixeqtl_pairs <- joint_df %>% filter(model=='MatrixeQTL') %>% select(gene, phenotype) %>% unique()
en_pairs <- joint_df %>% filter(model=='EN', population!='ALL') %>% select(gene, phenotype) %>% unique()

unique_mashr <- anti_join(mashr_pairs, matrixeqtl_pairs) %>% anti_join(en_pairs) %>% summarise(n_gene_pairs=n()) %>% mutate(model='MASHR')
unique_matrixeqtl <- anti_join(matrixeqtl_pairs, mashr_pairs) %>% anti_join(en_pairs) %>% summarise(n_gene_pairs=n()) %>% mutate(model='MatrixeQTL')
unique_en <- anti_join(en_pairs, matrixeqtl_pairs) %>% anti_join(mashr_pairs) %>% summarise(n_gene_pairs=n()) %>% mutate(model='EN')

joint_df_of_uniques <- rbind(unique_en, unique_mashr, unique_matrixeqtl)
unq_col <- ggplot(joint_df_of_uniques, aes(x=model, y=n_gene_pairs, fill=model)) + geom_col(position='dodge') +
  scale_fill_viridis_d() + xlab('MESA population') + ylab('Number of unique gene-trait pairs') + labs(fill='Model')


## total discoveries
n_pairs_per_pop_tissue <- joint_df %>% filter(population!='ALL') %>% group_by(model, tissue, population) %>% select(gene, gene_name, phenotype) %>% unique() %>% summarise(n_pairs=n())
n_pairs_per_pop <- joint_df %>% filter(population!='ALL') %>% group_by(model, population) %>% select(gene, gene_name, phenotype) %>% unique() %>% summarise(n_pairs=n())
n_pairs_per_pop$model <- gsub('mashr', 'MASHR', n_pairs_per_pop$model)
n_pairs_per_pop$model <- gsub('matrixeQTL', 'MatrixeQTL', n_pairs_per_pop$model)
n_pairs_per_pop$model <- gsub('elasticnet_unfiltered', 'EN', n_pairs_per_pop$model)

total_col <- ggplot(n_pairs_per_pop, aes(x=population, y=n_pairs, fill=model)) + geom_col(position='dodge') +
  scale_fill_viridis_d() + xlab('MESA population') + ylab('Number of gene-trait pairs') + labs(fill='Model')

ggarrange(total_col, unq_col, ncol=1, legend='right', common.legend=T)
ggsave('PAGE_PanUKBB_n_sig_pairs_colgraph.pdf', height=6, width=8)
ggsave('PAGE_PanUKBB_n_sig_pairs_colgraph.png', height=6, width=8)
ggsave('PAGE_PanUKBB_n_sig_pairs_colgraph.tiff', height=6, width=8)

## adding unique/non-unique info to dataframe
unique_mashr <- anti_join(mashr_pairs, matrixeqtl_pairs) %>% anti_join(en_pairs) %>% mutate(model='MASHR', unique='Y')
unique_matrixeqtl <- anti_join(matrixeqtl_pairs, mashr_pairs) %>% anti_join(en_pairs) %>% mutate(model='MatrixeQTL', unique='Y')
unique_en <- anti_join(en_pairs, matrixeqtl_pairs) %>% anti_join(mashr_pairs) %>% mutate(model='EN', unique='Y')
updated_joint_df <- rbind(unique_en, unique_mashr, unique_matrixeqtl) %>% full_join(joint_df) %>% 
  select(gene, gene_name, tissue, population, model, phenotype, page_zscore, page_pvalue, page_effectsize, panukbb_zscore, panukbb_pvalue, panukbb_effectsize, unique) %>%
  arrange(tissue, population, model, phenotype)
updated_joint_df$unique <- if_else(is.na(updated_joint_df$unique), 'N', 'Y')
fwrite(updated_joint_df, 'PAGE_PanUKBB_significant_results_updt.txt', col.names=T, quote=F, sep='\t')

## adding GWAS Catalog studies to dataframe
gwas_catalog <- fread('../alternative.gz') %>% clean_names() %>% select(disease_trait, mapped_trait, reported_gene_s, mapped_gene, study_accession)
genes_in_spredxican <- updated_joint_df %>% filter(population!='ALL') %>% select(gene_name, phenotype) %>% unique()

for (g in unique(genes_in_spredxican$gene_name)){
  tmp1 <- gwas_catalog %>% filter(grepl(g, reported_gene_s)) %>% select(disease_trait, mapped_trait, reported_gene_s, study_accession) %>% rename(genes_reported=reported_gene_s)
  tmp2 <- gwas_catalog %>% filter(grepl(g, mapped_gene)) %>% select(disease_trait, mapped_trait, mapped_gene, study_accession) %>% rename(genes_reported=mapped_gene)
  tmp_concat <- rbind(tmp1, tmp2) %>% unique()
  
  if(exists('gwas_catalog_filtered')){
    gwas_catalog_filtered <- rbind(gwas_catalog_filtered, tmp_concat)
  } else { gwas_catalog_filtered <- tmp_concat }
}

for (i in c(1:nrow(genes_in_spredxican))){
  gwas_subset <- gwas_catalog_filtered %>% filter(grepl(genes_in_spredxican$gene_name[i], genes_reported))
  
  if (genes_in_spredxican$phenotype[i]=='C-reactive_protein'){
    gwas_subset <- gwas_subset %>% filter(grepl('C-reactive protein', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Platelet_count'){
    gwas_subset <- gwas_subset %>% filter(grepl('platelet count', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='WBC_count'){
    gwas_subset <- gwas_subset %>% filter(grepl('monocyte count|lymphocyte count|neutrophil count|eosinophil count|basophil count', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Fasting_blood_glucose'){
    gwas_subset <- gwas_subset %>% filter(grepl('fasting blood glucose', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='HDL'){
    gwas_subset <- gwas_subset %>% filter(grepl('HDL cholesterol levels|high density lipoprotein cholesterol measurement', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='LDL'){
    gwas_subset <- gwas_subset %>% filter(grepl('LDL cholesterol levels|low density lipoprotein cholesterol measurement', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Height'){
    gwas_subset <- gwas_subset %>% filter(grepl('body height', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Hemoglobin_A1c'){
    gwas_subset <- gwas_subset %>% filter(grepl('HbA1c measurement', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Total_cholesterol'){
    gwas_subset <- gwas_subset %>% filter(grepl('total cholesterol', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Triglycerides'){
    gwas_subset <- gwas_subset %>% filter(grepl('triglyceride measurement', mapped_trait)) %>% select(study_accession) %>% unique()
  } else if (genes_in_spredxican$phenotype[i]=='Mean_corpuscular_hemoglobin'){
    gwas_subset <- gwas_subset %>% filter(grepl('mean corpuscular hemoglobin', mapped_trait)) %>% select(study_accession) %>% unique()
  } 
  
  if (nrow(gwas_subset)>0){
    gwas_subset <- gwas_subset$study_accession %>% t() %>% as.data.frame() 
    gwas_subset <- gwas_subset %>% unite('GWAS_catalog', c(colnames(gwas_subset)), sep=', ')

    small_updated_joint_df <- updated_joint_df %>% filter(gene_name==genes_in_spredxican$gene_name[i], phenotype==genes_in_spredxican$phenotype[i])
    small_updated_joint_df <- small_updated_joint_df %>% mutate(GWAS_catalog=gwas_subset) %>% as.data.table()
    
    if (exists('compiled_joint_df_wGWAScatalogstudies')){
      compiled_joint_df_wGWAScatalogstudies <- rbind(compiled_joint_df_wGWAScatalogstudies, small_updated_joint_df) %>% as.data.table()
    } else { compiled_joint_df_wGWAScatalogstudies <- small_updated_joint_df }
  }
}

updated_joint_df <- left_join(updated_joint_df, compiled_joint_df_wGWAScatalogstudies) %>% filter(population!='ALL')
updated_joint_df$GWAS_catalog <- updated_joint_df$GWAS_catalog %>% replace_na('none')
fwrite(updated_joint_df, 'PAGE_PanUKBB_significant_results_updt_wGWASCatalog.txt', col.names=T, quote=F, sep='\t')
fwrite(updated_joint_df, 'PAGE_PanUKBB_significant_results_updt_wGWASCatalog.csv', col.names=T, quote=F, sep=',')

novel_assocs <- updated_joint_df %>% filter(GWAS_catalog=='none')
sum_novel_assocs <- novel_assocs %>% select(gene_name, population, model, phenotype) %>% unique() %>% group_by(population, model) %>% summarise(novel_associations=n())
ggplot(sum_novel_assocs, aes(x=population, y=novel_associations, fill=model)) + geom_col(position='dodge') + 
  scale_fill_viridis_d(name='Model') + xlab('Populations') + ylab('Number of possible novel genes-trait pairs')
ggsave('PAGE_PanUKBB_new_assocs_per_model_colgraph.png', height=4, width=6)
ggsave('PAGE_PanUKBB_new_assocs_per_model_colgraph.pdf', height=4, width=6)
ggsave('PAGE_PanUKBB_new_assocs_per_model_colgraph.tiff', height=4, width=6)

reported_assocs <- updated_joint_df %>% filter(GWAS_catalog!='none') %>% select(gene_name, phenotype) %>% unique()
reported_assocs_permodel <- updated_joint_df %>% filter(GWAS_catalog!='none') %>% select(gene_name, population, model, phenotype) %>% 
  unique() %>% group_by(population, model) %>% summarise(n_pairs=n())
ggplot(reported_assocs_permodel, aes(x=population, y=n_pairs, fill=model)) + geom_col(position='dodge') + 
  scale_fill_viridis_d(name='Model') + xlab('Populations') + ylab('Number of known GWAS catalog associations')
ggsave('PAGE_PanUKBB_known_assocs_per_model_colgraph.png', height=4, width=6)
ggsave('PAGE_PanUKBB_known_assocs_per_model_colgraph.pdf', height=4, width=6)
ggsave('PAGE_PanUKBB_known_assocs_per_model_colgraph.tiff', height=4, width=6)

# upset plot
pop_traits <- updated_joint_df %>% select(gene, phenotype) %>% unique()
mashr_pop_upset_df <- data.frame(MASHR_AFA=as.numeric(), MASHR_EUR=as.numeric(), MASHR_CHN=as.numeric(), MASHR_HIS=as.numeric())
matrixeqtl_pop_upset_df <- data.frame(MatrixeQTL_AFA=as.numeric(), MatrixeQTL_EUR=as.numeric(), MatrixeQTL_CHN=as.numeric(), MatrixeQTL_HIS=as.numeric())
en_pop_upset_df <- data.frame(EN_AFA=as.numeric(), EN_EUR=as.numeric(), EN_CHN=as.numeric(), EN_HIS=as.numeric())
for (i in 1:nrow(pop_traits)){
  mashr_tmp_df <- data.frame(MASHR_AFA=0, MASHR_EUR=0, MASHR_CHN=0, MASHR_HIS=0)
  matrixeqtl_tmp_df <- data.frame(MatrixeQTL_AFA=0, MatrixeQTL_EUR=0, MatrixeQTL_CHN=0, MatrixeQTL_HIS=0)
  en_tmp_df <- data.frame(EN_AFA=0, EN_EUR=0, EN_CHN=0, EN_HIS=0)
  
  subset_joint_df <- updated_joint_df %>% filter(gene==pop_traits$gene[i], phenotype==pop_traits$phenotype[i])
  
  if (nrow(filter(subset_joint_df, model=='MASHR', population=='AFA'))>0){ 
    mashr_afa <- 1 
  } else { mashr_afa <- 0}
  if (nrow(filter(subset_joint_df, model=='MASHR', population=='EUR'))>0){ 
    mashr_eur <- 1 
  } else { mashr_eur <- 0}
  if (nrow(filter(subset_joint_df, model=='MASHR', population=='CHN'))>0){ 
    mashr_chn <- 1 
  } else { mashr_chn <- 0}
  if (nrow(filter(subset_joint_df, model=='MASHR', population=='HIS'))>0){ 
    mashr_his <- 1 
  } else { mashr_his <- 0}
  
  if (nrow(filter(subset_joint_df, model=='MatrixeQTL', population=='AFA'))>0){ 
    matrixeqtl_afa <- 1 
  } else { matrixeqtl_afa <- 0}
  if (nrow(filter(subset_joint_df, model=='MatrixeQTL', population=='EUR'))>0){ 
    matrixeqtl_eur <- 1 
  } else { matrixeqtl_eur <- 0}
  if (nrow(filter(subset_joint_df, model=='MatrixeQTL', population=='CHN'))>0){ 
    matrixeqtl_chn <- 1 
  } else { matrixeqtl_chn <- 0}
  if (nrow(filter(subset_joint_df, model=='MatrixeQTL', population=='HIS'))>0){ 
    matrixeqtl_his <- 1 
  } else { matrixeqtl_his <- 0}
  
  if (nrow(filter(subset_joint_df, model=='EN', population=='AFA'))>0){ 
    en_afa <- 1 
  } else { en_afa <- 0}
  if (nrow(filter(subset_joint_df, model=='EN', population=='EUR'))>0){ 
    en_eur <- 1 
  } else { en_eur <- 0}
  if (nrow(filter(subset_joint_df, model=='EN', population=='CHN'))>0){ 
    en_chn <- 1 
  } else { en_chn <- 0}
  if (nrow(filter(subset_joint_df, model=='EN', population=='HIS'))>0){ 
    en_his <- 1 
  } else { en_his <- 0}
  
  mashr_tmp_df[1,1] <- mashr_afa
  mashr_tmp_df[1,2] <- mashr_eur
  mashr_tmp_df[1,3] <- mashr_chn
  mashr_tmp_df[1,4] <- mashr_his
  
  matrixeqtl_tmp_df[1,1] <- matrixeqtl_afa
  matrixeqtl_tmp_df[1,2] <- matrixeqtl_eur
  matrixeqtl_tmp_df[1,3] <- matrixeqtl_chn
  matrixeqtl_tmp_df[1,4] <- matrixeqtl_his
  
  en_tmp_df[1,1] <- en_afa
  en_tmp_df[1,2] <- en_eur
  en_tmp_df[1,3] <- en_chn
  en_tmp_df[1,4] <- en_his
  
  mashr_pop_upset_df <- rbind(mashr_pop_upset_df, mashr_tmp_df)
  matrixeqtl_pop_upset_df <- rbind(matrixeqtl_pop_upset_df, matrixeqtl_tmp_df)
  en_pop_upset_df <- rbind(en_pop_upset_df, en_tmp_df)
  
}
pdf('mashr_pop_upset_plot.pdf', width=8, height=6.5)
mashr_pop_upset_df %>% upset(point.size=5.5, line.size=3, text.scale=3, nintersects=NA, keep.order=T, sets=c('MASHR_AFA','MASHR_EUR','MASHR_HIS','MASHR_CHN'), empty.intersections='on', mainbar.y.max=20)
dev.off() 
png('mashr_pop_upset_plot.png', width=8, height=6.5, units='in', res=300)
mashr_pop_upset_df %>% upset(point.size=5.5, line.size=3, text.scale=3, nintersects=NA, keep.order=T, sets=c('MASHR_AFA','MASHR_EUR','MASHR_HIS','MASHR_CHN'), empty.intersections='on', mainbar.y.max=20)
dev.off() 
tiff('mashr_pop_upset_plot.tiff', width=8, height=6.5, units='in', res=300)
mashr_pop_upset_df %>% upset(point.size=5.5, line.size=3, text.scale=3, nintersects=NA, keep.order=T, sets=c('MASHR_AFA','MASHR_EUR','MASHR_HIS','MASHR_CHN'), empty.intersections='on', mainbar.y.max=20)
dev.off() 

pdf('en_pop_upset_plot.pdf', width=8, height=6.5)
en_pop_upset_df %>% upset(point.size=5.5, line.size=3, text.scale=3, nintersects=NA, keep.order=T, sets=c('EN_AFA','EN_EUR','EN_HIS','EN_CHN'), empty.intersections='on', mainbar.y.max=20)
dev.off() 
png('en_pop_upset_plot.png', width=8, height=6.5, units='in', res=300)
en_pop_upset_df %>% upset(point.size=5.5, line.size=3, text.scale=3, nintersects=NA, keep.order=T, sets=c('EN_AFA','EN_EUR','EN_HIS','EN_CHN'), empty.intersections='on', mainbar.y.max=20)
dev.off() 
tiff('en_pop_upset_plot.tiff', width=8, height=6.5, units='in', res=300)
en_pop_upset_df %>% upset(point.size=5.5, line.size=3, text.scale=3, nintersects=NA, keep.order=T, sets=c('EN_AFA','EN_EUR','EN_HIS','EN_CHN'), empty.intersections='on', mainbar.y.max=20)
dev.off() 
