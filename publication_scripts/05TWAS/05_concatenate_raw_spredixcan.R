library(data.table)
library(tidyverse)
library(janitor)
'%&%' = function(a,b) paste (a,b,sep='')
setwd('/home/daniel/SPrediXcan/WGS_SPredixcan')

page_phenos <- c('Height','QRS_duration','C-reactive_protein','BMI','Chronic_kidney','Smoking','Coffee','Diastolic_blood_pressure','Glomerular_filtration_rate','End-stage_renal_disease','Fasting_blood_glucose','Fasting_blood_insulin','Hemoglobin_A1c','HDL','LDL','Hypertension','Mean_corpuscular_hemoglobin','Platelet_count','PR_interval','QT_interval','Systolic_blood_pressure','Total_cholesterol','Triglycerides','TypeII_diabetes','WBC_count','Waist-hip_ratio-50','Waist-hip_ratio-51','Waist-hip_ratio-52')
# reading PAGE S-PrediXcan outputs
for (tis in c('PBMC','Mono','Tcell')){
  for (pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & pop=='CHN'){
      next
    } else {
      for (m in c('mashr','matrixeQTL', 'elasticnet_unfiltered', 'TIGAR', 'JTI')){
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

panukbb_phenos <- c('BMI_calcuated','BMI_estimated','Coffee','C-reactive_protein','Diastolic_blood_pressure_auto','Diastolic_blood_pressure_manual','Fasting_blood_glucose','Glomerular_filtration_rate_cystain_C','Glomerular_filtration_rate_serum_creatinine_and_cystain_C','Glomerular_filtration_rate_serum_creatinine','HDL','Height_sitting','Height_standing','Hemoglobin_A1c','Hypertension_noncancer','Hypertension','LDL','Mean_corpuscular_hemoglobin','Platelet','PR_interval','QRS_duration','Smoking','Systolic_blood_pressure_auto','Systolic_blood_pressure_manual','Total_cholesterol','Triglycerides','TypeII_diabetes','Waist-hip_ratio_hip_circumference','Waist-hip_ratio_waist_circumference','WBC')
# reading PanUKBB S-PrediXcan outputs
for (tis in c('PBMC','Mono','Tcell')){
  for (pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & pop=='CHN'){
      next
    } else {
      for (m in c('mashr','matrixeQTL', 'elasticnet_unfiltered', 'TIGAR', 'JTI')){
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

# merging both spredixcan results
compiled_page_spredixcan <- compiled_page_spredixcan %>% mutate(study='PAGE')
compiled_panukbb_spredixcan <- compiled_panukbb_spredixcan %>% mutate(study='PanUKBB')

# reading significant non-zero heritability estimates
h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01)

# filtering spredixcan outputs so it only contains genes in the h2 df
for (tis in c('PBMC','Mono','Tcell')){
  genes.in.tissue <- h2estimates %>% filter(tissue==tis) %>% select(gene) %>% unique()
  page.filt <- compiled_page_spredixcan %>% filter(tissue==tis, gene %in% genes.in.tissue$gene)
  panukbb.filt <- compiled_panukbb_spredixcan %>% filter(tissue==tis, gene %in% genes.in.tissue$gene)
  
  if (exists('page_compiled_filtered')){
    page_compiled_filtered <- rbind(page_compiled_filtered, page.filt)
  } else {page_compiled_filtered <- page.filt}
  
  if (exists('panukbb_compiled_filtered')){
    panukbb_compiled_filtered <- rbind(panukbb_compiled_filtered, panukbb.filt)
  } else {panukbb_compiled_filtered <- panukbb.filt}
}

# compile files
compiled_both_studies <- rbind(page_compiled_filtered, panukbb_compiled_filtered)
compiled_both_studies$model <- gsub('mashr', 'MASHR', compiled_both_studies$model)
compiled_both_studies$model <- gsub('matrixeQTL', 'MatrixeQTL', compiled_both_studies$model)
compiled_both_studies$model <- gsub('elasticnet_unfiltered', 'EN', compiled_both_studies$model)

# save
fwrite(compiled_both_studies, 'SPrediXcan_PAGE_PanUKBB_raw_results.csv', quote=F, sep=',')