library(data.table)
library(tidyverse)
library(viridis)
'%&%' = function(a,b) paste(a,b,sep='')
setwd('/home/daniel/Geuvadis/WGS_predixcan')

for (tissue in c('PBMC','Mono','Tcell')){
  if (tissue=='PBMC'){
    #loading geuvadis predictions
    GeuALL_AFA_mashr<-fread('Geuvadis.ALL_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_EUR_mashr<-fread('Geuvadis.ALL_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_CHN_mashr<-fread('Geuvadis.ALL_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_HIS_mashr<-fread('Geuvadis.ALL_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_AFA_matrixeqtl<-fread('Geuvadis.ALL_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_EUR_matrixeqtl<-fread('Geuvadis.ALL_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_CHN_matrixeqtl<-fread('Geuvadis.ALL_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_HIS_matrixeqtl<-fread('Geuvadis.ALL_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_AFA_mashr<-fread('Geuvadis.CEU_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_EUR_mashr<-fread('Geuvadis.CEU_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_CHN_mashr<-fread('Geuvadis.CEU_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_HIS_mashr<-fread('Geuvadis.CEU_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_AFA_matrixeqtl<-fread('Geuvadis.CEU_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_EUR_matrixeqtl<-fread('Geuvadis.CEU_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_CHN_matrixeqtl<-fread('Geuvadis.CEU_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_HIS_matrixeqtl<-fread('Geuvadis.CEU_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_AFA_mashr<-fread('Geuvadis.FIN_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_EUR_mashr<-fread('Geuvadis.FIN_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_CHN_mashr<-fread('Geuvadis.FIN_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_HIS_mashr<-fread('Geuvadis.FIN_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_AFA_matrixeqtl<-fread('Geuvadis.FIN_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_EUR_matrixeqtl<-fread('Geuvadis.FIN_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_CHN_matrixeqtl<-fread('Geuvadis.FIN_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_HIS_matrixeqtl<-fread('Geuvadis.FIN_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_mashr<-fread('Geuvadis.GBR_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_mashr<-fread('Geuvadis.GBR_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_CHN_mashr<-fread('Geuvadis.GBR_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_mashr<-fread('Geuvadis.GBR_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_matrixeqtl<-fread('Geuvadis.GBR_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_matrixeqtl<-fread('Geuvadis.GBR_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_CHN_matrixeqtl<-fread('Geuvadis.GBR_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_matrixeqtl<-fread('Geuvadis.GBR_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_AFA_mashr<-fread('Geuvadis.TSI_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_EUR_mashr<-fread('Geuvadis.TSI_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_CHN_mashr<-fread('Geuvadis.TSI_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_HIS_mashr<-fread('Geuvadis.TSI_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_AFA_matrixeqtl<-fread('Geuvadis.TSI_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_EUR_matrixeqtl<-fread('Geuvadis.TSI_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_CHN_matrixeqtl<-fread('Geuvadis.TSI_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_HIS_matrixeqtl<-fread('Geuvadis.TSI_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_mashr<-fread('Geuvadis.YRI_PBMC.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_mashr<-fread('Geuvadis.YRI_PBMC.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_CHN_mashr<-fread('Geuvadis.YRI_PBMC.CHN_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_mashr<-fread('Geuvadis.YRI_PBMC.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_matrixeqtl<-fread('Geuvadis.YRI_PBMC.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_matrixeqtl<-fread('Geuvadis.YRI_PBMC.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_CHN_matrixeqtl<-fread('Geuvadis.YRI_PBMC.CHN_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_matrixeqtl<-fread('Geuvadis.YRI_PBMC.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_AFA_ENunf<-fread('Geuvadis.ALL_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuALL_EUR_ENunf<-fread('Geuvadis.ALL_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuALL_CHN_ENunf<-fread('Geuvadis.ALL_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuALL_HIS_ENunf<-fread('Geuvadis.ALL_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_AFA_ENunf<-fread('Geuvadis.CEU_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_EUR_ENunf<-fread('Geuvadis.CEU_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_CHN_ENunf<-fread('Geuvadis.CEU_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_HIS_ENunf<-fread('Geuvadis.CEU_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_AFA_ENunf<-fread('Geuvadis.FIN_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_EUR_ENunf<-fread('Geuvadis.FIN_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_CHN_ENunf<-fread('Geuvadis.FIN_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_HIS_ENunf<-fread('Geuvadis.FIN_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_ENunf<-fread('Geuvadis.GBR_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_ENunf<-fread('Geuvadis.GBR_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_CHN_ENunf<-fread('Geuvadis.GBR_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_ENunf<-fread('Geuvadis.GBR_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_AFA_ENunf<-fread('Geuvadis.TSI_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_EUR_ENunf<-fread('Geuvadis.TSI_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_CHN_ENunf<-fread('Geuvadis.TSI_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_HIS_ENunf<-fread('Geuvadis.TSI_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_ENunf<-fread('Geuvadis.YRI_PBMC.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC AFA', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_ENunf<-fread('Geuvadis.YRI_PBMC.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC EUR', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_CHN_ENunf<-fread('Geuvadis.YRI_PBMC.CHN_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC CHN', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_ENunf<-fread('Geuvadis.YRI_PBMC.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop='PBMC HIS', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    
    #binding rows by method
    mashr_spearmans<-bind_rows(GeuALL_AFA_mashr, GeuALL_EUR_mashr, GeuALL_CHN_mashr, GeuALL_HIS_mashr,
                               GeuCEU_AFA_mashr, GeuCEU_EUR_mashr, GeuCEU_CHN_mashr, GeuCEU_HIS_mashr,
                               GeuFIN_AFA_mashr, GeuFIN_EUR_mashr, GeuFIN_CHN_mashr, GeuFIN_HIS_mashr,
                               GeuGBR_AFA_mashr, GeuGBR_EUR_mashr, GeuGBR_CHN_mashr, GeuGBR_HIS_mashr,
                               GeuTSI_AFA_mashr, GeuTSI_EUR_mashr, GeuTSI_CHN_mashr, GeuTSI_HIS_mashr,
                               GeuYRI_AFA_mashr, GeuYRI_EUR_mashr, GeuYRI_CHN_mashr, GeuYRI_HIS_mashr)
    matrixeqtl_spearmans<-bind_rows(GeuALL_AFA_matrixeqtl, GeuALL_EUR_matrixeqtl, GeuALL_CHN_matrixeqtl, GeuALL_HIS_matrixeqtl,
                                    GeuCEU_AFA_matrixeqtl, GeuCEU_EUR_matrixeqtl, GeuCEU_CHN_matrixeqtl, GeuCEU_HIS_matrixeqtl,
                                    GeuFIN_AFA_matrixeqtl, GeuFIN_EUR_matrixeqtl, GeuFIN_CHN_matrixeqtl, GeuFIN_HIS_matrixeqtl,
                                    GeuGBR_AFA_matrixeqtl, GeuGBR_EUR_matrixeqtl, GeuGBR_CHN_matrixeqtl, GeuGBR_HIS_matrixeqtl,
                                    GeuTSI_AFA_matrixeqtl, GeuTSI_EUR_matrixeqtl, GeuTSI_CHN_matrixeqtl, GeuTSI_HIS_matrixeqtl,
                                    GeuYRI_AFA_matrixeqtl, GeuYRI_EUR_matrixeqtl, GeuYRI_CHN_matrixeqtl, GeuYRI_HIS_matrixeqtl)
    ENunf_spearmans<-bind_rows(GeuALL_AFA_ENunf, GeuALL_EUR_ENunf, GeuALL_CHN_ENunf, GeuALL_HIS_ENunf,
                               GeuCEU_AFA_ENunf, GeuCEU_EUR_ENunf, GeuCEU_CHN_ENunf, GeuCEU_HIS_ENunf,
                               GeuFIN_AFA_ENunf, GeuFIN_EUR_ENunf, GeuFIN_CHN_ENunf, GeuFIN_HIS_ENunf,
                               GeuGBR_AFA_ENunf, GeuGBR_EUR_ENunf, GeuGBR_CHN_ENunf, GeuGBR_HIS_ENunf,
                               GeuTSI_AFA_ENunf, GeuTSI_EUR_ENunf, GeuTSI_CHN_ENunf, GeuTSI_HIS_ENunf,
                               GeuYRI_AFA_ENunf, GeuYRI_EUR_ENunf, GeuYRI_CHN_ENunf, GeuYRI_HIS_ENunf)
    #joining everything in a single df
    joint_spearmans <- rbind(ENunf_spearmans, mashr_spearmans, matrixeqtl_spearmans)
    
    #also join by gene
    mashr_spearmans <- mashr_spearmans %>% rename(mashr_estimate = estimate)
    matrixeqtl_spearmans <- matrixeqtl_spearmans %>% rename(matrixeqtl_estimate = estimate)
    ENunf_spearmans <- ENunf_spearmans %>% rename(ENunf_estimate = estimate)
    joint_spearmans_by_gene <- full_join(matrixeqtl_spearmans, mashr_spearmans, by=c("gene_id","model_pop","Pop")) %>% 
      full_join(ENunf_spearmans, by=c("gene_id","model_pop","Pop")) %>%
      select("gene_id", "model_pop", "Pop", "matrixeqtl_estimate", "mashr_estimate", "ENunf_estimate")
    
    #heatmap w/o intersection of genes
    non_intersect_median_df <- joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(non_intersect_median_df, 'PBMC_median_spearmans_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = non_intersect_median_df, aes(x = model, y = Pop, fill = median_estimate)) + facet_wrap(~model_pop) +
      geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis (all genes)')
    ggsave('Geuvadis_PBMC_median_heatmap.pdf', height=5, width=6)
    ggsave('Geuvadis_PBMC_median_heatmap.png', height=5, width=6)
    
    #getting the intersection across models per MESA pop
    for (pop in c('AFA','EUR','HIS','CHN')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      matrixeqtl_genes <- matrixeqtl_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% intersect(matrixeqtl_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else if (pop=='HIS'){
        his_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else {
        chn_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      }
    }
    
    #heatmap w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions, chn_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, 'PBMC_median_spearmans_intersection_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_median_df, aes(x = model, y = Pop, fill = median_estimate)) + facet_wrap(~model_pop) +
      geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis (intersection of genes)')
    ggsave('Geuvadis_PBMC_median_heatmap_intersection.pdf', height=5, width=6)
    ggsave('Geuvadis_PBMC_median_heatmap_intersection.png', height=5, width=6)
    
    #wilcox rank sum test after intersection
    for (m in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (p in c('ALL', 'FIN', 'CEU', 'GBR', 'TSI', 'YRI')){
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
    fwrite(wilcox.tests.pvalues, 'wilcox_PBMC_pvalues_intersection.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)

    #getting the intersection across models per MESA pop - this time, excluding MatrixeQTL
    for (pop in c('AFA','EUR','HIS','CHN')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else if (pop=='HIS'){
        his_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else {
        chn_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      }
    }
    
    #heatmap w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions, chn_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, 'PBMC_median_spearmans_intersection_noMatrix_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_median_df, aes(x = model, y = Pop, fill = median_estimate)) + facet_wrap(~model_pop) +
      geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis (intersection of genes)')
    ggsave('Geuvadis_PBMC_median_heatmap_intersection_noMatrix.pdf', height=5, width=6)
    ggsave('Geuvadis_PBMC_median_heatmap_intersection_noMatrix.png', height=5, width=6)
    
    #wilcox rank sum test after intersection w/o MatrixeQTL
    for (m in c('AFA', 'EUR', 'HIS', 'CHN')){
      for (p in c('ALL', 'FIN', 'CEU', 'GBR', 'TSI', 'YRI')){
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
    fwrite(wilcox.tests.pvalues, 'wilcox_PBMC_pvalues_intersection_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
  } else {
    GeuALL_AFA_mashr<-fread('Geuvadis.ALL_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_EUR_mashr<-fread('Geuvadis.ALL_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_HIS_mashr<-fread('Geuvadis.ALL_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='ALL', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuALL_AFA_matrixeqtl<-fread('Geuvadis.ALL_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_EUR_matrixeqtl<-fread('Geuvadis.ALL_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_HIS_matrixeqtl<-fread('Geuvadis.ALL_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='ALL', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_AFA_mashr<-fread('Geuvadis.CEU_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_EUR_mashr<-fread('Geuvadis.CEU_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_HIS_mashr<-fread('Geuvadis.CEU_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='CEU', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_AFA_matrixeqtl<-fread('Geuvadis.CEU_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_EUR_matrixeqtl<-fread('Geuvadis.CEU_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_HIS_matrixeqtl<-fread('Geuvadis.CEU_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='CEU', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_AFA_mashr<-fread('Geuvadis.FIN_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_EUR_mashr<-fread('Geuvadis.FIN_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_HIS_mashr<-fread('Geuvadis.FIN_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='FIN', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_AFA_matrixeqtl<-fread('Geuvadis.FIN_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_EUR_matrixeqtl<-fread('Geuvadis.FIN_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_HIS_matrixeqtl<-fread('Geuvadis.FIN_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='FIN', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_mashr<-fread('Geuvadis.GBR_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_mashr<-fread('Geuvadis.GBR_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_mashr<-fread('Geuvadis.GBR_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='GBR', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_matrixeqtl<-fread('Geuvadis.GBR_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_matrixeqtl<-fread('Geuvadis.GBR_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_matrixeqtl<-fread('Geuvadis.GBR_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='GBR', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_AFA_mashr<-fread('Geuvadis.TSI_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_EUR_mashr<-fread('Geuvadis.TSI_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_HIS_mashr<-fread('Geuvadis.TSI_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='TSI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_AFA_matrixeqtl<-fread('Geuvadis.TSI_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_EUR_matrixeqtl<-fread('Geuvadis.TSI_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_HIS_matrixeqtl<-fread('Geuvadis.TSI_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='TSI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_mashr<-fread('Geuvadis.YRI_'%&% tissue %&%'.AFA_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_mashr<-fread('Geuvadis.YRI_'%&% tissue %&%'.EUR_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_mashr<-fread('Geuvadis.YRI_'%&% tissue %&%'.HIS_mashr_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='YRI', model='Mashr') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_matrixeqtl<-fread('Geuvadis.YRI_'%&% tissue %&%'.AFA_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_matrixeqtl<-fread('Geuvadis.YRI_'%&% tissue %&%'.EUR_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_matrixeqtl<-fread('Geuvadis.YRI_'%&% tissue %&%'.HIS_matrixeqtl_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='YRI', model='MatrixeQTL') %>% select(-c(test_R2_avg, p.value))
    GeuALL_AFA_ENunf<-fread('Geuvadis.ALL_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuALL_EUR_ENunf<-fread('Geuvadis.ALL_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuALL_HIS_ENunf<-fread('Geuvadis.ALL_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='ALL', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_AFA_ENunf<-fread('Geuvadis.CEU_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_EUR_ENunf<-fread('Geuvadis.CEU_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuCEU_HIS_ENunf<-fread('Geuvadis.CEU_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='CEU', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_AFA_ENunf<-fread('Geuvadis.FIN_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_EUR_ENunf<-fread('Geuvadis.FIN_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuFIN_HIS_ENunf<-fread('Geuvadis.FIN_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='FIN', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_AFA_ENunf<-fread('Geuvadis.GBR_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_EUR_ENunf<-fread('Geuvadis.GBR_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuGBR_HIS_ENunf<-fread('Geuvadis.GBR_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='GBR', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_AFA_ENunf<-fread('Geuvadis.TSI_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_EUR_ENunf<-fread('Geuvadis.TSI_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuTSI_HIS_ENunf<-fread('Geuvadis.TSI_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='TSI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_AFA_ENunf<-fread('Geuvadis.YRI_'%&% tissue %&%'.AFA_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' AFA', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_EUR_ENunf<-fread('Geuvadis.YRI_'%&% tissue %&%'.EUR_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' EUR', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    GeuYRI_HIS_ENunf<-fread('Geuvadis.YRI_'%&% tissue %&%'.HIS_ENunf_expression_spearman_correlation.txt',header=T) %>% mutate(model_pop=tissue %&%' HIS', Pop='YRI', model='ENunf') %>% select(-c(test_R2_avg, p.value))
    
    #binding rows by method
    mashr_spearmans<-bind_rows(GeuALL_AFA_mashr, GeuALL_EUR_mashr, GeuALL_HIS_mashr,
                               GeuCEU_AFA_mashr, GeuCEU_EUR_mashr, GeuCEU_HIS_mashr,
                               GeuFIN_AFA_mashr, GeuFIN_EUR_mashr, GeuFIN_HIS_mashr,
                               GeuGBR_AFA_mashr, GeuGBR_EUR_mashr, GeuGBR_HIS_mashr,
                               GeuTSI_AFA_mashr, GeuTSI_EUR_mashr, GeuTSI_HIS_mashr,
                               GeuYRI_AFA_mashr, GeuYRI_EUR_mashr, GeuYRI_HIS_mashr)
    matrixeqtl_spearmans<-bind_rows(GeuALL_AFA_matrixeqtl, GeuALL_EUR_matrixeqtl, GeuALL_HIS_matrixeqtl,
                                    GeuCEU_AFA_matrixeqtl, GeuCEU_EUR_matrixeqtl, GeuCEU_HIS_matrixeqtl,
                                    GeuFIN_AFA_matrixeqtl, GeuFIN_EUR_matrixeqtl, GeuFIN_HIS_matrixeqtl,
                                    GeuGBR_AFA_matrixeqtl, GeuGBR_EUR_matrixeqtl, GeuGBR_HIS_matrixeqtl,
                                    GeuTSI_AFA_matrixeqtl, GeuTSI_EUR_matrixeqtl, GeuTSI_HIS_matrixeqtl,
                                    GeuYRI_AFA_matrixeqtl, GeuYRI_EUR_matrixeqtl, GeuYRI_HIS_matrixeqtl)
    ENunf_spearmans<-bind_rows(GeuALL_AFA_ENunf, GeuALL_EUR_ENunf, GeuALL_HIS_ENunf,
                               GeuCEU_AFA_ENunf, GeuCEU_EUR_ENunf, GeuCEU_HIS_ENunf,
                               GeuFIN_AFA_ENunf, GeuFIN_EUR_ENunf, GeuFIN_HIS_ENunf,
                               GeuGBR_AFA_ENunf, GeuGBR_EUR_ENunf, GeuGBR_HIS_ENunf,
                               GeuTSI_AFA_ENunf, GeuTSI_EUR_ENunf, GeuTSI_HIS_ENunf,
                               GeuYRI_AFA_ENunf, GeuYRI_EUR_ENunf, GeuYRI_HIS_ENunf)
    #joining everything
    joint_spearmans <- rbind(ENunf_spearmans, mashr_spearmans, matrixeqtl_spearmans)
    
    #also join by gene
    mashr_spearmans <- mashr_spearmans %>% rename(mashr_estimate = estimate)
    matrixeqtl_spearmans <- matrixeqtl_spearmans %>% rename(matrixeqtl_estimate = estimate)
    ENunf_spearmans <- ENunf_spearmans %>% rename(ENunf_estimate = estimate)
    joint_spearmans_by_gene <- full_join(matrixeqtl_spearmans, mashr_spearmans, by=c("gene_id","model_pop","Pop")) %>% 
      full_join(ENunf_spearmans, by=c("gene_id","model_pop","Pop")) %>%
      select("gene_id", "model_pop", "Pop", "matrixeqtl_estimate", "mashr_estimate", "ENunf_estimate")
    
    #heatmap w/o intersection of genes
    non_intersect_median_df <- joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(non_intersect_median_df, tissue %&%'_median_spearmans_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = non_intersect_median_df, aes(x = model, y = Pop, fill = median_estimate)) + facet_wrap(~model_pop) +
      geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis (all genes)')
    ggsave('Geuvadis_'%&% tissue %&%'_median_heatmap.pdf', height=3, width=8)
    ggsave('Geuvadis_'%&% tissue %&%'_median_heatmap.png', height=3, width=8)
    
    #getting the intersection across models per MESA pop
    for (pop in c('AFA','EUR','HIS')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      matrixeqtl_genes <- matrixeqtl_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% intersect(matrixeqtl_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      } else {
        his_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes) %>% drop_na()
      }
    }
    
    #heatmap w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, tissue %&%'_median_spearmans_intersection_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_median_df, aes(x = model, y = Pop, fill = median_estimate)) + facet_wrap(~model_pop) +
      geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis (intersection of genes)')
    ggsave('Geuvadis_'%&% tissue %&%'_median_heatmap_intersection.pdf', height=3, width=8)
    ggsave('Geuvadis_'%&% tissue %&%'_median_heatmap_intersection.png', height=3, width=8)
    
    #wilcox rank sum test after intersection
    for (m in c('AFA', 'EUR', 'HIS')){
      for (p in c('ALL', 'FIN', 'CEU', 'GBR', 'TSI', 'YRI')){
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
    fwrite(wilcox.tests.pvalues, 'wilcox_'%&% tissue %&%'_pvalues_intersection.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
    
    #getting the intersection across models per MESA pop - this time, excluding MatrixeQTL
    for (pop in c('AFA','EUR','HIS')){
      mashr_genes <- mashr_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      ENunf_genes <- ENunf_spearmans %>% filter(grepl(pop, model_pop)) %>% drop_na() %>% select(gene_id) %>% unique()
      intersection_genes <- intersect(mashr_genes, ENunf_genes) %>% pull() %>% unique()
      
      if (pop=='AFA'){
        afa_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else if (pop=='EUR'){
        eur_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      } else {
        his_predictions <- joint_spearmans %>% filter(grepl(pop, model_pop), gene_id %in% intersection_genes, !model=='MatrixeQTL') %>% drop_na()
      }
    }
    
    #heatmap w/ intersection of genes
    intersect_joint_spearmans <- rbind(afa_predictions, eur_predictions, his_predictions)
    intersect_median_df <- intersect_joint_spearmans %>% select(model_pop, Pop, estimate, model) %>%
      group_by(model_pop, Pop, model) %>% drop_na() %>% summarize(n = n(), median_estimate = median(estimate, na.rm = T))
    fwrite(intersect_median_df, tissue %&%'_median_spearmans_intersection_noMatrix_df.txt', col.names = T, sep = ' ')
    
    ggplot(data = intersect_median_df, aes(x = model, y = Pop, fill = median_estimate)) + facet_wrap(~model_pop) +
      geom_tile() + scale_fill_viridis() + labs(title='Median Spearman correlation in Geuvadis (intersection of genes)')
    ggsave('Geuvadis_'%&% tissue %&%'_median_heatmap_intersection_noMatrix.pdf', height=3, width=6)
    ggsave('Geuvadis_'%&% tissue %&%'_median_heatmap_intersection_noMatrix.png', height=3, width=6)
    
    #wilcox rank sum test after intersection w/o MatrixeQTL
    for (m in c('AFA', 'EUR', 'HIS')){
      for (p in c('ALL', 'FIN', 'CEU', 'GBR', 'TSI', 'YRI')){
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
    fwrite(wilcox.tests.pvalues, 'wilcox_'%&% tissue %&%'_pvalues_intersection_noMatrixeQTL.txt', col.names = T, sep = ' ')
    rm(wilcox.tests.pvalues)
  }
}
