#!/bin/bash
#first arg = trait
#second arg = tissue
#third arg = pop
#fourth arg = meta/AFR from Pan_UK_Biobank

# PanUKBB

#Elastic Net
/home/daraujo1/anaconda3/bin/python /home/daraujo1/MetaXcan/software/SPrediXcan.py \
--model_db_path /home/daraujo1/data/MESA_local_ancestry/models/$2_$3_EN_baseline.db \
--covariance /home/daraujo1/data/MESA_local_ancestry/models/$2_$3_EN_baseline_covariances.txt \
--gwas_file /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/04_spredixcan/Pan_UK_Biobank/raw_summary_statistics/Pan.UK.Biobank.$1.processed.lifted_b38.txt \
--snp_column snp_id_38 --effect_allele_column alt --non_effect_allele_column ref --se_column se_$4 --beta_column beta_$4 \
--keep_non_rsid --output_file /home/daraujo1/data/Summary_statistics/PanUKBB/LA_SPrediXcan/Pan.UK.Biobank_$4_EN_$1_$2_$3.csv

#MASHR
/home/daraujo1/anaconda3/bin/python /home/daraujo1/MetaXcan/software/SPrediXcan.py \
--model_db_path /home/daraujo1/data/MESA_local_ancestry/models/$2_$3_mashr_baseline.db \
--covariance /home/daraujo1/data/MESA_local_ancestry/models/$2_$3_mashr_baseline_covariances.txt \
--gwas_file /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/04_spredixcan/Pan_UK_Biobank/raw_summary_statistics/Pan.UK.Biobank.$1.processed.lifted_b38.txt \
--snp_column snp_id_38 --effect_allele_column alt --non_effect_allele_column ref --se_column se_$4 --beta_column beta_$4 \
--keep_non_rsid --output_file /home/daraujo1/data/Summary_statistics/PanUKBB/LA_SPrediXcan/Pan.UK.Biobank_$4_mashr_$1_$2_$3.csv

#MatrixeQTL
/home/daraujo1/anaconda3/bin/python /home/daraujo1/MetaXcan/software/SPrediXcan.py \
--model_db_path /home/daraujo1/data/MESA_local_ancestry/models/$2_$3_matrixeQTL_baseline.db \
--covariance /home/daraujo1/data/MESA_local_ancestry/models/$2_$3_matrixeQTL_baseline_covariances.txt \
--gwas_file /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/04_spredixcan/Pan_UK_Biobank/raw_summary_statistics/Pan.UK.Biobank.$1.processed.lifted_b38.txt \
--snp_column snp_id_38 --effect_allele_column alt --non_effect_allele_column ref --se_column se_$4 --beta_column beta_$4 \
--keep_non_rsid --output_file /home/daraujo1/data/Summary_statistics/PanUKBB/LA_SPrediXcan/Pan.UK.Biobank_$4_matrixeQTL_$1_$2_$3.csv
