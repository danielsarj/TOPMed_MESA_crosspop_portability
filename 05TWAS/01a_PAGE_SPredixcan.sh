#!/bin/bash
#first arg = trait
#second arg = tissue
#third arg = pop

## PAGE

#Elastic Net
/home/daraujo1/anaconda3/bin/python /home/daraujo1/MetaXcan/software/SPrediXcan.py \
--model_db_path /home/cnguyen11/data/topmed_expression_whole_genome/04_elasticnet/outputs/dbs/$2_$3_PCAIR_elasticnet_models_unfiltered.db \
--covariance  /home/cnguyen11/data/topmed_expression_whole_genome/04_elasticnet/outputs/dbs/$2_$3_PCAIR_elasticnet_models_unfiltered_covariances.txt \
--gwas_file /home/daraujo1/data/Summary_statistics/WojcikG/harmonised_files/WojcikG_$1.harmonised.SNPsIDmod.txt.gz \
--snp_column SNP_hg38 --effect_allele_column Effect-allele --non_effect_allele_column Other-allele --se_column SE --beta_column Beta \
--keep_non_rsid --output_file /home/daraujo1/data/Summary_statistics/WojcikG/SPrediXcan/WojcikG_EN_$1_$2_$3.csv

#MASHR
/home/daraujo1/anaconda3/bin/python /home/daraujo1/MetaXcan/software/SPrediXcan.py \
--model_db_path /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/04_spredixcan/spredixcan_input/$2_$3_mashr_baseline.db \
--covariance  /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/03_covariance/covariance_files/$2_$3_mashr_baseline_WG_covariances.txt  \
--gwas_file /home/daraujo1/data/Summary_statistics/WojcikG/harmonised_files/WojcikG_$1.harmonised.SNPsIDmod.txt.gz \
--snp_column SNP_hg38 --effect_allele_column Effect-allele --non_effect_allele_column Other-allele --se_column SE --beta_column Beta \
--keep_non_rsid --output_file /home/daraujo1/data/Summary_statistics/WojcikG/SPrediXcan/WojcikG_MASHR_$1_$2_$3.csv

#MatrixeQTL
/home/daraujo1/anaconda3/bin/python /home/daraujo1/MetaXcan/software/SPrediXcan.py \
--model_db_path /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/04_spredixcan/spredixcan_input/$2_$3_matrixeQTL_baseline.db \
--covariance  /home/cnguyen11/data/topmed_expression_whole_genome/06_SPrediXcan/03_covariance/covariance_files/$2_$3_matrixeQTL_baseline_WG_covariances.txt  \
--gwas_file /home/daraujo1/data/Summary_statistics/WojcikG/harmonised_files/WojcikG_$1.harmonised.SNPsIDmod.txt.gz \
--snp_column SNP_hg38 --effect_allele_column Effect-allele --non_effect_allele_column Other-allele --se_column SE --beta_column Beta \
--keep_non_rsid --output_file /home/daraujo1/data/Summary_statistics/WojcikG/SPrediXcan/WojcikG_MatrixeQTL_$1_$2_$3.csv
