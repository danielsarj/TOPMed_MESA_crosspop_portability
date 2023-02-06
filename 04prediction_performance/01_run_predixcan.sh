#!/bin/bash
# ${1} = Geuvadis population
# ${2} = MESA tissue
# ${3} = MESA population
source ~/.bashrc
conda activate imlabtools

python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/cnguyen11/data/topmed_expression_whole_genome/04_elasticnet/outputs/dbs/${1}_${2}_PAGE_PCAIR_baseline_models_unfiltered.db \
--text_genotypes /home/daraujo1/data/METS/vcfs/dosages/METS.imputed.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/METS/vcfs/dosages/samples.txt \
--prediction_output /home/daraujo1/data/METS/WGS_PrediXcan/METS.imputed.${1}.${2}_ENunf_predict.txt \
--prediction_summary_output /home/daraujo1/data/METS/WGS_PrediXcan/METS.imputed.${1}.${2}_ENunf_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/TOPMed_MESA/mashr/mashr_db/WGS_files/mashr_models/${1}_${2}_mashr_baseline.db \
--text_genotypes /home/daraujo1/data/METS/vcfs/dosages/METS.imputed.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/METS/vcfs/dosages/samples.txt \
--prediction_output /home/daraujo1/data/METS/WGS_PrediXcan/METS.imputed.${1}.${2}_mashr_predict.txt \
--prediction_summary_output /home/daraujo1/data/METS/WGS_PrediXcan/METS.imputed.${1}.${2}_mashr_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/TOPMed_MESA/mashr/matrixeQTL/WGS_files/baseline_models/${1}_${2}_matrixeQTL_baseline.db \
--text_genotypes /home/daraujo1/data/METS/vcfs/dosages/METS.imputed.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/METS/vcfs/dosages/samples.txt \
--prediction_output /home/daraujo1/data/METS/WGS_PrediXcan/METS.imputed.${1}.${2}_matrixeqtl_predict.txt \
--prediction_summary_output /home/daraujo1/data/METS/WGS_PrediXcan/METS.imputed.${1}.${2}_matrixeqtl_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

