#!/bin/bash
# ${1} = Geuvadis population
# ${2} = MESA tissue
# ${3} = MESA population
source ~/.bashrc
conda activate imlabtools

#MatrixeQTL
python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/TOPMed_MESA/mashr/matrixeQTL/WGS_files/baseline_models/${2}_${3}_matrixeQTL_baseline.db \
--text_genotypes /home/daraujo1/data/Geuvadis/dosages/Geuvadis.${1}.plink.updated.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/Geuvadis/${1}_ids_formatted.txt \
--prediction_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_matrixeQTL_noth2filt_predict.txt \
--prediction_summary_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_matrixeQTL_noth2filt_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

#MASHR
python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/TOPMed_MESA/mashr/mashr_db/WGS_files/mashr_models/${2}_${3}_mashr_baseline.db \
--text_genotypes /home/daraujo1/data/Geuvadis/dosages/Geuvadis.${1}.plink.updated.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/Geuvadis/${1}_ids_formatted.txt \
--prediction_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_mashr_noth2filt_predict.txt \
--prediction_summary_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_mashr_noth2filt_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

#Elastic Net
python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/Geuvadis/WGS_predixcan/${2}_${3}_PAGE_PCAIR_baseline_models_unfiltered.db \
--text_genotypes /home/daraujo1/data/Geuvadis/dosages/Geuvadis.${1}.plink.updated.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/Geuvadis/${1}_ids_formatted.txt \
--prediction_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_EN_noth2filt_predict.txt \
--prediction_summary_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_EN_noth2filt_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

#TIGAR
python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/TOPMed_MESA/TIGAR/dbs/${2}_${3}_TIGAR_reduced_1e-04.db \
--text_genotypes /home/daraujo1/data/Geuvadis/dosages/Geuvadis.${1}.plink.updated.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/Geuvadis/${1}_ids_formatted.txt \
--prediction_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_TIGAR_noth2filt_predict.txt \
--prediction_summary_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_TIGAR_noth2filt_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

#JTI
python3 /home/daraujo1/MetaXcan/software/Predict.py \
--model_db_path /home/daraujo1/data/TOPMed_MESA/MR-JTI/dbs/${2}_${3}_JTI_baseline.db \
--text_genotypes /home/daraujo1/data/Geuvadis/dosages/Geuvadis.${1}.plink.updated.chr*.dosage.txt.gz \
--text_sample_ids /home/daraujo1/data/Geuvadis/${1}_ids_formatted.txt \
--prediction_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_JTI_predict.txt \
--prediction_summary_output /home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis.${1}_${2}.${3}_JTI_summary_predict.txt \
--verbosity 9 --throw --model_db_snp_key varID

