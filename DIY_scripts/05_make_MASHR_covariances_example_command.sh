#!/bin/bash

for pop in GBR YRI
do
  for chr in 22
  do
    Rscript ~/github_mashr_project/DIY_scripts/05_make_MASHR_covariances.R \
    -d ~/github_mashr_project/sample_data/dosages/GEUVADIS_${pop}_chr${chr}_dosage_unfiltered.txt.gz \
    -g ~/github_mashr_project/sample_data/gene_annotation.txt \
    -t ${pop}_MASHR \
    -c ${chr} \
    -m ~/github_mashr_project/sample_data/MASHR_models/${pop}_MASHR_weights.txt.gz \
    -o ~/github_mashr_project/sample_data/MASHR_models \
    -w 1000000
  done
done