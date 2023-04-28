#!/bin/bash

for pop in GBR YRI
do
  for chr in 22
  do
    Rscript DIY_scripts/05_make_MASHR_covariances.R \
    -d sample_data/dosages/GEUVADIS_${pop}_chr${chr}_dosage_unfiltered.txt.gz \
    -g sample_data/gene_annotation.txt \
    -t ${pop}_MASHR \
    -c ${chr} \
    -m sample_data/MASHR_models/${pop}_MASHR_weights.txt.gz \
    -o sample_data/MASHR_models \
    -w 1000000
  done
done
