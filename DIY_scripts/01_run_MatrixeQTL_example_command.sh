#!/bin/bash

for pop in GBR YRI
do
  for chr in 22
  do
    Rscript /DIY_scripts/01_run_MatrixeQTL.R \
    -d sample_data/dosages/GEUVADIS_${pop}_chr${chr}_dosage_filtered.txt.gz \
    -e sample_data/exp_tbls/GEUVADIS_${pop}_expression.txt.gz \
    -g sample_data/gene_annotation.txt \
    -t ${pop}_chr${chr}_cis \
    -o sample_data/cis_ES/ \
    -w 1000000
  done
done
