#!/bin/bash

for pop in GBR YRI
do
  for chr in 22
  do
    Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/05_make_MASHR_covariances.R \
    -d ~/TOPMed_MESA_crosspop_portability/sample_data/dosages/GEUVADIS_${pop}_chr${chr}_dosage_unfiltered.txt.gz \
    -g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
    -t ${pop}_MASHR \
    -c ${chr} \
    -m ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_models/${pop}_MASHR_weights.txt.gz \
    -o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_models \
    -w 1000000
  done
done