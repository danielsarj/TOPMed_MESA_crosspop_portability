#!/bin/bash

for pop in GBR YRI
do
  Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/01_run_MatrixeQTL.R \
  -d ~/TOPMed_MESA_crosspop_portability/sample_data/dosages/GEUVADIS_${pop}_chr22_dosage_filtered.txt.gz \
  -e ~/TOPMed_MESA_crosspop_portability/sample_data/exp_tbls/GEUVADIS_${pop}_expression.txt.gz \
  -g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
  -t ${pop}_chr22_cis \
  -o ~/TOPMed_MESA_crosspop_portability/sample_data/cis_ES/ \
  -w 1000000
done