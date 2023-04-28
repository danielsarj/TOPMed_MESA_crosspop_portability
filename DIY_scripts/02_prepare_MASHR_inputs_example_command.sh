#!/bin/bash

Rscript DIY_scripts/02_prepare_MASHR_inputs.R \
-l sample_data/chr22_files.txt \
-g sample_data/gene_annotation.txt \
-c 22 \
-o sample_data/MASHR_inputs
