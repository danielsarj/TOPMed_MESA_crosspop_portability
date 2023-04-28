#!/bin/bash

Rscript DIY_scripts/03_run_MASHR.R \
-i sample_data/MASHR_inputs \
-g sample_data/gene_annotation.txt \
-o sample_data/MASHR_outputs
