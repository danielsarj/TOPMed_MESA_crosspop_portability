#!/bin/bash

Rscript DIY_scripts/04_make_MASHR_db.R \
-f sample_data/MASHR_outputs \
-c GBR-YRI \
-g sample_data/gene_annotation.txt \
-o sample_data/MASHR_models
