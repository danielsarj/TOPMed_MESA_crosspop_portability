#!/bin/bash

Rscript ~/github_mashr_project/DIY_scripts/04_make_MASHR_db.R \
-f ~/github_mashr_project/sample_data/MASHR_outputs \
-c GBR-YRI \
-g ~/github_mashr_project/sample_data/gene_annotation.txt \
-o ~/github_mashr_project/sample_data/MASHR_models