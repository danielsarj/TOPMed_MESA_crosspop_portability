#!/bin/bash

Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/04_make_MASHR_db.R \
-f ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_outputs \
-c GBR-YRI \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_models