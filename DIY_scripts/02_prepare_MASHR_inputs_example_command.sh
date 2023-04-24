#!/bin/bash

Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/02_prepare_MASHR_inputs.R \
-l ~/TOPMed_MESA_crosspop_portability/sample_data/chr22_files.txt \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-c 22 \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_inputs