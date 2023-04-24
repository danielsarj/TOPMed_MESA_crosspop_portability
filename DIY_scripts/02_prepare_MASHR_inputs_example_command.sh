#!/bin/bash

Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/02_prepare_MASHR_inputs.R \
-l ~/TOPMed_MESA_crosspop_portability/sample_data/pop_files.txt \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_inputs