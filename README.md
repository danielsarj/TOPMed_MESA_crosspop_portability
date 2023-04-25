# Optimizing models for use in multi-ethnic TWAS using MASHR

Most of the transcriptome prediction models are built using European-descent individuals' data. As a consequence, they are less accurate when applied to different populations. To improve cross-population transcriptome prediction accuracy, we built a pipeline that leverages cis-eQTLs effect size estimates across distinct populations using [multivariate adaptive shrinkage](https://doi.org/10.1038/s41588-018-0268-8). 

To build our transcriptome prediction models, we used data from the [TOPMed Multi-Ethnic Study of Atherosclerosis](https://doi.org/10.1093/aje/kwf113). Our models are designed to work with the [MetaXcan](https://github.com/hakyimlab/MetaXcan) family of TWAS tools and are freely available for anyone to use. They can be found in our [Zenodo](https://doi.org/10.5281/zenodo.7551845) repository. 

For more information about model design and performance, please refer to the research article "_Multivariate adaptive shrinkage improves cross-population transcriptome prediction for transcriptome-wide association studies in underrepresented populations_" by Daniel S. Araujo _et al_. [bioRxiv](https://doi.org/10.1101/2023.02.09.527747)

## Table of Contents

* [Do-it-yourself](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#do-it-yourself)
  * [Pipeline requirements](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#pipeline-requirements)
  * [Step 1. Estimate unadjusted cis-eQTL effect sizes](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#step-1-estimate-unadjusted-cis-eqtl-effect-sizes)
  * [Step 2. Prepare MASHR inputs](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#step-2-prepare-mashr-inputs)
  * [Step 3. Run MASHR](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#step-3-run-mashr)
  * [Step 4. Make model files for MetaXcan](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#step-4-make-model-files-for-metaxcan)
  * [Step 5. Make model covariance files for MetaXcan](https://github.com/danielsarj/TOPMed_MESA_crosspop_portability#step-5-make-model-covariance-files-for-metaxcan)

## Do-it-yourself 

Although we provide our own set of MASHR-based transcriptome models built using TOPMed MESA for use in TWAS, we encourage users to build their own set of transcriptome models in case they want to. For this, we supply a set of scripts that can be found in `./DIY_scripts/`. We also provide a small test data, found in `./sample_data/`, so users can easily learn the format of input files for the pipeline we designed. The test data is a small subset of the [Geuvadis dataset](https://doi.org/10.1038%2Fnature12531).

Note, we chose to provide our pipeline in a set of scripts (instead of a single one) to simplify debugging by users in case one step is not working as intended. Additionally, this also facilitates running some steps in parallel.

### Pipeline requirements

The pipeline is built on R and requires the following R packages:

* [argparse](https://cran.r-project.org/web/packages/argparse/index.html) 
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
* [mashr](https://cran.r-project.org/web/packages/mashr/index.html)
* [MatrixEQTL](https://cran.r-project.org/web/packages/MatrixEQTL/index.html)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
* [RSQLite](https://cran.r-project.org/web/packages/RSQLite/index.html)
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
* [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
___

### Step 1. Estimate unadjusted cis-eQTL effect sizes

The first step is to estimate cis-eQTL effect sizes individually for each condition. For this, use `./DIY_scripts/01_run_MatrixeQTL.R`. We recommend running this step in parallel for each population and chromosome. The arguments to run the first step of the pipeline are:

* `-d`, `--snpdosage`: Path to the SNP dosage file. 
* `-e`, `--geneexpression`: Path to the gene expression file. 
* `-g`, `--geneannotation`: Path to the gene annotation file. 
* `-t`, `--tag`: Prefix for the output file. 
* `-o`, `--outputdir`: Output directory path. By default, it is `./`. 
* `-w`, `--window`: Size of the cis-SNP window from the TSS and TES. By default, it is `1e6`. 

### Input files:
1. SNP dosage file: it is a space- or tab-separated file, containing 5+N columns, in which N is the number of individuals in the dataset. Each row is a different SNP, and from the 6th column forward, the columns contain each sample's dosage concerning the SNP. We highly recommend removing ambiguous and multi-allelic SNPs. It is okay if a dosage file for one population contains SNPs missing in other populations. Each population should have its own dosage file.

chr | snp_ID | pos | ref_allele | alt_allele | Sample1
--- | --- | --- | --- | --- | ---
22 | chr22:26329350:T:C | 26329350 | T | C | 0
22 | chr22:26329903:C:T | 26329903 | C | T | 0

2. Gene expression file: it is a space- or tab-separated file, containing 1+N columns, in which N is the number of individuals in the dataset. Each row is a different gene, and from the 2nd column forward, the columns contain each sample's expression levels. We recommend providing the residual of normalized expression after adjusting for covariates. Additionally, samples should be in the same order as the ones in the dosage file. This file should only contain genes that are found in every population gene expression file. Each population should have its own gene expression file.

gene_id | Sample1
--- | ---
ENSG00000128203 | -0.929926753044128
ENSG00000100099 | -0.747495412826538

3. Gene annotation file: it is a space- or tab-separated file containing 6 columns. Each row is a different gene and contains information about their ID, name, TSS, TES, and type. The ID under the `gene_id` column should match the ID in the column of the same name from the gene expression file. This file should only contain genes that are found in every population gene expression file.  

chr | gene_id | gene_name | start | end | gene_type
--- | --- | --- | --- | --- | ---
22 | ENSG00000128203 | ASPHD2 | 26429260 | 26445015 | protein_coding
22 | ENSG00000100099 | HPS4 | 26443423 | 26483837 | protein_coding

#### Example command

```
Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/01_run_MatrixeQTL.R \
-d ~/TOPMed_MESA_crosspop_portability/sample_data/dosages/GEUVADIS_GBR_chr22_dosage_filtered.txt.gz \
-e ~/TOPMed_MESA_crosspop_portability/sample_data/exp_tbls/GEUVADIS_GBR_expression.txt.gz \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-t GBR_chr22_cis \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/cis_ES/ \
-w 1000000
```

### Output files:

The first step will produce output files (one per population/chromosome) that are space-separated and contain 7 columns. The first 6 columns are default Matrix eQTL output columns, and the 7th contains the standard error (required to run MASHR).

snps | gene | statistic | pvalue | FDR | beta | SE
--- | --- | --- | --- | --- | --- | ---
chr22:26379178 | ENSG00000261188 | 4.17859262847936 | 7.14867054260006e-05 | 0.493258267439404 | 1.85382924858658 | 0.443649193259887
chr22:26358282 | ENSG00000261188 | 3.56862236345645 | 0.00059619153356162 | 0.703727948206992 | 0.373151478008807 | 0.104564574226169

___

### Step 2. Prepare MASHR inputs

The second step consists of preparing the input files required to run MASHR. For this, use `./DIY_scripts/02_prepare_MASHR_inputs.R`. We recommend running this step in parallel for each chromosome. The arguments to run the second step of the pipeline are:

* `-l`, `--listinputs`: Path to file containing inputs to be read for each condition. 
* `-g`, `--geneannotation`: Path to the gene annotation file. 
* `-c`, `--chromosome`: Chromosome number to be analyzed. 
* `-o`, `--output`: Output directory path. 

### Input files:
1. List of inputs: it is a space- or tab-separated file, containing 3 columns and no header. For each row, provide the condition code, the path to the first step's output, and the path to the dosage file (the same used for step 1). 

GBR | ~/TOPMed_MESA_crosspop_portability/sample_data/cis_ES/GBR_chr22_cis.txt | ~/TOPMed_MESA_crosspop_portability/sample_data/dosages/GEUVADIS_GBR_chr22_dosage_filtered.txt.gz
--- | --- | ---
YRI | ~/TOPMed_MESA_crosspop_portability/sample_data/cis_ES/YRI_chr22_cis.txt | ~/TOPMed_MESA_crosspop_portability/sample_data/dosages/GEUVADIS_YRI_chr22_dosage_filtered.txt.gz

2. Gene annotation file: it is a space- or tab-separated file containing 6 columns. Each row is a different gene and contains information about their ID, name, TSS, TES, and type. 

chr | gene_id | gene_name | start | end | gene_type
--- | --- | --- | --- | --- | ---
22 | ENSG00000128203 | ASPHD2 | 26429260 | 26445015 | protein_coding
22 | ENSG00000100099 | HPS4 | 26443423 | 26483837 | protein_coding

#### Example command

```
Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/02_prepare_MASHR_inputs.R \
-l ~/TOPMed_MESA_crosspop_portability/sample_data/chr22_files.txt \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-c 22 \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_inputs
```

### Output files:

The second step will produce output files (two per gene in the dataset) that are space-separated and contain 3+N columns, in which N is the number of populations. Output files ending in `_beta.txt` contain, as the name suggests, the effect sizes for each SNP-gene pair for each population. Effect sizes equal to 0 mean that the SNP in the respective row is missing from the population's dosage file.

snps | snp_ID | gene | GBR_beta | YRI_beta
--- | --- | --- | --- | ---
chr22:26329350 | chr22:26329350:T:C | ENSG00000100104 | -0.000942147025196598 | 0.00481675332681294
chr22:26329511 | chr22:26329511:A:C | ENSG00000100104 | 0 | 0.290500545943225

Output files ending in `_SE.txt` contain, as the name suggests, the standard errors for each SNP-gene pair for each population. Standard errors equal to 10 mean that the SNP in the respective row is missing from the population's dosage file.

snps | snp_ID | gene | GBR_SE | YRI_SE
--- | --- | --- | --- | ---
chr22:26329350 | chr22:26329350:T:C | ENSG00000100104 | -0.121966506680216 | 0.178397616871198
chr22:26329511 | chr22:26329511:A:C | ENSG00000100104 | 10 | 0.389835210194495

___

### Step 3. Run MASHR

The third step consists of adjusting the effect sizes using MASHR. For this, use `./DIY_scripts/03_run_MASHR.R`. The arguments to run the third step of the pipeline are:

* `-i`, `--input`: Path of the directory with MASHR input files made using the second step of the pipeline. 
* `-g`, `--geneannotation`: Path to the gene annotation file. 
* `-o`, `--output`: Output directory path. 

### Input files:
1. Gene annotation file: it is a space- or tab-separated file containing 6 columns. Each row is a different gene and contains information about their ID, name, TSS, TES, and type. 

chr | gene_id | gene_name | start | end | gene_type
--- | --- | --- | --- | --- | ---
22 | ENSG00000128203 | ASPHD2 | 26429260 | 26445015 | protein_coding
22 | ENSG00000100099 | HPS4 | 26443423 | 26483837 | protein_coding

#### Example command

```
Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/03_run_MASHR.R \
-i ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_inputs \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_outputs
```

### Output files:

The third step will produce output files (three per gene in the dataset) that are space-separated and contain 3+N columns, in which N is the number of populations. Output files ending in `_MASHR_beta.txt` contain, as the name suggests, the MASHR-adjusted effect sizes for each SNP-gene pair for each population. 

snps | snp_ID | gene | GBR_beta | YRI_beta
--- | --- | --- | --- | ---
chr22:26329350 | chr22:26329350:T:C | ENSG00000100104 | 1.83290685636514e-05 | 7.41946765078412e-05
chr22:26329511 | chr22:26329511:A:C | ENSG00000100104 | 0.000765449390867857 | 0.00459603564276789

Output files ending in `_MASHR_SD.txt` contain, as the name suggests, the standard deviations for each SNP-gene pair for each population. 

snps | snp_ID | gene | GBR_SD | YRI_SD
--- | --- | --- | --- | ---
chr22:26329350 | chr22:26329350:T:C | ENSG00000100104 | 0.00900736394846143 | 0.0235644643803819
chr22:26329511 | chr22:26329511:A:C | ENSG00000100104 | 0.0141805513785063 | 0.0557317748086211

Output files ending in `_MASHR_lfsr.txt` contain, as the name suggests, the local false sign rates for each SNP-gene pair for each population. Local false sign rates measure the confidence in the sign of each effect rather than confidence in each effect being non-zero [[ref](https://doi.org/10.1093/biostatistics/kxw041)]. 

snps | snp_ID | gene | GBR_lfsr | YRI_lfsr
--- | --- | --- | --- | ---
chr22:26329350 | chr22:26329350:T:C | ENSG00000100104 | 0.982997249517182 | 0.980835958151674
chr22:26329511 | chr22:26329511:A:C | ENSG00000100104 | 0.97180709286982 | 0.963355769953236

___

### Step 4. Make model files for MetaXcan

The fourth step consists of making the transcriptome prediction models in the format required by MetaXcan family tools. For this, use `./DIY_scripts/04_make_MASHR_db.R`. The arguments to run the fourth step of the pipeline are:

* `-f`, `--filesdirectory`: Path of the directory with files containing MASHR outputs generated by the third step of our pipeline. 
* `-g`, `--geneannotation`: Path to the gene annotation file. 
* `-c`, `--codes`: Conditions code used, separated by a hyphen ("-"). 
* `-o`, `--outpath`: Output directory path. 

### Input files:
1. Gene annotation file: it is a space- or tab-separated file containing 6 columns. Each row is a different gene and contains information about their ID, name, TSS, TES, and type. 

chr | gene_id | gene_name | start | end | gene_type
--- | --- | --- | --- | --- | ---
22 | ENSG00000128203 | ASPHD2 | 26429260 | 26445015 | protein_coding
22 | ENSG00000100099 | HPS4 | 26443423 | 26483837 | protein_coding

#### Example command

```
Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/04_make_MASHR_db.R \
-f ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_outputs \
-c GBR-YRI \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_models
```

### Output files:

The fourth step will produce output files (three per population). Two of them are space-separated text files, and one is a `.db` file. Output files ending in `_MASHR.db` are the files that should be used to run MetaXcan tools. Output files ending in `_MASHR_weights.txt` contain, as the name suggests, the MASHR-adjusted effect sizes of each SNP in the MASHR models. The format follows the one required by MetaXcan tools. 

gene | rsid | varID | ref_allele | eff_allele | weight
--- | --- | --- | --- | --- | ---
ENSG00000100099 | chr22:26400509 | chr22:26400509:A:G | A | G | 0.157689114274992
ENSG00000100104 | chr22:26393189 | chr22:26393189:A:G | A | G | 0.0263773669384809

Output files ending in `_MASHR_summaries.txt` contain, as the name suggests, the summaries of each gene in the MASHR models. The format follows the one required by MetaXcan tools. The last three columns are empty, as we did not implement a cross-validation step in our pipeline. This does not interfere with the usability of the models in MetaXcan. 

gene | genename | n.snps.in.model | pred.perf.R2 | pred.perf.pval | pred.perf.qval
--- | --- | --- | --- | --- | ---
ENSG00000128203 | ASPHD2 | 1 | | | 
ENSG00000100099 | HPS4 | 1 | | | 

___

### Step 5. Make model covariance files for MetaXcan

The fifth step consists of making the SNP-SNP covariance files required to run TWAS using GWAS summary statistics data with MetaXcan family tools. For this, use `./DIY_scripts/05_make_MASHR_covariances.R`. We recommend running this step in parallel for each population and chromosome, and then concatenating all outputs together per population. The arguments to run the fifth step of the pipeline are:

* `-d`, `--snpdosage`: Path to the SNP dosage file (unfiltered). 
* `-g`, `--geneannotation`: Path to the gene annotation file. 
* `-t`, `--tag`: Prefix for the output file. 
* `-c`, `--chromosome`: Chromosome number to be analyzed. 
* `-m`, `--model`: Path to model weight text file for which covariance will be computed (generated by step 4 of our pipeline). 
* `-o`, `--outputdir`: Output directory path. By default, it is `./`. 
* `-w`, `--window`: Size of the cis-SNP window from the TSS and TES. By default, it is `1e6`. 

### Input files:
1. Unfiltered SNP dosage file: it is a space- or tab-separated file, containing 5+N columns, in which N is the number of individuals in the dataset. Each row is a different SNP, and from the 6th column forward, the columns contain each sample's dosage concerning the SNP. We highly recommend removing ambiguous and multi-allelic SNPs. This file should contain all SNPs before any QC filtering, such as MAF as HWE. Each population should have its own dosage file.

chr | snp_ID | pos | ref_allele | alt_allele | Sample1
--- | --- | --- | --- | --- | ---
22 | chr22:26329350:T:C | 26329350 | T | C | 0
22 | chr22:26329511:A:C | 26329511 | A | C | 0

2. Gene annotation file: it is a space- or tab-separated file containing 6 columns. Each row is a different gene and contains information about their ID, name, TSS, TES, and type. 

chr | gene_id | gene_name | start | end | gene_type
--- | --- | --- | --- | --- | ---
22 | ENSG00000128203 | ASPHD2 | 26429260 | 26445015 | protein_coding
22 | ENSG00000100099 | HPS4 | 26443423 | 26483837 | protein_coding

3. MASHR model weights file: it is a space-separated file containing 6 columns. It is generated by step 4 of our pipeline. As the name suggests, it contains the MASHR-adjusted effect sizes of each SNP in the MASHR models. The format follows the one required by MetaXcan tools. 

gene | rsid | varID | ref_allele | eff_allele | weight
--- | --- | --- | --- | --- | ---
ENSG00000100099 | chr22:26400509 | chr22:26400509:A:G | A | G | 0.157689114274992
ENSG00000100104 | chr22:26393189 | chr22:26393189:A:G | A | G | 0.0263773669384809

#### Example command

```
Rscript ~/TOPMed_MESA_crosspop_portability/DIY_scripts/05_make_MASHR_covariances.R \
-d ~/TOPMed_MESA_crosspop_portability/sample_data/dosages/GEUVADIS_GBR_chr22_dosage_unfiltered.txt.gz \
-g ~/TOPMed_MESA_crosspop_portability/sample_data/gene_annotation.txt \
-t GBR_MASHR \
-c 22 \
-m ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_models/GBR_MASHR_weights.txt.gz \
-o ~/TOPMed_MESA_crosspop_portability/sample_data/MASHR_models \
-w 1000000
```

### Output files:

The fifth step will produce a single output file per run. It is a space-separated text file with 4 columns with the SNP-SNP covariance values in the 4th column. The format follows the one required by MetaXcan tools. It is required to run MetaXcan tools using GWAS summary statistics. 

GENE | RSID1 | RSID2|  VALUE
--- | --- | --- | ---
ENSG00000128203 | chr22:26369201 | chr22:26369201 | 0.461833105335157
ENSG00000100099 | chr22:26400509 | chr22:26400509 | 0.238030095759234
