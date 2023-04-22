# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/06Elastic_net

setwd('/home/daraujo1/scratch/TOPMed_MESA/Elastic_Net')
source('01a_WGS_EN_baseline.R')
'%&%' <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly=TRUE)
tissue <- argv[1]
pop <- argv[2]
chrom <- argv[3]

snp_annot_file <- '/home/daraujo1/data/TOPMed_MESA/EN_inputs/' %&% tissue %&% '/annotation_files/' %&% pop %&% '/' %&% tissue %&% '.' %&% pop %&% '.' %&% chrom %&% '.annotation_nodup.txt.gz' #no dup snps!
gene_annot_file <- '/home/daraujo1/data/TOPMed_MESA/EN_inputs/' %&% tissue %&% '/annotation_files/annotation_all_chr_genes_ENSG_mod.txt' #file with all genes
genotype_file <- '/home/daraujo1/data/TOPMed_MESA/EN_inputs/' %&% tissue %&% '/dosage_files/' %&% pop %&% '/' %&% tissue %&% '.' %&% pop %&% '.' %&% chrom %&% '.dosage_nodup_v3.txt.gz' #no dup snps!
expression_file <- '/home/daraujo1/data/TOPMed_MESA/EN_inputs/' %&% tissue %&% '/adjusted_expression/gene_level/chr' %&% chrom %&% '/RNASeq_TOPMed_' %&% pop %&% '_ln_adjAgeSex_mean_rank-inverse_adj10PCs_MEQTL_sorted_geneidsupdated.txt.gz' #expresion matrix 
prefix <- '/home/daraujo1/data/TOPMed_MESA/EN_outputs/' %&% tissue %&% '/' %&% pop %&% '/' %&% tissue %&% '_' %&% pop %&% '_chr' %&% chrom %&% '_base' #done

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)

