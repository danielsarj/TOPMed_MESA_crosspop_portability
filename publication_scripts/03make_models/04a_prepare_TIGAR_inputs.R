library(tidyverse)
library(data.table)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep='')
setwd('/home/daniel/mashr/TIGAR/TOPMed_MESA')

# arguments
parser <- ArgumentParser()
parser$add_argument("--tis", help="MESA tissue")
parser$add_argument("--pop", help="MESA population")
args <- parser$parse_args()

# make expression tables & id files
if (args$pop != 'EUR'){
  init_exp <- fread('/home/chris/topmed_expression_whole_genome/expression/'%&% args$tis %&%'_expression_'%&% args$pop %&%'age_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10_transposed.txt',
                  header=T)
  } else {
    init_exp <- fread('/home/chris/topmed_expression_whole_genome/expression/'%&% args$tis %&%'_expression_CAUage_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10_transposed.txt',
                      header=T)
  }
gene_anno <- fread('/home/daniel/gencode_annotation/gene_annotation_v38_60649.txt') %>% select(-gene_type) %>% filter(chr!='M', chr!='Y')
gene_anno$chr <- as.numeric(gene_anno$chr)
full_exp <- inner_join(gene_anno, init_exp, by=c('gene_id'))
colnames(full_exp)[1:5] <- c('CHROM', 'TargetID', 'GeneName', 'GeneStart', 'GeneEnd')
full_exp <- full_exp %>% select('CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName', everything())
ids <- colnames(full_exp)[6:ncol(full_exp)] %>% as.data.frame()
fwrite(full_exp, 'inputs/exp/'%&% args$tis %&%'_'%&% args$pop %&%'_TIGAR_expression_tbl.txt', col.names=T, quote=F, sep='\t')
fwrite(ids, 'inputs/ids/'%&% args$tis %&%'_'%&% args$pop %&%'_TIGAR_sampleIDs.txt', col.names=F, quote=F, sep='\t')


# making dosage file per chr
for (chrom in c(1:23, 25)){
  ids <- fread('inputs/ids/'%&% args$tis %&%'_'%&% args$pop %&%'_TIGAR_sampleIDs.txt') %>% pull() %>% as.character()
  genotype_file <- fread('/home/daniel/MESA_genotypes_subset/WGS_files/dosages/'%&% args$pop %&%'.'%&% args$tis %&%'.chr'%&% chrom %&%'_genotype_uniq.page.txt.gz', header=T)
  genotype_anno <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/' %&% args$pop %&% '.' %&% args$tis %&% '.WG_noDup_locations.txt.gz', header=T) %>%
    filter(chr==chrom) %>% filter(!grepl("A:T", varID)) %>% filter(!grepl("T:A", varID)) %>% filter(!grepl("C:G", varID)) %>% filter(!grepl("G:C", varID))
  dosage_file <- inner_join(genotype_anno, genotype_file, by=c('rsid'='snp_ID'))
  colnames(dosage_file)[1:6] <- c('ID', 'CHROM', 'POS', 'varID', 'REF', 'ALT')
  dosage_file <- dosage_file %>% select(-varID) %>% select(CHROM, POS, ID, REF, ALT, all_of(ids))

  fwrite(dosage_file, 'inputs/dosage/'%&% args$tis %&%'_'%&% args$pop %&%'_chr' %&% chrom %&% '_TIGAR_dosage.txt', col.names=T, quote=F, sep='\t')
  system('bgzip inputs/dosage/'%&% args$tis %&%'_'%&% args$pop %&%'_chr' %&% chrom %&% '_TIGAR_dosage.txt')
  system('tabix -f -S1 -s1 -b2 -e2 inputs/dosage/'%&% args$tis %&%'_'%&% args$pop %&%'_chr' %&% chrom %&% '_TIGAR_dosage.txt.gz')
  
  # making TIGAR command line
  command <- c('../TIGAR_Model_Train.sh --model DPR --gene_exp inputs/exp/'%&% args$tis %&%'_'%&% args$pop %&%'_TIGAR_expression_tbl.txt --train_sampleID inputs/ids/'%&% args$tis %&%'_'%&% args$pop %&%'_TIGAR_sampleIDs.txt --chr ' %&% chrom %&% ' --genofile inputs/dosage/'%&% args$tis %&%'_'%&% args$pop %&%'_chr' %&% chrom %&% '_TIGAR_dosage.txt.gz --genofile_type dosage --maf 0.01 --hwe 0.000001 --cvR2 0 --dpr 1 --ES fixed --thread 1 --out_dir outputs/'%&% args$tis %&%'/'%&% args$pop %&%' --TIGAR_dir ../')
  system('bash ' %&% command)
}
