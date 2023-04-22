suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep='')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
parser$add_argument('--pop', help='population whose files will be analyzed')
parser$add_argument('--pcs', help='number of PCs in file name')
args <- parser$parse_args()

page <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/PAGE_SNPlist.txt.gz') %>% pull(SNP_hg38)
snp_location_snps <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/' %&% args$pop %&% '.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz', header=T) %>% 
  filter(!grepl("A:T", varID)) %>% filter(!grepl("T:A", varID)) %>% filter(!grepl("C:G", varID)) %>% filter(!grepl("G:C", varID)) %>% pull(rsid) # removes ambiguous variants

matrixeqtl <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/output/cis_eQTLs_' %&% args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs %&% '_WG_all_cis.txt.gz', header=T)
matrixeqtl_new <- matrixeqtl %>% filter(snps %in% snp_location_snps) %>% filter(snps %in% page)
matrixeqtl_new$SE <- (matrixeqtl_new$beta/abs(matrixeqtl_new$statistic)) #https://stats.stackexchange.com/questions/337070/compute-standard-error-from-beta-p-value-sample-size-and-the-number-of-regres - compute SE which is necessary for mashr

fwrite(matrixeqtl_new, '/home/daniel/mashr/matrixeQTL/WGS_files/output/cis_eQTLs_' %&% args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs %&% '_WG_noDup.wSE_cis_PAGEintersect-new.txt', sep='\t', quote=F)
system('gzip /home/daniel/mashr/matrixeQTL/WGS_files/output/cis_eQTLs_' %&% args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs %&% '_WG_noDup.wSE_cis_PAGEintersect-new.txt')