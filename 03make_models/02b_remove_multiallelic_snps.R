# loading libraries
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep='')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
parser$add_argument('--pop', help='population whose files will be analyzed')
parser$add_argument('--pcs', help='number of PCs in file name')
args <- parser$parse_args()

print(args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs)

matrixeqtl <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/output/cis_eQTLs_' %&% args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs %&% '_WG_all_cis.txt.gz', header=T)
matrixeqtl_unique_snps <- matrixeqtl %>% select(snps, gene) 
matrixeqtl_unique_snps <- matrixeqtl_unique_snps[!(duplicated(matrixeqtl_unique_snps) | duplicated(matrixeqtl_unique_snps, fromLast=T)),] # removes multiallelic snps 
matrixeqtl_new <- left_join(matrixeqtl_unique_snps, matrixeqtl)
fwrite(matrixeqtl_new, '/home/daniel/mashr/matrixeQTL/WGS_files/output/cis_eQTLs_' %&% args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs %&% '_WG_noDup_cis.txt', sep='\t', quote=F)
system('gzip /home/daniel/mashr/matrixeQTL/WGS_files/output/cis_eQTLs_' %&% args$tissue %&% '_' %&% args$pop %&% '_PCs_' %&% args$pcs %&% '_WG_noDup_cis.txt')

#snp_location <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/' %&% args$pop %&% '.' %&% args$tissue %&% '.WG_locations.txt.gz', header=T)
#snp_location_unique_snps <- snp_location %>% select(rsid) 
#snp_location_unique_snps <- snp_location_unique_snps[!(duplicated(snp_location_unique_snps) | duplicated(snp_location_unique_snps, fromLast=T)),]
#snp_location_new <- left_join(snp_location_unique_snps, snp_location)
#fwrite(snp_location_new, '/home/daniel/mashr/matrixeQTL/WGS_files/' %&% args$pop %&% '.' %&% args$tissue %&% '.WG_noDup_locations.txt', sep=' ', quote=F)