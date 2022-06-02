suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
"%&%" <- function(a,b) paste(a,b, sep='')
setwd('/home/daniel/mashr/matrixeQTL/WGS_files/output')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
args <- parser$parse_args()

#loading matrixeQTL files
afa <- fread('cis_eQTLs_' %&% args$tissue %&% '_AFA_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz') %>%
  inner_join(fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/AFA.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz'), by=c('snps'='rsid'))
eur <- fread('cis_eQTLs_' %&% args$tissue %&% '_EUR_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz') %>%
  inner_join(fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/EUR.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz'), by=c('snps'='rsid'))
his <- fread('cis_eQTLs_' %&% args$tissue %&% '_HIS_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz') %>%
  inner_join(fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/HIS.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz'), by=c('snps'='rsid'))
if (args$tissue=='PBMC'){
  chn <- fread('cis_eQTLs_PBMC_CHN_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz') %>%
    inner_join(fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/CHN.PBMC.WG_noDup_locations.txt.gz'), by=c('snps'='rsid'))
  df_list <- list(afa, eur, chn, his)
} else {df_list <- list(afa, eur, his)} 

#get top SNP per pop
for (df in df_list){
  tmp <- df %>% group_by(gene) %>% slice(which.min(pvalue))
  tmp <- tmp %>% select(gene, snps, varID, refAllele, effectAllele)
  if (exists('final.df')){
    final.df <- rbind(final.df, tmp)
  } else {final.df <- tmp}
}

#making final output files
final.df <- final.df %>% arrange(gene) %>% unique()
table <- table(final.df$gene) %>% reshape2::melt()
colnames(table) <- c('gene','n_snps')

afa_g <- afa %>% select(gene) %>% unique()
eur_g <- eur %>% select(gene) %>% unique()
his_g <- his %>% select(gene) %>% unique()
if (args$tissue=='PBMC'){
  chn_g <- chn %>% select(gene) %>% unique()
  shared_g <- inner_join(afa_g, eur_g) %>% inner_join(chn_g) %>% inner_join(his_g)
  full_g <- full_join(afa_g, eur_g) %>% full_join(chn_g) %>% full_join(his_g)
} else {
  shared_g <- inner_join(afa_g, eur_g) %>% inner_join(his_g)
  full_g <- full_join(afa_g, eur_g) %>% full_join(his_g)
} 

final.df.shared <- inner_join(final.df, shared_g)
table.shared <- table(final.df.shared$gene) %>% reshape2::melt()
colnames(table.shared) <- c('gene','n_snps')

fwrite(final.df, 'top_' %&% args$tissue %&% '_sigSNPs_all_genes-new.txt', sep='\t')
fwrite(final.df.shared, 'top_' %&% args$tissue %&% '_sigSNPs_shared_genes-new.txt', sep='\t')
fwrite(table, 'top_' %&% args$tissue %&% '_sigSNPs_all_genes_table-new.txt', sep='\t')
fwrite(table.shared, 'top_' %&% args$tissue %&% '_sigSNPs_shared_genes_table-new.txt', sep='\t')
#fwrite(afa, 'cis_eQTLs_' %&% args$tissue %&% '_AFA_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt', sep='\t', quote=F)
#fwrite(eur, 'cis_eQTLs_' %&% args$tissue %&% '_EUR_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt', sep='\t', quote=F)
#fwrite(his, 'cis_eQTLs_' %&% args$tissue %&% '_HIS_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt', sep='\t', quote=F)
#if (args$tissue=='PBMC'){
#  fwrite(chn, 'cis_eQTLs_PBMC_CHN_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt', sep='\t', quote=F)
#}