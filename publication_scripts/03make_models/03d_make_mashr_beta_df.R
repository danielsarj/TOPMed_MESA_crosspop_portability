suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep='')
setwd('/home/daniel/mashr/mashr_dfs/WGS_files/outputs')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
args <- parser$parse_args()

sigSNPs <- fread('top_' %&% args$tissue %&% '_mashr_sigSNPs_shared_genes.txt.gz')
gene_list <- list.files(args$tissue,pattern='allSNPs') %>% substr(1,15) %>% unique()

for (gene in gene_list){
  tmp <- fread(args$tissue %&% '/' %&% gene  %&% '_allSNPs_mashr_beta.txt.gz') %>% inner_join(sigSNPs)
  if (exists('final.df')){
    final.df <- rbind(final.df, tmp)
  } else {final.df <- tmp}
}

fwrite(final.df, '/home/daniel/mashr/mashr_db/WGS_files/top_' %&% args$tissue %&% '_mashr_sigSNPs-betas_shared_genes.txt', sep='\t', quote=F)

if (args$tissue=='PBMC'){
  final.df <- final.df %>% filter(!(AFA_beta==0 & EUR_beta==0 & HIS_beta==0 & CHN_beta==0))
} else {final.df <- final.df %>% filter(!(AFA_beta==0 & EUR_beta==0 & HIS_beta==0))}

fwrite(final.df, '/home/daniel/mashr/mashr_db/WGS_files/top_' %&% args$tissue %&% '_mashr_nonzerosSNPs-betas_shared_genes.txt', sep='\t', quote=F)

final.table <- final.df %>% select(gene, snps)
final.table <- table(final.table$gene) %>% reshape2::melt()
colnames(final.table) <- c('gene','n_snps')

fwrite(final.table, '/home/daniel/mashr/mashr_db/WGS_files/top_' %&% args$tissue %&% '_mashr_nonzerosSNPs-betas_shared_genes_table.txt', sep='\t', quote=F)