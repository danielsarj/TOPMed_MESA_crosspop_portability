suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(reshape2))
"%&%" <- function(a,b) paste(a,b, sep='')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
args <- parser$parse_args()

setwd('/home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue)
gene_list <- list.files('.',pattern='allSNPs') %>% substr(1,15) %>% unique()

for (gene in gene_list){
  print(gene)
  lfsr <- fread(gene  %&% '_allSNPs_mashr_lfsr.txt.gz')
  afa_lfsr <- lfsr %>% select(gene, snps, varID, AFA_beta) %>% rename(lfsr=AFA_beta)
  eur_lfsr <- lfsr %>% select(gene, snps, varID, EUR_beta) %>% rename(lfsr=EUR_beta)
  his_lfsr <- lfsr %>% select(gene, snps, varID, HIS_beta) %>% rename(lfsr=HIS_beta)
  if (args$tissue=='PBMC'){
    chn_lfsr <- lfsr %>% select(gene, snps, varID, CHN_beta) %>% rename(lfsr=CHN_beta)
    df_list <- list(afa_lfsr, eur_lfsr, his_lfsr, chn_lfsr)
  } else {
    df_list <- list(afa_lfsr, eur_lfsr, his_lfsr)
  }
  
  #get top SNP per pop
  for (df in df_list){
    tmp <- df %>% slice(which.min(lfsr))
    if (exists('final.df')){
      final.df <- rbind(final.df, tmp)
    } else {final.df <- tmp}
  }
}

final.df <- final.df %>% select(-lfsr) %>% arrange(gene) %>% unique()
table <- table(final.df$gene) %>% reshape2::melt()
colnames(table) <- c('gene','n_snps')
fwrite(final.df, '/home/daniel/mashr/mashr_dfs/WGS_files/outputs/top_' %&% args$tissue %&% '_mashr_sigSNPs_shared_genes.txt', sep='\t')
fwrite(table, '/home/daniel/mashr/mashr_dfs/WGS_files/outputs/top_' %&% args$tissue %&% '_mashr_sigSNPs_shared_genes_table.txt', sep='\t')