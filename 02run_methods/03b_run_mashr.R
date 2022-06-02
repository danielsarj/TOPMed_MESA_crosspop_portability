suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(argparse))
suppressMessages(library(mashr))
'%&%' = function(a,b) paste (a,b,sep='')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
parser$add_argument('--mode', help='mode used to prepare mashr input dataframes. expects either sigSNPs or allSNPs')
parser$add_argument('--gfile', help='gene text file number. expects a number between 1 and 10')
args <- parser$parse_args()

setwd('/home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue)
gene_list <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/output/top_' %&% args$tissue %&% '_sigSNPs_shared_genes_table_split' %&% args$gfile %&% '-new.txt.gz') %>% pull(gene)

for (gene in gene_list){
  print(gene)
  beta <- fread(gene %&% '_' %&% args$mode %&% '_beta.txt.gz', header=T, stringsAsFactors=F) %>% select(contains('beta')) %>% as.matrix()
  se <- fread(gene %&% '_' %&% args$mode %&% '_SE.txt.gz', header=T, stringsAsFactors=F) %>% select(contains('SE')) %>% as.matrix() %>% abs() 
  #use abs() as recommended by mashr:
  #'Both Bhat and Shat are zero (or near zero) for some input data. Please check your input. 
  #If it is expected please set Shat to a positive number to avoid numerical issues;'
  df <- fread(gene %&% '_' %&% args$mode %&% '_beta.txt.gz', header=T, stringsAsFactors=F) %>% select(gene, snps, varID)
  
  if (nrow(df)>1){
  #set up main data object
  data = mash_set_data(beta, se)
  
  #covariance matrixes
  data.c = cov_canonical(data) #canonical
  data.pca = cov_pca(data,min(ncol(beta),nrow(beta))) #pca
  data.ed = cov_ed(data, data.pca) #data-driven
  
  #fit model
  m = mash(data, Ulist=c(data.ed,data.c)) #also computes posterior summaries
  
  #posterior summaries
  posterior_lfsr <- get_lfsr(m)
  posterior_lfsr <- cbind(df, posterior_lfsr)
  
  posterior_mean <- get_pm(m)
  posterior_mean <- cbind(df, posterior_mean)

  posterior_sd <- get_psd(m)
  posterior_sd <- cbind(df, posterior_sd)
  
  #write output
  fwrite(posterior_lfsr,file='/home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue %&% '/' %&% gene %&% '_' %&% args$mode %&% '_mashr_lfsr.txt',quote=F,sep='\t')
  fwrite(posterior_mean,file='/home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue %&% '/' %&% gene %&% '_' %&% args$mode %&% '_mashr_beta.txt',quote=F,sep='\t')
  fwrite(posterior_sd,file='/home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue %&% '/' %&% gene %&% '_' %&% args$mode %&% '_mashr_SD.txt',quote=F,sep='\t')
  system('gzip /home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue %&% '/' %&% gene %&% '_' %&% args$mode %&% '_mashr_lfsr.txt')
  system('gzip /home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue %&% '/' %&% gene %&% '_' %&% args$mode %&% '_mashr_beta.txt')
  system('gzip /home/daniel/mashr/mashr_dfs/WGS_files/outputs/' %&% args$tissue %&% '/' %&% gene %&% '_' %&% args$mode %&% '_mashr_SD.txt')
  
  }
}