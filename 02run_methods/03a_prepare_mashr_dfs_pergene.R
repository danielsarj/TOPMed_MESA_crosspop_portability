suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(argparse))
setwd('/home/daniel/mashr/matrixeQTL/WGS_files/output')
'%&%' = function(a,b) paste (a,b,sep='')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
parser$add_argument('--mode', help='way to prepare mashr input dataframes. either sigSNPs to build only using the sigSNPs per populations, or allSNPs to use all SNPs per gene across all pops')
parser$add_argument('--gfile', help='gene text file number. expects a number between 1 and 10')
args <- parser$parse_args()

#if (args$mode=='sigSNPs'){
  #loading input files
#  sigSNPs <- fread('top_' %&% args$tissue %&% '_sigSNPs_shared_genes.txt.gz')
#  AFA_matrixeqtl <- fread('cis_eQTLs_' %&% args$tissue %&% '_AFA_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt.gz') %>% right_join(sigSNPs) %>% select(gene, snps, varID, beta, SE)
#  EUR_matrixeqtl <- fread('cis_eQTLs_' %&% args$tissue %&% '_EUR_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt.gz') %>% right_join(sigSNPs) %>% select(gene, snps, varID, beta, SE)
#  HIS_matrixeqtl <- fread('cis_eQTLs_' %&% args$tissue %&% '_HIS_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt.gz') %>% right_join(sigSNPs) %>% select(gene, snps, varID, beta, SE)
#  if (args$tissue=='PBMC'){
#    CHN_matrixeqtl <- fread('cis_eQTLs_PBMC_CHN_PCs_10_WG_noDup.wSE_cis_PAGEintersect.txt.gz') %>% right_join(sigSNPs) %>% select(gene, snps, varID, beta, SE)
#  } 

#  #list of genes
#  gene_list <- sigSNPs %>% pull(gene) %>% unique()
  
#  #df to keep track of wrong things lol
#  something_wrong.df=data.frame(gene=as.character(), tissue=as.character(), time=as.character())
  
  #take a gene from the list, gets betas and SEs for all pops
  #and writes the output files
#  for (genes in gene_list){
#    cat('Making ' %&% args$tissue %&% ' - ' %&% genes %&% ' input data frames\n')
#    AFA_filt <- filter(AFA_matrixeqtl, gene == genes) %>% rename(AFA_beta=beta, AFA_SE=SE)
#    EUR_filt <- filter(EUR_matrixeqtl, gene == genes) %>% rename(EUR_beta=beta, EUR_SE=SE)
#    HIS_filt <- filter(HIS_matrixeqtl, gene == genes) %>% rename(HIS_beta=beta, HIS_SE=SE)
#    if (args$tissue=='PBMC'){
#      CHN_filt <- filter(CHN_matrixeqtl, gene == genes) %>% rename(CHN_beta=beta, CHN_SE=SE)
#      allpops_beta <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% full_join(CHN_filt) %>% select(-contains('SE')) %>% arrange(snps)
#     allpops_se <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% full_join(CHN_filt) %>% select(-contains('beta')) %>% arrange(snps)
#   } else {
#     allpops_beta <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% select(-contains('SE')) %>% arrange(snps)
#     allpops_se <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% select(-contains('beta')) %>% arrange(snps)
#   }
    
    #writes final files
#   snp_list <- allpops_beta %>% pull(varID) %>% unique()
#   if (nrow(allpops_beta)==nrow(allpops_se) && nrow(allpops_beta)==length(snp_list)){
#     fwrite(allpops_beta, '/home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue %&% '/' %&% genes %&% '_' %&% args$mode %&% '_beta.txt', quote=F, sep='\t', na=0)
#     fwrite(allpops_se, '/home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue %&% '/' %&% genes %&% '_' %&% args$mode %&% '_SE.txt', quote=F, sep='\t', na=10)
#   } else {
#     something_wrong.df$gene <- genes
#     something_wrong.df$tissue <- args$tissue
#     something_wrong.df$time <- Sys.time()
#   }
# }
#}

if (args$mode=='allSNPs'){
  #loading input files
  AFA_matrixeqtl <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/AFA.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz', header=T) %>% 
    filter(!grepl("A:T", varID)) %>% filter(!grepl("T:A", varID)) %>% filter(!grepl("C:G", varID)) %>% filter(!grepl("G:C", varID)) %>% select(rsid, varID) %>% 
    rename(snps=rsid) %>% unique() %>% inner_join(fread('cis_eQTLs_' %&% args$tissue %&% '_AFA_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz'), by=c('snps'='snps')) %>% arrange(gene) %>% select(gene, snps, varID, beta, SE) %>% unique()
  EUR_matrixeqtl <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/EUR.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz', header=T) %>% 
    filter(!grepl("A:T", varID)) %>% filter(!grepl("T:A", varID)) %>% filter(!grepl("C:G", varID)) %>% filter(!grepl("G:C", varID)) %>% select(rsid, varID) %>% 
    rename(snps=rsid) %>% unique() %>% inner_join(fread('cis_eQTLs_' %&% args$tissue %&% '_EUR_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz'), by=c('snps'='snps')) %>% arrange(gene) %>% select(gene, snps, varID, beta, SE) %>% unique()
  HIS_matrixeqtl <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/HIS.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz', header=T) %>% 
    filter(!grepl("A:T", varID)) %>% filter(!grepl("T:A", varID)) %>% filter(!grepl("C:G", varID)) %>% filter(!grepl("G:C", varID)) %>% select(rsid, varID) %>% 
    rename(snps=rsid) %>% unique() %>% inner_join(fread('cis_eQTLs_' %&% args$tissue %&% '_HIS_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz'), by=c('snps'='snps')) %>% arrange(gene) %>% select(gene, snps, varID, beta, SE) %>% unique()
  if (args$tissue=='PBMC'){
    CHN_matrixeqtl <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/CHN.PBMC.WG_noDup_locations.txt.gz', header=T) %>% 
      filter(!grepl("A:T", varID)) %>% filter(!grepl("T:A", varID)) %>% filter(!grepl("C:G", varID)) %>% filter(!grepl("G:C", varID)) %>% select(rsid, varID) %>% 
      rename(snps=rsid) %>% unique() %>% inner_join(fread('cis_eQTLs_PBMC_CHN_PCs_10_WG_noDup.wSE_cis_PAGEintersect-new.txt.gz'), by=c('snps'='snps')) %>% arrange(gene) %>% select(gene, snps, varID, beta, SE) %>% unique()
  }  
  
  #list of genes
  gene_list <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/output/top_' %&% args$tissue %&% '_sigSNPs_shared_genes_table_split' %&% args$gfile %&% '-new.txt.gz') %>% pull(gene) %>% unique()  
  
  #df to keep track of wrong things lol
  something_wrong.df=data.frame(gene=as.character(), tissue=as.character(), time=as.character())
  
  #take a gene from the list, gets betas and SEs for all pops
  #and writes the output files
  for (genes in gene_list){
    cat('Making ' %&% args$tissue %&% ' - ' %&% genes %&% ' input data frames\n')
    AFA_filt <- filter(AFA_matrixeqtl, gene == genes) %>% rename(AFA_beta=beta, AFA_SE=SE)
    EUR_filt <- filter(EUR_matrixeqtl, gene == genes) %>% rename(EUR_beta=beta, EUR_SE=SE)
    HIS_filt <- filter(HIS_matrixeqtl, gene == genes) %>% rename(HIS_beta=beta, HIS_SE=SE)
    if (args$tissue=='PBMC'){
      CHN_filt <- filter(CHN_matrixeqtl, gene == genes) %>% rename(CHN_beta=beta, CHN_SE=SE)
      allpops_beta <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% full_join(CHN_filt) %>% select(-contains('SE')) %>% arrange(snps)
      allpops_se <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% full_join(CHN_filt) %>% select(-contains('beta')) %>% arrange(snps)
    } else {
      allpops_beta <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% select(-contains('SE')) %>% arrange(snps)
      allpops_se <- full_join(AFA_filt, EUR_filt) %>% full_join(HIS_filt) %>% select(-contains('beta')) %>% arrange(snps)
    }
    
    #writes final files
    snp_list <- allpops_beta %>% pull(varID) %>% unique()
    if (nrow(allpops_beta)==nrow(allpops_se) && nrow(allpops_beta)==length(snp_list)){
    fwrite(allpops_beta, '/home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue %&% '/' %&% genes %&% '_' %&% args$mode %&% '_beta.txt', quote=F, sep='\t', na=0)
    fwrite(allpops_se, '/home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue %&% '/' %&% genes %&% '_' %&% args$mode %&% '_SE.txt', quote=F, sep='\t', na=10)
    system('gzip /home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue %&% '/' %&% genes %&% '_' %&% args$mode %&% '_beta.txt')
    system('gzip /home/daniel/mashr/mashr_dfs/WGS_files/inputs/' %&% args$tissue %&% '/' %&% genes %&% '_' %&% args$mode %&% '_SE.txt')
    } else {
      something_wrong.df$gene <- genes
      something_wrong.df$tissue <- args$tissue
      something_wrong.df$time <- Sys.time()
    }
  }
}
 
if (nrow(something_wrong.df)!=0){
  fwrite(something_wrong.df, '/home/daniel/mashr/mashr_dfs/WGS_files/inputs/something_wrong_df_' %&% args$mode %&% 'input-mashr.txt', append=T)
}