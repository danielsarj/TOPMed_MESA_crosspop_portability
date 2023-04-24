# Loading libraries and defining arguments
suppressMessages(library(tidyverse))
suppressMessages(library(RSQLite))
suppressMessages(library(data.table))
suppressMessages(library(readr))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
parser <- ArgumentParser()
parser$add_argument('-f', '--filesdirectory', help='path of the directory with files containing MASHR outputs')
parser$add_argument('-c', '--codes', help='conditions code used, separated by a hyphen ("-")')
parser$add_argument('--pop', help='population whose files will be analyzed')
args <- parser$parse_args()

# Change working directory to where input files are 
setwd('/home/daniel/github_mashr_project/sample_data/MASHR_outputs')

# Get conditions codes
codes <- c('GBR-YRI') %>% str_split(pattern='-') %>% unlist()

# Get gene names in the MASHR output files directory
gene_list <- list.files('.') %>% substr(1,15) %>% unique()

# Figure out what is the most significant SNP per gene, per pop
for (working_gene in gene_list){
  print('INFO: Assessing top SNPs for ' %&% working_gene)
  
  # Read data frame containing LFSRs
  lfsr <- fread(working_gene %&% '_MASHR_lfsr.txt.gz', header=T)
  
  # Initialize empty list
  lfsr_dfs <- list()
  
  for (i in 1:length(codes)){
    lfsr_dfs[[i]] <- lfsr %>% select(gene, snps, snp_ID, contains(codes[i]))
    
  }
  
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



out_prefix<- '/home/daniel/mashr/mashr_db/WGS_files/mashr_models/' %&% args$tissue %&% '_' %&% args$pop %&% '_mashr_baseline'
#sigSNPs <- fread('top_' %&% args$tissue %&% '_mashr_sigSNPs_shared_genes.txt.gz')
SNPs_anno <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/MetaPop.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz')
#sigSNPs_g <- fread('top_' %&% tissue %&% '_mashr_allSNPs_shared_genes_table.txt.gz')

### reads mashr output file
#reading input
mashr_in <- fread('/home/daniel/mashr/mashr_db/WGS_files/top_' %&% args$tissue %&% '_mashr_sigSNPs-betas_shared_genes.txt.gz', stringsAsFactors=F) %>% 
  select(gene, snps, varID, contains(args$pop)) %>% rename(beta=contains(args$pop)) %>%
  left_join(SNPs_anno) %>% arrange(gene)

### making final files
#weights table
weights <- mashr_in %>% select('gene','snps','varID','refAllele','effectAllele','beta') %>% unique() %>% na_if(0) %>% drop_na()
weights <- rename(weights, 
                  gene=gene,
                  rsid=snps,
                  varID=varID,
                  ref_allele=refAllele, 
                  eff_allele=effectAllele,
                  weight=beta)
genes_list <- weights %>% pull(gene) %>% unique()
genes_table <- table(weights$gene) %>% as.data.frame() %>% rename(gene=Var1, n_snps=Freq)

#extra table
model_summaries <- select(read.table('/home/daniel/gencode_annotation/gene_annotation_v38_60649.txt', header=T, stringsAsFactors=F), c('gene_id','gene_name')) %>% unique() %>% filter(gene_id %in% genes_list) %>% inner_join(genes_table, by=c('gene_id'='gene'))
model_summaries$rho_avg_squared <- rep(NA, nrow(model_summaries))
model_summaries$zscore_pval <- rep(NA, nrow(model_summaries))
model_summaries$zscore_qval <- rep(NA, nrow(model_summaries))
model_summaries <- rename(model_summaries, 
                          gene=gene_id, 
                          genename=gene_name,
                          n.snps.in.model=n_snps,
                          pred.perf.R2=rho_avg_squared,
                          pred.perf.pval=zscore_pval,
                          pred.perf.qval=zscore_qval)
#final files
fwrite(model_summaries, out_prefix %&% '_summaries.txt', col.names=T, quote=F, sep=' ')
fwrite(weights, out_prefix %&% '_weights.txt', col.names=T, quote=F, sep=' ')
conn <- dbConnect(drv = driver, out_prefix %&% '.db')
dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbDisconnect(conn)