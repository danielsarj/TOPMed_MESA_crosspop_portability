# Loading libraries and defining arguments
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(RSQLite))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
parser <- ArgumentParser()
parser$add_argument('-f', '--filesdirectory', help='path of the directory with files containing MASHR outputs')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-c', '--codes', help='conditions code used, separated by a hyphen ("-")')
parser$add_argument('-o', '--outpath', help='output directory path')
args <- parser$parse_args()

# working directory to where mashr results files are 
mashr_dir = args$filesdirectory

# Get conditions codes
codes <- args$codes %>% str_split(pattern='-') %>% unlist()

# Get gene names in the MASHR output files directory
gene_list <- list.files(mashr_dir %&% '/') %>% substr(1,15) %>% unique()

# Figure out what is the most significant SNP per gene, per pop
print('INFO: Assessing top SNPs per conditions for each gene')
for (working_gene in gene_list){
  # Read data frame containing LFSRs
  lfsr <- fread(mashr_dir %&% "/" %&% working_gene %&% '_MASHR_lfsr.txt.gz', header=T)
  
  # Initialize empty list
  lfsr_dfs <- list()
  
  # Get lfsr per condition
  for (i in 1:length(codes)){
    lfsr_dfs[[i]] <- lfsr %>% select(gene, snps, snp_ID, contains(codes[i]))
  }

  # Get top SNP per condition
  for (j in 1:length(codes)){
    tmp <- lfsr_dfs[[j]] 
    colnames(tmp)[4] <- 'lfsr'
    tmp <- tmp %>% slice(which.min(lfsr))
    if (exists('top_SNPs_df')){
      top_SNPs_df <- rbind(top_SNPs_df, tmp)
    } else {top_SNPs_df <- tmp}
  }
}

# Remove duplicates from the top SNP df
top_SNPs_df <- top_SNPs_df %>% select(-lfsr) %>% unique()

# Get MASHR-adjusted betas for the top SNPs
print('INFO: Making condition-specific transcriptome models')
for (c in codes){
  print('INFO: Current condition code is ' %&% c)
  
  # Get betas for each gene
  for (working_gene in gene_list){
    mashr_in <- fread(mashr_dir %&% "/" %&% working_gene %&% '_MASHR_beta.txt.gz') %>% select(gene, snps, snp_ID, contains(c)) %>% 
      inner_join(top_SNPs_df, by=c('gene', 'snps', 'snp_ID'))
    if (exists('weights_df')){
      weights_df <- rbind(weights_df, mashr_in)
    } else {weights_df <- mashr_in}
  }
  
  # Make column with REF and ALT alleles
  weights_df <- weights_df %>% 
    mutate(refAllele=substr(snp_ID, nchar(snp_ID)-2, nchar(snp_ID)-2), effectAllele=substr(snp_ID, nchar(snp_ID), nchar(snp_ID))) %>% 
    rename(beta=contains(c)) %>% filter(beta!=0)
  
  # Get weight data frame into PrediXcan format
  weights_df <- weights_df %>% select('gene','snps','snp_ID','refAllele','effectAllele','beta') %>%
    rename(gene=gene, rsid=snps, varID=snp_ID, ref_allele=refAllele, eff_allele=effectAllele, weight=beta) 
  
  # Get number of SNPs per model 
  genes_table <- table(weights_df$gene) %>% as.data.frame() %>% rename(gene=Var1, n_snps=Freq)
  
  # Make summary table for df file
  model_summaries <- fread(args$geneannotation, header=T, stringsAsFactors=F) %>% 
    select(gene_id, gene_name) %>% unique() %>% filter(gene_id %in% gene_list) %>% 
    inner_join(genes_table, by=c('gene_id'='gene'))
  model_summaries$rho_avg_squared <- rep(NA, nrow(model_summaries))
  model_summaries$zscore_pval <- rep(NA, nrow(model_summaries))
  model_summaries$zscore_qval <- rep(NA, nrow(model_summaries))
  model_summaries <- model_summaries %>% rename(gene=gene_id, genename=gene_name, n.snps.in.model=n_snps, pred.perf.R2=rho_avg_squared,
                            pred.perf.pval=zscore_pval, pred.perf.qval=zscore_qval)
  
  # Make final files
  fwrite(model_summaries, args$outpath %&% '/' %&% c %&%'_MASHR_summaries.txt', col.names=T, quote=F, sep=' ')
  fwrite(weights_df, args$outpath %&% '/' %&% c %&%'_MASHR_weights.txt', col.names=T, quote=F, sep=' ')
  conn <- dbConnect(drv = driver, args$outpath %&% '/' %&% c %&%'_MASHR.db')
  dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
  dbExecute(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbWriteTable(conn, 'weights', weights_df, overwrite = TRUE)
  dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbDisconnect(conn)
  rm(weights_df)
  print('INFO: Successfully made transcriptome prediction model files for '%&% c)
}
