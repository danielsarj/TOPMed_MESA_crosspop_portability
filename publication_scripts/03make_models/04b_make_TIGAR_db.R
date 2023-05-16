library(tidyverse)
library(data.table)
library(argparse)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
setwd('/home/daniel/mashr/TIGAR/TOPMed_MESA')

# arguments
parser <- ArgumentParser()
parser$add_argument("--tis", help="MESA tissue")
parser$add_argument("--pop", help="MESA population")
args <- parser$parse_args()

# set prefix name for model
out_prefix<- 'models/' %&% args$tis %&% '_' %&% args$pop %&% '_TIGAR_reduced_noh2_'

# getting h2 genes
heritable_genes <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% 
  filter(h2-2*se > 0.01, tissue==args$tis) %>% select(gene) %>% unique() %>% pull()

# read in weights for all chromossomes
for (i in c(1:23, 25)){
  tmp <- fread('outputs/'%&% args$tis %&%'/'%&% args$pop %&%'/DPR_CHR'%&% i %&%'/CHR'%&% i %&%'_DPR_train_eQTLweights.txt.gz') %>% 
    select(TargetID, POS, snpID, REF, ALT, beta)
  tmp$POS <- 'chr' %&% i %&%':'%&% tmp$POS
  tmp$snpID <- 'chr' %&% tmp$snpID
  colnames(tmp) <- c('gene', 'rsid', 'varID', 'ref_allele', 'eff_allele', 'weight')
  
  if (exists('final_weights_table')){
    final_weights_table <- rbind(final_weights_table, tmp)
  } else { final_weights_table <- tmp }
  
  tmp <- fread('outputs/'%&% args$tis %&%'/'%&% args$pop %&%'/DPR_CHR'%&% i %&%'/CHR'%&% i %&%'_DPR_train_GeneInfo.txt') %>% select(TargetID, GeneName)
 
  if (exists('final_extra_table')){
    final_extra_table <- rbind(final_extra_table, tmp)
  } else { final_extra_table <- tmp }
}

# filtering gene list
final_weights_table <- final_weights_table %>% filter(gene %in% heritable_genes) 
final_extra_table <- final_extra_table %>% filter(TargetID %in% heritable_genes) 

# reducing tables
for (j in c(1e-4)){
  red_final_weights_table <- final_weights_table %>% filter(weight >= j) 
  snps_p_gene <- red_final_weights_table %>% group_by(gene) %>% summarise(n.snps.in.model=n())
  red_final_extra_table <- inner_join(final_extra_table, snps_p_gene, by=c('TargetID'='gene')) %>% mutate(pred.perf.R2=NA, pred.perf.pval=NA, pred.perf.qval=NA)
  colnames(red_final_extra_table) <- c('gene', 'genename', 'n.snps.in.model', 'pred.perf.R2', 'pred.perf.pval', 'pred.perf.qval')

  # final files
  fwrite(red_final_extra_table, out_prefix %&% j %&%'_summaries.txt', col.names=T, quote=F, sep=' ')
  fwrite(red_final_weights_table, out_prefix %&% j %&%'_weights.txt', col.names=T, quote=F, sep=' ')
  conn <- dbConnect(drv = driver, out_prefix %&% j %&%'.db')
  dbWriteTable(conn, 'extra', red_final_extra_table, overwrite = TRUE)
  dbExecute(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbWriteTable(conn, 'weights', red_final_weights_table, overwrite = TRUE)
  dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbDisconnect(conn)
}