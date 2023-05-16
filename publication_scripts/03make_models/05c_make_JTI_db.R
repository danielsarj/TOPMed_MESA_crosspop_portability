library(tidyverse)
library(data.table)
library(argparse)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
setwd('/home/daniel/mashr/MR-JTI/TOPMed_MESA')

# arguments
parser <- ArgumentParser()
parser$add_argument("--tis", help="MESA tissue")
parser$add_argument("--pop", help="MESA population")
args <- parser$parse_args()

# set prefix name for model
out_prefix<- 'models/' %&% args$tis %&% '_' %&% args$pop %&% '_JTI_baseline'

# getting h2 genes
heritable_genes <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% 
  filter(h2-2*se > 0.01, tissue==args$tis) %>% select(gene) %>% unique() %>% pull()

# getting list of files to parse through
list_of_files <- list.files(path = 'outputs/' %&% args$tis %&% '/' %&% args$pop, pattern = '.txt')

# read in weights for all genes
for (f in list_of_files){
  tmp <- fread('outputs/' %&% args$tis %&% '/' %&% args$pop %&% '/' %&% f) %>% select(gene, chr_bp, ref_allele, counted_allele, weight)
  tmp$chr_bp <- 'chr' %&% tmp$chr_bp %>% gsub('_', ':', .)
  tmp$varID <- tmp$chr_bp %&% ':' %&% tmp$ref_allele %&% ':' %&% tmp$counted_allele
  colnames(tmp) <- c('gene', 'rsid', 'ref_allele', 'eff_allele', 'weight', 'varID')
  tmp <- tmp %>% select(gene, rsid, varID, ref_allele, eff_allele, weight)
  
  if (exists('final_weights_table')){
    final_weights_table <- rbind(final_weights_table, tmp)
  } else { final_weights_table <- tmp }
}

# filtering gene list
final_weights_table <- final_weights_table %>% filter(gene %in% heritable_genes) 

# making extra table
genes_list <- final_weights_table %>% pull(gene) %>% unique()
genes_table <- table(final_weights_table$gene) %>% as.data.frame() %>% rename(gene=Var1, n_snps=Freq)
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

# making final files
fwrite(model_summaries, out_prefix %&% '_summaries.txt', col.names=T, quote=F, sep=' ')
fwrite(final_weights_table, out_prefix %&% '_weights.txt', col.names=T, quote=F, sep=' ')
conn <- dbConnect(drv = driver, out_prefix %&% '.db')
dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
dbWriteTable(conn, 'weights', final_weights_table, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbDisconnect(conn)