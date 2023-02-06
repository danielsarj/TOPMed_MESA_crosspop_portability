suppressMessages(library(tidyverse))
suppressMessages(library(RSQLite))
suppressMessages(library(data.table))
suppressMessages(library(readr))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
setwd('/home/daniel/mashr/mashr_dfs/WGS_files/outputs')

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type whose files will be analyzed')
parser$add_argument('--pop', help='population whose files will be analyzed')
args <- parser$parse_args()

h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01, tissue==args$tissue) %>% select(gene) %>% unique()

out_prefix<- '/home/daniel/mashr/final_models/' %&% args$tissue %&% '_' %&% args$pop %&% '_mashr_baseline'
SNPs_anno <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/MetaPop.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz')

### reads mashr output file
#reading input
mashr_in <- fread('/home/daniel/mashr/mashr_db/WGS_files/top_' %&% args$tissue %&% '_mashr_sigSNPs-betas_shared_genes.txt.gz', stringsAsFactors=F) %>% 
  filter(gene %in% h2estimates$gene) %>% select(gene, snps, varID, contains(args$pop)) %>% rename(beta=contains(args$pop)) %>%
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