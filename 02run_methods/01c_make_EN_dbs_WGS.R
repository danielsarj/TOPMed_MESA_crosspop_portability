# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/06Elastic_net

library(dplyr)
library(RSQLite)
library(qvalue)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
populations<-c("ALL")

for (pop in populations){
  out_prefix<-'/home/daraujo1/data/TOPMed_MESA/EN_outputs/EN_WGS_ALL/PBMC_' %&% pop %&% '_PCAIR_baseline_models'

  model_summaries <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_base_chr1_model_summaries.txt', header = T, stringsAsFactors = F)
  tiss_summary <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_base_chr1_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
  weights <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_base_chr1_weights.txt', header = T, stringsAsFactors = F)
  covariances<-read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_base_chr1_covariances.txt', header = T, stringsAsFactors = F)
  n_samples <- tiss_summary$n_samples

  cat("Reading in results per chr \n")
  for ( chr in c(2:23, 25)){
    cat(chr,"\n")
    tmp <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_model_summaries.txt', header = T, stringsAsFactors = F)
    model_summaries<-rbind.data.frame(model_summaries,tmp)
    tmp2 <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
    tiss_summary<-rbind.data.frame(tiss_summary,tmp2)
    tmp3 <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_weights.txt', header = T, stringsAsFactors = F)
    weights<-rbind.data.frame(weights,tmp3)
    tmp4<-read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_covariances.txt', header = T, stringsAsFactors = F)
    covariances<-rbind.data.frame(covariances,tmp4)
  }

  weights <- rename(weights,gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
  sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = "PBMC")
  construction <- tiss_summary %>% select(cv_seed)
  pvalues<-model_summaries$zscore_pval
  qvalues<-tryCatch(qvalue(pvalues), error = function(cond) {message('Error'); message(geterrmessage()); list()})
  model_summaries <- rename(model_summaries,
                            gene = gene_id,
                            genename = gene_name,
                            n.snps.in.window = n_snps_in_window,
                            n.snps.in.model = n_snps_in_model,
                            pred.perf.R2 = rho_avg_squared,
                            pred.perf.pval = zscore_pval
                            )# %>% mutate(pred.perf.qval = qvalues$qvalues)
  if (length(qvalues) == 0){
    model_summaries <- model_summaries %>% mutate(pred.perf.qval = 0)
  } else {
    model_summaries <- model_summaries %>% mutate(pred.perf.qval = qvalues$qvalues)
  }

  conn <- dbConnect(drv = driver, out_prefix %&% '_unfiltered.db')
  dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  dbWriteTable(conn, 'construction', construction, overwrite = TRUE)

  fwrite(covariances,out_prefix %&% '_unfiltered_covariances.txt',col.names = T,sep='\t',row.names=F,quote = F)

  ### Write only significant models
  conn3 <- dbConnect(drv = driver, out_prefix %&%  '_rho0.1_zpval0.05.db')
  signif_models<-filter(model_summaries, pred.perf.pval < 0.05 , rho_avg > 0.1)
  dbWriteTable(conn3, 'extra', signif_models, overwrite = TRUE)
  dbGetQuery(conn3, "CREATE INDEX gene_model_summary ON extra (gene)")
  weights_signif<- weights %>% filter ( gene %in% signif_models$gene)
  dbWriteTable(conn3, 'weights', weights_signif, overwrite = TRUE)
  dbGetQuery(conn3, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(conn3, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(conn3, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbWriteTable(conn3, 'sample_info', sample_info, overwrite = TRUE)
  dbWriteTable(conn3, 'construction', construction, overwrite = TRUE)

  covariances_signif<-covariances %>% filter(GENE %in% signif_models$gene)

  fwrite(covariances_signif,out_prefix %&% '_rho0.1_zpval0.05_covariances.txt',col.names = T,sep='\t',row.names=F,quote = F)
}

for (pop in populations){
  out_prefix<-'/home/daraujo1/data/TOPMed_MESA/EN_outputs/EN_WGS_ALL/PBMC_' %&% pop %&% '_PAGE_PCAIR_baseline_models'

  model_summaries <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_PAGE_base_chr1_model_summaries.txt', header = T, stringsAsFactors = F)
  tiss_summary <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_PAGE_base_chr1_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
  weights <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_PAGE_base_chr1_weights.txt', header = T, stringsAsFactors = F)
  covariances<-read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr1_PAGE_base_chr1_covariances.txt', header = T, stringsAsFactors = F)
  n_samples <- tiss_summary$n_samples

  cat("Reading in results per chr \n")
  for ( chr in c(2:23, 25)){
    cat(chr,"\n")
    tmp <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_PAGE_base_chr' %&% chr %&% '_model_summaries.txt', header = T, stringsAsFactors = F)
    model_summaries<-rbind.data.frame(model_summaries,tmp)
    tmp2 <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_PAGE_base_chr' %&% chr %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
    tiss_summary<-rbind.data.frame(tiss_summary,tmp2)
    tmp3 <- read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_PAGE_base_chr' %&% chr %&% '_weights.txt', header = T, stringsAsFactors = F)
    weights<-rbind.data.frame(weights,tmp3)
    tmp4<-read.table('/home/daraujo1/data/TOPMed_MESA/EN_inputs/EN_WGS_ALL/output/PBMC_' %&% pop %&% '_chr' %&% chr %&% '_PAGE_base_chr' %&% chr %&% '_covariances.txt', header = T, stringsAsFactors = F)
    covariances<-rbind.data.frame(covariances,tmp4)
  }

  weights <- rename(weights,gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
  sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = "PBMC")
  construction <- tiss_summary %>% select(cv_seed)
  pvalues<-model_summaries$zscore_pval
  qvalues<-tryCatch(qvalue(pvalues), error = function(cond) {message('Error'); message(geterrmessage()); list()})
  model_summaries <- rename(model_summaries,
                            gene = gene_id,
                            genename = gene_name,
                            n.snps.in.window = n_snps_in_window,
                            n.snps.in.model = n_snps_in_model,
                            pred.perf.R2 = rho_avg_squared,
                            pred.perf.pval = zscore_pval
                            )# %>% mutate(pred.perf.qval = qvalues$qvalues)
  if (length(qvalues) == 0){
    model_summaries <- model_summaries %>% mutate(pred.perf.qval = 0)
  } else {
    model_summaries <- model_summaries %>% mutate(pred.perf.qval = qvalues$qvalues)
  }

  conn <- dbConnect(drv = driver, out_prefix %&% '_unfiltered.db')
  dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
  dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
  dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
  dbWriteTable(conn, 'construction', construction, overwrite = TRUE)

  fwrite(covariances,out_prefix %&% '_unfiltered_covariances.txt',col.names = T,sep='\t',row.names=F,quote = F)

  ### Write only significant models
  conn3 <- dbConnect(drv = driver, out_prefix %&%  '_rho0.1_zpval0.05.db')
  signif_models<-filter(model_summaries, pred.perf.pval < 0.05 , rho_avg > 0.1)
  dbWriteTable(conn3, 'extra', signif_models, overwrite = TRUE)
  dbGetQuery(conn3, "CREATE INDEX gene_model_summary ON extra (gene)")
  weights_signif<- weights %>% filter ( gene %in% signif_models$gene)
  dbWriteTable(conn3, 'weights', weights_signif, overwrite = TRUE)
  dbGetQuery(conn3, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(conn3, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(conn3, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbWriteTable(conn3, 'sample_info', sample_info, overwrite = TRUE)
  dbWriteTable(conn3, 'construction', construction, overwrite = TRUE)

  covariances_signif<-covariances %>% filter(GENE %in% signif_models$gene)

  fwrite(covariances_signif,out_prefix %&% '_rho0.1_zpval0.05_covariances.txt',col.names = T,sep='\t',row.names=F,quote = F)
}
