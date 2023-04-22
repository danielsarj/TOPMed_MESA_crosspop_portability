library(tidyverse)
library(RSQLite)
library(qvalue)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep='')
driver <- dbDriver('SQLite')
setwd('/home/daniel/mashr/EN/output_files')

for (tis in c('Mono','PBMC','Tcell')){
  for (pop in c('AFA','EUR','HIS','CHN')){
    if (tis!='PBMC' & pop=='CHN'){
      next
    } else {
      
      h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01, tissue==tis) %>% select(gene) %>% unique()
  
      cat(tis,"\t",pop,"\n")
      out_prefix<-'/home/daniel/mashr/final_models/' %&% tis %&% '_' %&% pop %&% '_EN_baseline'
  
      model_summaries <- read.table(tis %&% '_' %&% pop %&% '_chr1_base_chr1_model_summaries.txt', header = T, stringsAsFactors = F)
      tiss_summary <- read.table(tis %&% '_' %&% pop %&% '_chr1_base_chr1_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
      weights <- read.table(tis %&% '_' %&% pop %&% '_chr1_base_chr1_weights.txt', header = T, stringsAsFactors = F)
      covariances<-read.table(tis %&% '_' %&% pop %&% '_chr1_base_chr1_covariances.txt', header = T, stringsAsFactors = F)
      n_samples <- tiss_summary$n_samples
  
      cat("Reading in results per chr \n")
      for (chr in c(c(2:23),25)){
        cat(chr,"\n")
        tmp <- read.table(tis %&% '_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_model_summaries.txt', header = T, stringsAsFactors = F)
        model_summaries<-rbind.data.frame(model_summaries,tmp)
        tmp2 <- read.table(tis %&% '_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F)
        tiss_summary<-rbind.data.frame(tiss_summary,tmp2)
        tmp3 <- read.table(tis %&% '_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_weights.txt', header = T, stringsAsFactors = F)
        weights<-rbind.data.frame(weights,tmp3)
        tmp4 <-read.table(tis %&% '_' %&% pop %&% '_chr' %&% chr %&% '_base_chr' %&% chr %&% '_covariances.txt', header = T, stringsAsFactors = F)
        covariances<-rbind.data.frame(covariances,tmp4)
      }
  
  
      weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
      weights <- weights %>% filter(gene %in% h2estimates$gene)
      covariances <- covariances %>% filter(GENE %in% h2estimates$gene)
      model_summaries <- model_summaries %>% filter(gene_id %in% h2estimates$gene) %>% filter(gene_id %in% weights$gene)
      sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = tis)
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
                            )
      
      if (length(qvalues) == 0){
        model_summaries <- model_summaries %>% mutate(pred.perf.qval = 0)
      } else {
        model_summaries <- model_summaries %>% mutate(pred.perf.qval = qvalues$qvalues)
      }
  
      conn <- dbConnect(drv = driver, out_prefix %&% '.db')
      dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
      dbGetQuery(conn, "CREATE INDEX gene_model_summary ON extra (gene)")
      dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
      dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
      dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
      dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
      dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
      dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
      fwrite(covariances, out_prefix %&% '_covariances.txt',col.names = T,sep=' ',row.names=F,quote = F)
      fwrite(weights, out_prefix %&% '_weights.txt',col.names = T,sep=' ',row.names=F,quote = F)
      fwrite(model_summaries, out_prefix %&% '_summaries.txt',col.names = T,sep=' ',row.names=F,quote = F)
    }
  }
}