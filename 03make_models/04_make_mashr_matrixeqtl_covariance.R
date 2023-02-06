# same script for both MASHR and MatrixeQTL models

library(data.table)
library(tidyverse)
library(reshape2)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep = "")

parser <- ArgumentParser()
parser$add_argument("--chr", help="chromosome whose files will be analyzed")
parser$add_argument("--pop", help="population whose files will be analyzed")
parser$add_argument("--tissue", help="tissue whose files will be analyzed")
args <- parser$parse_args()

get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- snp_annot_file_name %>%
    filter(!((refAllele == 'A' & effectAllele == 'T') |
               (refAllele == 'T' & effectAllele == 'A') |
               (refAllele == 'C' & effectAllele == 'G') |
               (refAllele == 'G' & effectAllele == 'C')) &
             !(is.na(rsid))) %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(rsid)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  cis_gt <- gt_df %>% select(one_of(intersect(snp_info$rsid, colnames(gt_df)))) %>% unique()
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language.
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}

get_genotype <- function(genotype_file_name) {
  n <- genotype_file_name$snp_ID
  gt_df <- genotype_file_name[,-1] %>% t() %>% as.data.frame()
  colnames(gt_df) <- n
  gt_df
}

do_covariance <- function(gene_id, cis_gt, rsids) {
  model_gt <- cis_gt[,rsids, drop=FALSE]
  colnames(model_gt) <- rsids
  geno_cov <- cov(model_gt)
  geno_cov[lower.tri(geno_cov)] <- NA
  cov_df <- melt(geno_cov, varnames = c("rsid1", "rsid2"), na.rm = TRUE) %>%
    mutate(gene=gene_id) %>%
    select(GENE=gene, RSID1=rsid1, RSID2=rsid2, VALUE=value) %>%
    arrange(GENE, RSID1, RSID2)
  cov_df
}

get_gene_annotation <- function(gene_annot_file_name, chrom){
  gene_df <- gene_annot_file_name %>% filter((chr == chrom))
  gene_df
}


####

snp_annot_file <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/snp_annotation/MetaPop.' %&% args$tissue %&% '.WG_noDup_locations.txt.gz', header=T, stringsAsFactors=F) %>% get_filtered_snp_annot()
gene_annot_file <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/gene_annotation_v34_60669.txt.gz', header=T, stringsAsFactors=F)
samples <- fread('/home/daniel/MESA_genotypes_subset/WGS_files/dosages/MESA_TOPMed_WGS_' %&% args$tissue %&% '_' %&% args$pop %&% '_unfilt.samples.txt', header=F, stringsAsFactors=F) %>% select(V1)
summary <- fread('/home/daniel/MESA_genotypes_subset/IDs_files/' %&% args$tissue %&% '_ALL_fulldata.txt',header=T) %>% right_join(samples, by=c('wgs_id'='V1')) %>% pull(sidno)
genotype_file <- fread('/home/daniel/MESA_genotypes_subset/WGS_files/dosages/MESA_TOPMed_WGS_' %&% args$tissue %&% '_' %&% args$pop %&% '_unfilt.chr' %&% args$chr %&% '.dosage.txt.gz', header=F, stringsAsFactors=F)
genotype_file <- setNames(genotype_file, c('chr', 'varID', 'pos', 'ref', 'alt', 'altfreq', as.character(summary)))
genotype_file <- mutate(genotype_file, snp_ID=paste(chr,':',pos, sep='')) %>% select(-c(chr, varID, pos, ref, alt, altfreq)) %>% select(snp_ID, everything())
unique_snps <- genotype_file %>% select(snp_ID)
unique_snps <- unique_snps[!(duplicated(unique_snps) | duplicated(unique_snps, fromLast=T)),]
genotype_file <- unique_snps %>% left_join(genotype_file)

prefix<-'/home/daniel/mashr/mashr_db/WGS_files/covariance_files/' %&% args$tissue %&% '_' %&% args$pop %&% '_mashr_baseline'
covariance_file <- prefix %&% '_chr' %&% args$chr %&% '_covariances.txt'
covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
write(covariance_col, file=covariance_file, ncol=4, sep=' ')

gene_annot <- get_gene_annotation(gene_annot_file, args$chr)
genes_in_chr <- gene_annot %>% pull(gene_id)
genes <- fread('/home/daniel/mashr/mashr_db/WGS_files/mashr_models/' %&% args$tissue %&% '_' %&% args$pop %&% '_mashr_baseline_summaries.txt.gz', fill=T) %>% pull(gene)
genes <- intersect(genes, genes_in_chr)

h2estimates <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% filter(h2-2*se > 0.01, tissue==args$tissue) %>% select(gene) %>% unique()
genes <- intersect(genes, h2estimates$gene)

for (i in 1:length(genes)){
  working_gene <- genes[i]
  print(working_gene)
  coords <- get_gene_coords(gene_annot, working_gene)
  gt_df <- get_genotype(genotype_file)
  cis_gt <- get_cis_genotype(gt_df, snp_annot_file, coords, 1000000)
  weighted_snps_info <- fread('/home/daniel/mashr/mashr_db/WGS_files/mashr_models/' %&% args$tissue %&% '_' %&% args$pop %&% '_mashr_baseline_weights.txt.gz', header=T, stringsAsFactors=F) %>% filter(gene==working_gene) 
  covariance_df <- do_covariance(working_gene, cis_gt, weighted_snps_info$rsid) %>% unique()
  write.table(covariance_df, file=covariance_file, append=T, quote=F, col.names=F, row.names=F, sep=' ')
}