# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/04pQTL

# Loading libraries and defining arguments and functions
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(argparse))
"%&%" <- function(a,b) paste(a,b, sep = "")
parser <- ArgumentParser()
parser$add_argument('-d', '--snpdosage', help='file path of the unfilterd SNP dosage file')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-t', '--tag', help='file tag for this run of samples')
parser$add_argument('-c', '--chromosome', help='chromosome to be analyzed')
parser$add_argument('-m', '--model', help='path to model weight file for which covariance will be computed')
parser$add_argument('-o', '--outputdir', help='output directory path', type='character', default='./')
parser$add_argument('-w','--window', help='maximum distance between snps to be considered cis', type='double', default=1e6)
args <- parser$parse_args()
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
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt))
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

### 
print('INFO: '%&% args$tag %&%' '%&% args$chromosome)
print('INFO: Loading files')
# Get genes to compute covariance for
gene_annot <- fread(args$geneannotation, header=T, stringsAsFactors=F) %>% filter(chr==args$chromosome)
genes_in_model <- fread(args$model, fill=T) %>% pull(gene) %>% unique()
list_genes <- intersect(gene_annot$gene_id, genes_in_model)

# Read dosage file (unfiltered!)
genotype_file <- fread(args$snpdosage, header=T, stringsAsFactors=F)

# Get SNP annotation
snp_annot_file <- genotype_file %>% select(chr, snp_ID, pos, ref_allele, alt_allele) %>% mutate(rsid = 'chr' %&% chr %&% ':' %&% pos) %>% select(rsid, chr, snp_ID, pos, ref_allele, alt_allele)
colnames(snp_annot_file) <- c('rsid', 'chr', 'varID', 'pos', 'refAllele', 'effectAllele')

# Remove extra columns from dosage file
genotype_file <- genotype_file %>% mutate(snp_ID = 'chr' %&% chr %&% ':' %&% pos) %>% select(-c(chr, pos, ref_allele, alt_allele))

# Make empty cov file
covariance_file <- args$outputdir %&% '/' %&% args$tag %&%'_chr'%&% args$chromosome %&% '_covariances.txt'
covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
write(covariance_col, file=covariance_file, ncol=4, sep=' ')

# Computing covariance
print('INFO: Computing SNP-SNP covariance')
for (i in 1:length(list_genes)){
  working_gene <- list_genes[i]
  print('INFO: gene ' %&% working_gene)
  coords <- get_gene_coords(gene_annot, working_gene)
  gt_df <- get_genotype(genotype_file)
  cis_gt <- get_cis_genotype(gt_df, snp_annot_file, coords, args$window)
  weighted_snps_info <- fread(args$model, header=T, stringsAsFactors=F) %>% filter(gene==working_gene) 
  covariance_df <- do_covariance(working_gene, cis_gt, weighted_snps_info$rsid) %>% unique()
  write.table(covariance_df, file=covariance_file, append=T, quote=F, col.names=F, row.names=F, sep=' ')
}