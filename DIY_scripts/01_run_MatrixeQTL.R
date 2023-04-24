# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/04pQTL

# Loading libraries and defining arguments
suppressMessages(library(argparse))
suppressMessages(library(MatrixEQTL))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
"%&%" <- function(a,b) paste(a,b, sep = "")
parser <- ArgumentParser()
parser$add_argument('-d', '--snpdosage', help='file path of the snp dosage file')
parser$add_argument('-e', '--geneexpression', help='file path of the gene expression file')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-t', '--tag', help='file tag for this run of samples')
parser$add_argument('-o', '--outputdir', help='file tag for this run of samples', type='character', default='./')
parser$add_argument('-w','--window', help='maximum distance between snps to be considered cis', type='double', default=1e6)
args <- parser$parse_args()

# Load genotype data
print('INFO: Loading dosage data')
snpspos <- fread(args$snpdosage, header=T) %>% mutate(rsid = 'chr' %&% chr %&% ':' %&% pos) %>% 
  select(rsid, chr, pos) # retrieve snp position info
tmp_snps <- fread(args$snpdosage, header=T) %>% mutate(snp_ID = 'chr' %&% chr %&% ':' %&% pos) %>% 
  dplyr::select(-c(chr, pos, ref_allele, alt_allele)) %>% as.matrix(rownames=1) # get dosage in the correct format to load into SlicedData()
snps <- SlicedData$new()
snps$fileDelimiter <- '' # space character
snps$fileSliceSize <- 2000 # read file in slices of 2,000 rows
suppressWarnings(snps$CreateFromMatrix(tmp_snps))
rm(tmp_snps) # free memory

# Load gene expression data
print('INFO: Loading gene expression data')
gene <- SlicedData$new()
gene$fileDelimiter <- '' # space character
gene$fileSkipRows <- 1 # one row of column labels
gene$fileSkipColumns <- 1 # one column of row labels
gene$fileSliceSize <- 2000 # read file in slices of 2,000 rows
suppressWarnings(gene$LoadFile(args$geneexpression))
genepos <- fread(args$geneannotation, header=T) %>%
  select(gene_id, chr, start, end) # retrieve gene TSS/TES info

# Run the analysis
print('INFO: Assessing cis-eQTLs effect sizes')
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  pvOutputThreshold = 0,
  pvOutputThreshold.cis = 1,
  useModel = modelLINEAR,
  cvrt = SlicedData$new(),
  errorCovariance = numeric(), 
  verbose = F,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = args$window,
  pvalue.hist = F,
  min.pv.by.genesnp = F,
  noFDRsaveMemory = F)

# Get SE info 
cis_es_out <- me$cis$eqtls %>% mutate(SE = beta/abs(statistic))

# Save results
print('INFO: Saving results')
fwrite(cis_es_out, args$outputdir %&% args$tag %&% '.txt', col.names=T, sep=' ')