# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/04pQTL

#devtools::install_github("argparse", "trevorld")
library(argparse)
suppressPackageStartupMessages(library(MatrixEQTL))
parser <- ArgumentParser()
parser$add_argument("-sg", "--snpgenotype", help="file path of the snp genotype file")
parser$add_argument("-sl", "--snplocation", help="file path of the snp locaiton file")
parser$add_argument("-ge", "--geneexpression", help="file path of the gene expression file")
parser$add_argument("-gl", "--genelocation", help="file path of the gene location file")
parser$add_argument("-t", "--tag", help="file tag for this run of samples")
parser$add_argument("-o", "--outputdir", help="file tag for this run of samples", type="character", default="./")
parser$add_argument("--cov", help="file path to covariates file for MEQTL")
parser$add_argument("--cis", help="threshold for writing cis snp significance", type="double", default=1)
parser$add_argument("--trans", help="threshold for writing trans snp significance", type="double", default=0)
parser$add_argument("--window", help="maximum distance between snps to be considered cis", type="double", default=1e6)

args <- parser$parse_args()

#strip = gregexpr(pattern = 'chr', text = args$genelocation, ignore.case = T)
CHRnum <- args$tag  #substr(args$genelocation, strip, nchar(args$genelocation))

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = args$snpgenotype;
snps_location_file_name = args$snplocation;

# Gene expression file name
expression_file_name = args$geneexpression;
gene_location_file_name = args$genelocation;

# Covariates file name
# Set to character() for no covariates
if (!is.null(args$cov)){
  covariates_file_name = args$cov;
} else {
  covariates_file_name = character();
}


# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = paste(args$outputdir, "/trans_eQTLs_", CHRnum, ".txt", sep = "");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = args$cis;
pvOutputThreshold_tra = args$trans;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = args$window;
str(cisDist)
## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = ""; # the TAB character
snps$fileOmitCharacters = "-1"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
print("gene")
## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = ""; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "-1"; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
if (!is.null(args$cov)){
  cvrt$LoadFile(covariates_file_name);
} 


## Run the analysis
print("SNP")
snpspos = read.table(snps_location_file_name, sep="",header = TRUE, stringsAsFactors = FALSE);
str(snpspos)
print("gene")
genepos = read.csv(gene_location_file_name, header = TRUE, sep='', stringsAsFactors = FALSE);
str(genepos)
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
unlink(output_file_name_tra);
unlink(output_file_name_cis);
library(data.table)

## Results:
CisOutput<- paste(args$outputdir, "/cis_eQTLs_", CHRnum, ".txt", sep = "")
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
fwrite(me$cis$eqtls, CisOutput, sep = '\t')


#TransOutput <- paste(args$outputdir, "/trans_eQTLs_", CHRnum, ".txt", sep = "") #must hardcode output path for this to work
#cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected QTLs:', '\n');
#str(me$all$eqtls)
#fwrite(me$all$eqtls, TransOutput, sep = '\t')

## Malske the histogram of local and distant p-values
PlotOutput = paste(args$outputdir, "/cis_hist_", CHRnum, ".pdf", sep = "") #must hardcode path here as well
pdf(PlotOutput)
plot(me)
dev.off()
