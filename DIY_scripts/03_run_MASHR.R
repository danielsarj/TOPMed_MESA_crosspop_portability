# Loading libraries and defining arguments
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))
suppressMessages(library(mashr))
'%&%' = function(a,b) paste (a,b,sep='')
parser <- ArgumentParser()
parser$add_argument('-i', '--input', help='path of the directory with input files')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-o', '--output', help='path of the output directory')
args <- parser$parse_args()

# Change working directory to where input files are 
setwd(args$input)

# Get list of genes
gene_list <- fread(args$geneannotation) %>% pull(gene_id) %>% unique()  

# Run MASHR for each gene at a time
for (working_gene in gene_list){
  print('INFO: Running MASHR with gene ' %&% working_gene)
  
  # Load beta and SE dfs as matrices
  beta <- fread(working_gene %&% '_beta.txt.gz', header=T, stringsAsFactors=F) %>% select(contains('beta')) %>% as.matrix()
  se <- fread(working_gene %&% '_SE.txt.gz', header=T, stringsAsFactors=F) %>% select(contains('SE')) %>% as.matrix() %>% abs() 
    #use abs() as recommended by mashr:
    #'Both Bhat and Shat are zero (or near zero) for some input data. Please check your input. 
    #If it is expected please set Shat to a positive number to avoid numerical issues;'
  df_anno <- fread(working_gene %&% '_beta.txt.gz', header=T, stringsAsFactors=F) %>% select(gene, snps, snp_ID)
  
  if (nrow(df_anno)==1){
    print('INFO: Not enough SNPs in the input data frame to run MASHR for ' %&% working_gene %&%'. Skipping gene')
    next
  } else {
    # Set up main MASHR data object
    data = mash_set_data(beta, se)
  
    # Get covariance matrices
    data.c = cov_canonical(data) # canonical
    data.pca = cov_pca(data, min(ncol(beta),nrow(beta))) # pca
    data.ed = cov_ed(data, data.pca) # data-driven
  
    # Fit model
    print('INFO: Fitting model')
    m = mash(data, Ulist=c(data.ed,data.c))
  
    # Get posterior summaries
    posterior_lfsr <- get_lfsr(m) # local false sign rate
    posterior_lfsr <- cbind(df_anno, posterior_lfsr)
    colnames(posterior_lfsr) <- gsub('_beta', '_lfsr', colnames(posterior_lfsr))
    posterior_mean <- get_pm(m) # new betas
    posterior_mean <- cbind(df_anno, posterior_mean)
    posterior_sd <- get_psd(m) # standard deviantion
    posterior_sd <- cbind(df_anno, posterior_sd)
    colnames(posterior_sd) <- gsub('_beta', '_SD', colnames(posterior_sd))

    # Write output
    fwrite(posterior_lfsr, file=args$output %&% '/' %&% working_gene %&% '_MASHR_lfsr.txt', quote=F, sep=' ')
    fwrite(posterior_mean, file=args$output %&% '/' %&% working_gene %&% '_MASHR_beta.txt', quote=F, sep=' ')
    fwrite(posterior_sd, file=args$output %&% '/' %&% working_gene %&% '_MASHR_SD.txt', quote=F, sep=' ')
    print('INFO: Successfully ran MASHR for ' %&% working_gene)
  }
}