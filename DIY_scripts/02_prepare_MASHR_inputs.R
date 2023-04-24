# Loading libraries and defining arguments
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(argparse))
'%&%' = function(a,b) paste (a,b,sep='')
parser <- ArgumentParser()
parser$add_argument('-l', '--listinputs', help='path to file containing inputs to be read for each condition')
parser$add_argument('-g', '--geneannotation', help='file path of the gene annotation file')
parser$add_argument('-c', '--chromosome', help='chromosome number to be analyzed')
parser$add_argument('-o', '--output', help='path of the output directory')
args <- parser$parse_args()

# Initialize lists of conditions
pop_codes <- list()
inputs <- list()

# Read file with population code, MatrixeQTL cis-eQTLs effect sizes, and dosage (per row)
files_to_analyze <- fread(args$listinputs, header=F)

# Read and save data into our lists
print('INFO: Reading MatrixeQTL and dosage files')
for (i in 1:nrow(files_to_analyze)){
  print('Current code: ' %&% files_to_analyze$V1[i])
  pop_codes[[i]] <- files_to_analyze$V1[i] # get code
  cis_es <- fread(files_to_analyze$V2[i], header=T, fill=T) %>% select(snps, gene, beta, SE) # get MatrixeQTL results
  dosage <- fread(files_to_analyze$V3[i], header=T, fill=T) %>% mutate(snps = 'chr' %&% chr %&% ':' %&% pos) %>%
    select(snps, snp_ID) # read in SNP annotation
  inputs[[i]] <- inner_join(dosage, cis_es, by=c('snps'='snps')) %>% suppressWarnings() %>% arrange(gene) %>% unique() %>% drop_na()
  print('Successfully read ' %&% files_to_analyze$V1[i] %&% ' input files')
}
rm(cis_es, dosage) # free memory

# Get list of genes
gene_list <- fread(args$geneannotation) %>% filter(chr==args$chromosome) %>% pull(gene_id) %>% unique()  

# Take a gene from the list, gets betas and SEs for all pops, and writes the output files
for (working_gene in gene_list){
  print('INFO: Making ' %&% working_gene %&% ' input data frames')
  
  # Initialize empty data frames
  beta_df <- data.frame()
  se_df <- data.frame()
    
  # Get data for each condition
  for (i in 1:length(pop_codes)){
    beta_tmp_df <- inputs[[i]] %>% filter(gene == working_gene) %>% select(-SE) # get beta info
    colnames(beta_tmp_df)[4] <- c(pop_codes[[i]] %&% '_beta')
    se_tmp_df <- inputs[[i]] %>% filter(gene == working_gene) %>% select(-beta) # get SE info
    colnames(se_tmp_df)[4] <- c(pop_codes[[i]] %&% '_SE')
    
    # Check if it is the first condition being read
    if (nrow(beta_df)==0){
      # If it is, just add data to final dfs
      beta_df <- beta_tmp_df
      se_df <- se_tmp_df
    } else {
      # If it is not, join it to the df 
      beta_df <- full_join(beta_df, beta_tmp_df, by=join_by(snps, snp_ID, gene)) %>% suppressWarnings() %>% arrange(snps)
      se_df <- full_join(se_df, se_tmp_df, by=join_by(snps, snp_ID, gene)) %>% suppressWarnings() %>% arrange(snps)
    }
  }

  # Write final files
  fwrite(beta_df, args$output %&% '/' %&% working_gene %&% '_beta.txt', quote=F, sep=' ', na=0)
  fwrite(se_df, args$output %&% '/' %&% working_gene %&% '_SE.txt', quote=F, sep=' ', na=10)
  print('INFO: Successfully made MASHR input files for ' %&% working_gene)
}