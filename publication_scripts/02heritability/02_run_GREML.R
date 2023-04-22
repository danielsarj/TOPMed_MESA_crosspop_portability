library(data.table)
library(tidyverse)
library(argparse)
'%&%' = function(a,b) paste (a,b,sep='')

parser <- ArgumentParser()
parser$add_argument("--pop", help="population whose files will be analyzed")
parser$add_argument("--tissue", help="tissue whose files will be analyzed")
parser$add_argument("--r2", help="LD r2 to use")
args <- parser$parse_args()

# reading IDs file
ids_file <- fread('/home/daniel/MESA_genotypes_subset/IDs_files/'%&% args$tissue %&%'_ALL_fulldata.txt') %>% filter(pop==args$pop) %>% select(sidno, wgs_id) %>% unique()
id_save <- cbind(ids_file$wgs_id, ids_file$wgs_id) %>% as.data.frame()
fwrite(id_save, '/home/daniel/MESA_heritability/plink_files/'%&% args$tissue %&%'_'%&% args$pop %&%'.txt', col.names=F, quote=F, sep=' ')
  
# reading gene list 
gene_list <- fread('/home/daniel/MESA_heritability/gene_annotation_v38_strandTSS.txt.gz') %>% select(chr, gene_id, start, end) %>% filter(chr!='M', chr!='Y')
gene_list$chr <- as.numeric(gene_list$chr)
gene_list$chr[gene_list$chr==25] <- 23

# reading expression tbl
gene_exp_tbl <- fread('/home/daniel/MESA_heritability/'%&% args$tissue %&%'_expression_'%&% args$pop %&%'age_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10.txt.gz') %>%
  left_join(ids_file) %>% right_join(id_save, by=c('wgs_id'='V1'))
genes_to_loop <- grep('ENSG', colnames(gene_exp_tbl), value=T) %>% intersect(gene_list$gene_id)

for (working_gene in genes_to_loop){
  print(args$tissue %&%' -- '%&% args$pop %&% ' -- ' %&% working_gene %&% ' -- start')
  chr_tss <- gene_list %>% filter(gene_id==working_gene) %>% mutate(left_tss=start-1000000, right_tss=end+1000000)
  if (chr_tss$left_tss <= 0){ chr_tss$left_tss <- 1 }
  chr_tss <- c(chr_tss$chr, chr_tss$left_tss, chr_tss$right_tss)
  
  # making plink bfiles w/ SNPs in 1MB window 
  system('plink --bfile /home/daniel/MESA_heritability/plink_files/MESA_TOPMed_WGSX.'%&% args$tissue %&%'.'%&% args$pop %&%'.'%&% as.character(args$r2) %&%'.hg38.pruned --chr '%&% chr_tss[1] %&%' --from-bp '%&% chr_tss[2] %&%' --to-bp '%&% chr_tss[3] %&%' --maf 0.01 --make-bed --out /home/daniel/MESA_heritability/plink_files/per_pop_per_gene/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%'')

  # making GRM
  if (chr_tss[1] != 23){
    system('gcta64 --bfile /home/daniel/MESA_heritability/plink_files/per_pop_per_gene/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%' --make-grm --out /home/daniel/MESA_heritability/GRMs/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%' --thread-num 10')
  } else {
    system('gcta64 --bfile /home/daniel/MESA_heritability/plink_files/per_pop_per_gene/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%' --make-grm-xchr --out /home/daniel/MESA_heritability/GRMs/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%' --thread-num 10')
  }

  # making pheno file
  pheno_file <- gene_exp_tbl %>% select(wgs_id, working_gene) %>% mutate(wgs_id_d = wgs_id) %>% select(wgs_id, wgs_id_d, working_gene)
  fwrite(pheno_file, '/home/daniel/MESA_heritability/pheno_files/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'.txt', col.names=F, quote=F, sep=' ')

  # estimating h2
  system('gcta64 --reml --grm /home/daniel/MESA_heritability/GRMs/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%' --pheno /home/daniel/MESA_heritability/pheno_files/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'.txt --out /home/daniel/MESA_heritability/greml_outputs/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_r'%&% as.character(args$r2) %&%'_constrained')
  system('gcta64 --reml --reml-no-constrain --grm /home/daniel/MESA_heritability/GRMs/MESA_TOPMed_'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_1Mbwindow_r'%&% as.character(args$r2) %&%' --pheno /home/daniel/MESA_heritability/pheno_files/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'.txt --out /home/daniel/MESA_heritability/greml_outputs/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_r'%&% as.character(args$r2) %&%'_noconstrained')
  
  # reading result
  if (file.exists('/home/daniel/MESA_heritability/greml_outputs/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_r'%&% as.character(args$r2) %&%'_noconstrained.hsq')==T){
    hsq <- scan('/home/daniel/MESA_heritability/greml_outputs/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_r'%&% as.character(args$r2) %&%'_noconstrained.hsq','character')
    res <- c(working_gene, hsq[14], hsq[15], hsq[25]) %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(res) <- c('gene', 'h2', 'se', 'pval')
  
    if (exists('final.hsq.df_nocon')){
      final.hsq.df_nocon <- rbind(final.hsq.df_nocon, res)
    } else {final.hsq.df_nocon <- res}
    print(args$tissue %&%' -- '%&% args$pop %&% ' -- ' %&% working_gene %&% ' -- successfully finished')
  } else {
    res <- c(working_gene) %>% as.data.frame() 
    colnames(res) <- c('gene')
    
    if (exists('final.hsq.df.wrong_nocon')){
      final.hsq.df.wrong_nocon <- rbind(final.hsq.df.wrong_nocon, res)
    } else {final.hsq.df.wrong_nocon <- res}
    print(args$tissue %&%' -- '%&% args$pop %&% ' -- ' %&% working_gene %&% ' -- went wrong')
  }
  
  if (file.exists('/home/daniel/MESA_heritability/greml_outputs/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_r'%&% as.character(args$r2) %&%'_constrained.hsq')==T){
    hsq <- scan('/home/daniel/MESA_heritability/greml_outputs/'%&% args$tissue %&%'_'%&% args$pop %&%'_'%&% working_gene %&%'_r'%&% as.character(args$r2) %&%'_constrained.hsq','character')
    res <- c(working_gene, hsq[14], hsq[15], hsq[25]) %>% as.data.frame() %>% t() %>% as.data.frame()
    colnames(res) <- c('gene', 'h2', 'se', 'pval')
    
    if (exists('final.hsq.df_con')){
      final.hsq.df_con <- rbind(final.hsq.df_con, res)
    } else {final.hsq.df_con <- res}
    print(args$tissue %&%' -- '%&% args$pop %&% ' -- ' %&% working_gene %&% ' -- successfully finished')
  } else {
    res <- c(working_gene) %>% as.data.frame() 
    colnames(res) <- c('gene')
    
    if (exists('final.hsq.df.wrong_con')){
      final.hsq.df.wrong_con <- rbind(final.hsq.df.wrong_con, res)
    } else {final.hsq.df.wrong_con <- res}
    print(args$tissue %&%' -- '%&% args$pop %&% ' -- ' %&% working_gene %&% ' -- went wrong')
  }
}

if (exists('final.hsq.df_con')){
  fwrite(final.hsq.df_con, '/home/daniel/MESA_heritability/'%&% args$tissue %&%'_'%&% args$pop %&%'_hsqestimates_r'%&% as.character(args$r2) %&%'_constrained.txt', col.names=T, sep=' ', quote=F)
}
if (exists('final.hsq.df.wrong_con')){
  fwrite(final.hsq.df.wrong_con, '/home/daniel/MESA_heritability/'%&% args$tissue %&%'_'%&% args$pop %&%'_genes_that_went_wrong_r'%&% as.character(args$r2) %&%'_constrained.txt', col.names=T, sep=' ', quote=F)
}
if (exists('final.hsq.df_nocon')){
  fwrite(final.hsq.df_nocon, '/home/daniel/MESA_heritability/'%&% args$tissue %&%'_'%&% args$pop %&%'_hsqestimates_r'%&% as.character(args$r2) %&%'_noconstrained.txt', col.names=T, sep=' ', quote=F)
}
if (exists('final.hsq.df.wrong_nocon')){
  fwrite(final.hsq.df.wrong_nocon, '/home/daniel/MESA_heritability/'%&% args$tissue %&%'_'%&% args$pop %&%'_genes_that_went_wrong_r'%&% as.character(args$r2) %&%'_noconstrained.txt', col.names=T, sep=' ', quote=F)
}