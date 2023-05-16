library(data.table)
library(tidyverse)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep='')
setwd('/home/daniel/mashr/MR-JTI/TOPMed_MESA')

# arguments
parser <- ArgumentParser()
parser$add_argument("--tis", help="MESA tissue")
parser$add_argument("--pop", help="MESA population")
args <- parser$parse_args()

# get list of heritable genes
heritable_genes <- fread('/home/daniel/MESA_heritability/plots/significant_h2estimates_noconstrained_r0.2.txt') %>% 
  filter(h2-2*se > 0.01, tissue==args$tis) %>% select(gene) %>% unique() %>% pull()
  
# load correlation matrix
corr_matrix <- fread('inputs/'%&% args$tis %&%'_transcriptome_wide_median_correlation.txt')
  
# load expression tables
afa_exp <- fread('/home/daniel/TOPMed_MESA_RNAseq/exp_tables/' %&% args$tis %&% '_expression_AFAage_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10_transposed.txt', header=T)
eur_exp <- fread('/home/daniel/TOPMed_MESA_RNAseq/exp_tables/' %&% args$tis %&% '_expression_EURage_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10_transposed.txt', header=T)
his_exp <- fread('/home/daniel/TOPMed_MESA_RNAseq/exp_tables/' %&% args$tis %&% '_expression_HISage_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10_transposed.txt', header=T)
if (args$tis == 'PBMC'){
    chn_exp <- fread('/home/daniel/TOPMed_MESA_RNAseq/exp_tables/PBMC_expression_CHNage_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10_transposed.txt', header=T)
    genes_to_loop <- dplyr::intersect(afa_exp$gene_id, eur_exp$gene_id) %>% dplyr::intersect(his_exp$gene_id) %>% dplyr::intersect(chn_exp$gene_id)
  } else {
    genes_to_loop <- dplyr::intersect(afa_exp$gene_id, eur_exp$gene_id) %>% dplyr::intersect(his_exp$gene_id)
}

# looping through every gene
for (g in genes_to_loop){
  print(args$tis %&% ' - ' %&% args$pop %&% ' - ' %&% g %&% ' - start')

  # get respective pairwise PearsOn correlation values
  if (args$tis == 'PBMC' & args$pop == 'AFA'){
      afa_corr = corr_matrix[1, 2] %>% pull()
      eur_corr = corr_matrix[2, 2] %>% pull()
      his_corr = corr_matrix[3, 2] %>% pull()
      chn_corr = corr_matrix[4, 2] %>% pull()
    } else if (args$tis != 'PBMC' & args$pop == 'AFA'){
      afa_corr = corr_matrix[1, 2] %>% pull()
      eur_corr = corr_matrix[2, 2] %>% pull()
      his_corr = corr_matrix[3, 2] %>% pull()
    } else if (args$tis == 'PBMC' & args$pop == 'EUR'){
      afa_corr = corr_matrix[1, 3] %>% pull()
      eur_corr = corr_matrix[2, 3] %>% pull()
      his_corr = corr_matrix[3, 3] %>% pull()
      chn_corr = corr_matrix[4, 3] %>% pull()
    } else if (args$tis != 'PBMC' & args$pop == 'EUR'){
      afa_corr = corr_matrix[1, 3] %>% pull()
      eur_corr = corr_matrix[2, 3] %>% pull()
      his_corr = corr_matrix[3, 3] %>% pull()
    } else if (args$tis == 'PBMC' & args$pop == 'HIS'){
      afa_corr = corr_matrix[1, 4] %>% pull()
      eur_corr = corr_matrix[2, 4] %>% pull()
      his_corr = corr_matrix[3, 4] %>% pull()
      chn_corr = corr_matrix[4, 4] %>% pull()
    } else if (args$tis != 'PBMC' & args$pop == 'HIS'){
      afa_corr = corr_matrix[1, 4] %>% pull()
      eur_corr = corr_matrix[2, 4] %>% pull()
      his_corr = corr_matrix[3, 4] %>% pull()
    } else if (args$tis == 'PBMC' & args$pop == 'CHN'){
      afa_corr = corr_matrix[1, 5] %>% pull()
      eur_corr = corr_matrix[2, 5] %>% pull()
      his_corr = corr_matrix[3, 5] %>% pull()
      chn_corr = corr_matrix[4, 5] %>% pull()
  }
      
  # load exp values for gene G, and add necessary columns
  afa_tmp_exp_table <- afa_exp %>% filter(gene_id==g) %>% select(-gene_id) %>% t() %>% as.data.frame() %>% mutate(tissue='AFA', sampleid=colnames(afa_exp)[2:ncol(afa_exp)], exp_w=afa_corr, dhs_w=1) %>% 
    rename(exp=V1) %>% select(tissue, sampleid, exp, exp_w, dhs_w)
  eur_tmp_exp_table <- eur_exp %>% filter(gene_id==g) %>% select(-gene_id) %>% t() %>% as.data.frame() %>% mutate(tissue='EUR', sampleid=colnames(eur_exp)[2:ncol(eur_exp)], exp_w=eur_corr, dhs_w=1) %>% 
    rename(exp=V1) %>% select(tissue, sampleid, exp, exp_w, dhs_w)
  his_tmp_exp_table <- his_exp %>% filter(gene_id==g) %>% select(-gene_id) %>% t() %>% as.data.frame() %>% mutate(tissue='HIS', sampleid=colnames(his_exp)[2:ncol(his_exp)], exp_w=his_corr, dhs_w=1) %>% 
    rename(exp=V1) %>% select(tissue, sampleid, exp, exp_w, dhs_w)
    
  if (args$tis != 'PBMC'){
    if (args$pop == 'AFA'){
        # join exp tables with AFA at first position
        final_exp_table <- rbind(afa_tmp_exp_table, eur_tmp_exp_table, his_tmp_exp_table)
      } else if (args$pop == 'EUR'){
        # join exp tables with EUR at first position
        final_exp_table <- rbind(eur_tmp_exp_table, afa_tmp_exp_table, his_tmp_exp_table)
      } else {
        # join exp tables with HIS at first position
        final_exp_table <- rbind(his_tmp_exp_table, afa_tmp_exp_table, eur_tmp_exp_table)
      }
    } else {
      # do same things w/ CHN if it's PBMC
      chn_tmp_exp_table <- chn_exp %>% filter(gene_id==g) %>% select(-gene_id) %>% t() %>% as.data.frame() %>% mutate(tissue='CHN', sampleid=colnames(chn_exp)[2:ncol(chn_exp)], exp_w=chn_corr, dhs_w=1) %>% 
        rename(exp=V1) %>% select(tissue, sampleid, exp, exp_w, dhs_w)
        
      if (args$pop == 'AFA'){
        # join exp tables with AFA at first position
        final_exp_table <- rbind(afa_tmp_exp_table, eur_tmp_exp_table, his_tmp_exp_table, chn_tmp_exp_table)
      } else if (args$pop == 'EUR'){
        # join exp tables with EUR at first position
        final_exp_table <- rbind(eur_tmp_exp_table, afa_tmp_exp_table, his_tmp_exp_table, chn_tmp_exp_table)
      } else if (args$pop == 'HIS') {
        # join exp tables with HIS at first position
        final_exp_table <- rbind(his_tmp_exp_table, afa_tmp_exp_table, eur_tmp_exp_table, chn_tmp_exp_table)
      } else {
        # join exp tables with CHN at first position
        final_exp_table <- rbind(chn_tmp_exp_table, afa_tmp_exp_table, eur_tmp_exp_table, his_tmp_exp_table)
      }
    }
      
    # write final exp table
    fwrite(final_exp_table, 'inputs/exp_tbls/'%&% args$tis %&%'_'%&% args$pop %&%'_'%&% g %&%'_exp_table.txt', sep='\t', col.names=T, quote=F)

    # making JTI command line
    command <- c('Rscript /home/daniel/mashr/MR-JTI/model_training/JTI/JTI.r --tissue='%&% args$pop %&%' --geneid='%&% g %&%' --genotype_path=/home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/plink_files/MESA_TOPMed_WGSX.'%&% args$tis %&%'.'%&% args$pop %&%'.hg38.filt --expression_path=/home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/exp_tbls/'%&% args$tis %&%'_'%&% args$pop %&%'_'%&% g %&%'_exp_table.txt --tmp_folder=/home/daniel/mashr/MR-JTI/TOPMed_MESA/tmp/'%&% args$tis %&%' --gencode_path /home/daniel/mashr/MR-JTI/TOPMed_MESA/inputs/gencode.v38.GRCh38.txt --out_path /home/daniel/mashr/MR-JTI/TOPMed_MESA/outputs/'%&% args$tis %&%'/'%&% args$pop)
    system(command)
    print(args$tis %&% ' - ' %&% args$pop %&% ' - ' %&% g %&% ' - end')
}