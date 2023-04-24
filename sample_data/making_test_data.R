library(tidyverse)
library(data.table)

gene_annot_file <- fread('/home/daniel/mashr/matrixeQTL/WGS_files/gene_annotation_v34_60669.txt.gz', header=T, stringsAsFactors=F) %>% 
  filter(chr == '22', start >= 26429260, end <= 26514568) %>% arrange(start)
fwrite(gene_annot_file, '/home/daniel/github_mashr_project/test_data/gene_annotation.txt', col.names=T, sep=' ', quote=F)
min_dosage_start <- 26429260 - 100000
max_dosage_end <- 26514568 + 100000

dosage_gbr <- fread('/home/daniel/Geuvadis/dosages/Geuvadis.GBR.plink.updated.chr22.dosage.txt.gz')
samples_gbr <- fread('/home/daniel/Geuvadis/GBR_ids_formatted.txt', header=F) %>% pull(V1)
colnames(dosage_gbr) <- c('chr', 'snp_ID', 'pos', 'ref_allele', 'alt_allele', 'maf', all_of(samples_gbr))

dosage_yri <- fread('/home/daniel/Geuvadis/dosages/Geuvadis.YRI.plink.updated.chr22.dosage.txt.gz')
samples_yri <- fread('/home/daniel/Geuvadis/YRI_ids_formatted.txt', header=F) %>% pull(V1)
colnames(dosage_yri) <- c('chr', 'snp_ID', 'pos', 'ref_allele', 'alt_allele', 'maf', all_of(samples_yri))

dosage <- full_join(dosage_gbr, dosage_yri, by=c('snp_ID'))
dosage <- dosage %>% dplyr::select(chr.x, snp_ID, pos.x, ref_allele.x, alt_allele.x, maf.x, maf.y, all_of(samples_gbr), all_of(samples_yri))
colnames(dosage)[1:7] <- c('chr', 'snp_ID', 'pos', 'ref_allele', 'alt_allele', 'maf_gbr', 'maf_yri')
dosage$chr <- 22
dosage <- dosage %>% filter(pos >= min_dosage_start, pos <= max_dosage_end, maf_gbr > 0 | maf_yri > 0)

dosage_unique_snps <- dosage %>% select(pos) 
dosage_unique_snps <- dosage_unique_snps[!(duplicated(dosage_unique_snps) | duplicated(dosage_unique_snps, fromLast=T)),] # removes multiallelic snps 
dosage <- left_join(dosage_unique_snps, dosage)

dosage_gbr_raw <- dosage %>% dplyr::select(chr, snp_ID, pos, ref_allele, alt_allele, all_of(samples_gbr))
fwrite(dosage_gbr_raw, '/home/daniel/github_mashr_project/sample_data/dosages/GEUVADIS_GBR_chr22_dosage_unfiltered.txt', col.names=T, quote=F, sep=' ')

dosage_gbr_filt <- dosage %>% filter(maf_gbr > 0) %>% dplyr::select(chr, snp_ID, pos, ref_allele, alt_allele, all_of(samples_gbr))
fwrite(dosage_gbr_filt, '/home/daniel/github_mashr_project/sample_data/dosages/GEUVADIS_GBR_chr22_dosage_filtered.txt', col.names=T, quote=F, sep=' ')

dosage_yri_raw <- dosage %>% dplyr::select(chr, snp_ID, pos, ref_allele, alt_allele, all_of(samples_yri))
fwrite(dosage_yri_raw, '/home/daniel/github_mashr_project/sample_data/dosages/GEUVADIS_YRI_chr22_dosage_unfiltered.txt', col.names=T, quote=F, sep=' ')

dosage_yri_filt <- dosage %>% filter(maf_yri > 0) %>% dplyr::select(chr, snp_ID, pos, ref_allele, alt_allele, all_of(samples_yri))
fwrite(dosage_yri_filt, '/home/daniel/github_mashr_project/sample_data/dosages/GEUVADIS_YRI_chr22_dosage_filtered.txt', col.names=T, quote=F, sep=' ')

exp_tbl_gbr <- fread('~/Geuvadis/rna_counts/GBR_10PCAIR_PF_adj_rinv_TOPMED_expression10_peer_factor_adjusted_filt.txt.gz') %>%
  filter(gene_id %in% gene_annot_file$gene_id)
fwrite(exp_tbl_gbr, '/home/daniel/github_mashr_project/sample_data/exp_tbls/GEUVADIS_GBR_expression.txt', col.names=T, quote=F, sep=' ')

exp_tbl_yri <- fread('~/Geuvadis/rna_counts/YRI_10PCAIR_PF_adj_rinv_TOPMED_expression10_peer_factor_adjusted_filt.txt.gz') %>%
  filter(gene_id %in% gene_annot_file$gene_id)
fwrite(exp_tbl_yri, '/home/daniel/github_mashr_project/sample_data/exp_tbls/GEUVADIS_YRI_expression.txt', col.names=T, quote=F, sep=' ')