## Script to run PC-Air

# script based on Ryan Schubert's PCAir scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/02PCAIR

# Loading libraries
library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(dplyr)
library(tibble)
library(data.table)
'%&%' = function(a,b) paste (a,b,sep='')

for (tissue in c('AFA','EUR','HIS','CHN')){
  # 01. Make GDS
  snpgdsBED2GDS(bed.fn="/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/MESA_auto_filt.bed",
                bim.fn="/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/MESA_auto_filt.bim",
                fam.fn="/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/MESA_auto_filt.fam",
                out.gdsfn="/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/00autosome.gds")
  
  # 02. King estimation
  gdsfile<-"/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/00autosome.gds"
  gds <- snpgdsOpen(gdsfile)
  king<-snpgdsIBDKING(gds)
  kingMat<-king$kinship
  snpgdsClose(gds)
  saveRDS(king,file="/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/King_matrix.RDS")
  
  # 03. LD prune
  gdsfile<-"/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/00autosome.gds"
  gds <- snpgdsOpen(gdsfile)
  snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
                            ld.threshold=sqrt(0.3), verbose=FALSE)
  pruned <- unlist(snpset, use.names=FALSE)
  saveRDS(pruned,"/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/pruned_set.RDS")
  snpgdsClose(gds)
  
  # 04. Run PCAIR
  gdsfile<-"/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/00autosome.gds"
  pruned<-readRDS("/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/pruned_set.RDS")
  king<-readRDS("/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/King_matrix.RDS")
  kingMat<-king$kinship
  colnames(kingMat)<-king$sample.id
  row.names(kingMat)<-king$sample.id
  geno <- GdsGenotypeReader(filename = gdsfile)
  genoData <- GenotypeData(geno)
  mypcair <- pcair(genoData, kinobj = kingMat, divobj = kingMat,
                   snp.include = pruned)
  png("/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/PCAIR_PC1_VS_PC2.png")
  plot(mypcair)
  dev.off()
  eigenvec<-mypcair$vectors %>% as.data.frame() %>% rownames_to_column(var="sample_id")
  str(eigenvec)
  val<-mypcair$values %>% as.data.frame()
  fwrite(eigenvec,"/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/PCAIR.eigenvec",col.names = T,row.names = F,sep='\t')
  fwrite(val,"/home/daniel/MESA_genotypes_subset/PCAir/"%&% tissue %&%"_pop/PCAIR.eigenval")
} 