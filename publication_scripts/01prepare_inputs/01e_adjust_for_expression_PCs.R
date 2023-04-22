library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
'%&%' = function(a,b) paste(a,b,sep="")

pops<-c("ALL")
tiss<-c("Mono")
ranks<-c(10,20,30)
for (rank in ranks) {
  cat(rank, '\n')
  for (tis in tiss){
    for ( pop in pops ){
      cat(pop,'\n')
      #expression <- as.data.frame(fread("/home/chris/topmed_expression_whole_genome/expression/" %&% tis %&% "_expression_" %&% pop %&% "age_sex_adj_PC10.txt",header=T))  
      #expression <- fread("/home/chris/topmed_expression_whole_genome/expression/Mono_expression_" %&% pop %&% "age_sex_adj.txt",header=T)
      expression <- fread("/home/chris/topmed_expression_whole_genome/expression/" %&% tis %&% "_expression_" %&% pop %&% "age_sex_adj_rinv_PC10.txt",header=T)
      #print(head(expression))
  
      #pcs <- fread("/home/chris/topmed_expression_whole_genome/02_QC_PCAIR/" %&% pop %&% "/Mono/QC/PCA/unmerged_pca.eigenvec",header=TRUE)
      #pcs <- fread("/home/chris/topmed_expression_whole_genome/02_QC_PCAIR/" %&% pop %&% "/Mono/PCAIR/PCAIR.eigenvec",header=TRUE)
      pcs <- fread("/home/chris/topmed_expression_whole_genome/expression/PCA/" %&% tis %&% "_" %&% pop %&% "age_sex_adj_PC10_expression_PCs" %&% rank %&% ".txt",header=TRUE)
      #print(head(pcs))
 
      #pcmat <- as.matrix(pcs[,-1:-2])
      print(str(pcs))
      if(rank==10){
      pcmat <- as.matrix(pcs[,2:11])
      }
      if(rank==20){
      pcmat <- as.matrix(pcs[,2:21])
      }
      if(rank==30){
      pcmat <- as.matrix(pcs[,2:31])
      }
      #pcmat_20 <- as.matrix(pcs[,2:21])
      #pcmat_30 <- as.matrix(pcs[,2:31])
      print(head(pcmat))
  
      print(length(pcs$sample_id))
      print(length(expression$sidno))
      finaldf <-expression %>% dplyr::filter(sidno %in% pcs$sample_id) 
      sidno <-finaldf %>% select("sidno")
      finaldf<-finaldf %>% column_to_rownames(var="sidno")
      # str(finaldf)
      # 
      str(sidno)
      lmfunc <- function(x){resid(lm(x ~ pcmat))}
      adjmat <- apply(finaldf, 2, lmfunc)
      adjmat<- cbind.data.frame(sidno,adjmat)
      str(adjmat)

      fwrite(adjmat, file = "/home/chris/topmed_expression_whole_genome/expression/" %&% tis %&% "_expression_" %&% pop %&% "age_sex_adj_rinv_PC10_expression_PCs" %&% rank %&% ".txt",quote=F,sep="\t")
    }
  }
}
