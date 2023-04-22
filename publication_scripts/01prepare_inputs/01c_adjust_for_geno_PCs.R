library(dplyr)
library(data.table)
library(tibble)
'%&%' = function(a,b) paste(a,b,sep="")

rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

pops<-c("AFA","ALL","CAU","HIS","CHN")
for ( pop in pops ){
  cat(pop,'\n')
  expression <- fread("/home/chris/topmed_expression_whole_genome/expression/PBMC_expression_" %&% pop %&% "age_sex_adj.txt",header=T)

  #print(head(expression))
  
  #pcs <- fread("/home/chris/topmed_expression_whole_genome/02_QC_PCAIR/" %&% pop %&% "/PBMC/QC/PCA/unmerged_pca.eigenvec",header=TRUE)
  pcs <- fread("/home/chris/topmed_expression_whole_genome/02_QC_PCAIR/" %&% pop %&% "/PBMC/PCAIR/PCAIR.eigenvec",header=TRUE)
  
  #print(head(pcs))
 
  #pcmat <- as.matrix(pcs[,-1:-2])
  pcmat <- as.matrix(pcs[,2:11])
  
  #print(head(pcmat))
  
  print(length(pcs$sample_id))
  print(length(expression$sidno))
  finaldf <-expression %>% dplyr::filter(sidno %in% pcs$sample_id) 
  ordered_pcs<-finaldf[,1:ncol(pcs)] %>% select(-"sidno")
  sidno <-finaldf %>% select("sidno")
  finaldf<-finaldf %>% column_to_rownames(var="sidno")
  # str(finaldf)
  # 
  #print(colnames(finaldf))
  #ordered_pcs<-finaldf[,1:ncol(pcs)] %>% select(-"sidno")
  
  lmfunc <- function(x){resid(lm(x ~ .,data=ordered_pcs))}
  adjmat <- apply(finaldf, 2, lmfunc)
  adjmat <- apply(adjmat, 2, rankinv)
  adjmat<- cbind.data.frame(sidno,adjmat)
  str(adjmat)

  fwrite(adjmat, file = "/home/chris/topmed_expression_whole_genome/expression/PBMC_expression_" %&% pop %&% "age_sex_adj_rinv_PC10.txt",quote=F,sep="\t")
}

