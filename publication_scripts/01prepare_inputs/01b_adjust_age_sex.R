## 00 make data adjusted for age sex and exam

# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/03adjust_expression

#################################################
# SET UP ENVIRONMENT
#################################################

library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(viridis)
"%&%" = function(a,b) paste(a,b,sep="")

out.dir = "/home/daniel/MESA_genotypes_subset/Elastic_Net/adjusted_expression/"


################################################
# READ PROTEIN DATA
################################################

protdf <- fread("/home/daniel/TOPMed_MESA_RNAseq/tpm_files/tpm_subsets/TOPMed_MESA_RNAseq_WithDemographicInfo_mean01filt_add106.txt",header=TRUE)
protdf <- as_tibble(protdf)

afa1df <- filter(protdf,pop=="AFA",exam==1,complete.cases(protdf))
afa5df <- filter(protdf,pop=="AFA",exam==5,complete.cases(protdf))
eur1df <- filter(protdf,pop=="EUR",exam==1,complete.cases(protdf))
eur5df <- filter(protdf,pop=="EUR",exam==5,complete.cases(protdf))
chn1df <- filter(protdf,pop=="CHN",exam==1,complete.cases(protdf))
chn5df <- filter(protdf,pop=="CHN",exam==5,complete.cases(protdf))
his1df <- filter(protdf,pop=="HIS",exam==1,complete.cases(protdf))
his5df <- filter(protdf,pop=="HIS",exam==5,complete.cases(protdf))
all1df <- filter(protdf,exam==1,complete.cases(protdf)); #print(all1df[!complete.cases(all1df),1])
all5df <- filter(protdf,exam==5,complete.cases(protdf)); #print(all5df[!complete.cases(all5df),1])
print(all1df)


###############################################
# ADJUST
###############################################

library(preprocessCore) #has normalize.quantiles function
rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

adjmatlist <- list() #list to store adjmat's
adjdflist <- list() #list to store adjdf's
for(exam in c("1","5")){
  for ( pop in c("AFA","HIS","ALL","EUR"))
  {
    df <- get(tolower(pop) %&% exam %&% "df")
    #print(all(is.na(df)))
    #print(df)
    rawmat <- as.matrix(df[,6:ncol(df)]); #print(dim(rawmat)
    #print(rawmat)
    #print(all(is.na(rawmat)))
    logmat <- log(rawmat); #print(dim(logmat)) #natural log transform
    #print(all(is.na(logmat)))
    lmfunc <- function(x){resid(lm(x ~ df$age + df$gender_code))} #get residuals of prot after adj for age & sex
    adjmat <- apply(logmat, 2, lmfunc) ; #print(dim(adjmat)) #apply lmfunc to each column of logmat
    name <- pop %&% exam
    adjmatlist[[name]] <- adjmat
    adjdf <- cbind(df[,1], adjmat)
    adjdflist[[name]] <- adjdf
  }
}

elementMean <- function(my.list) { 
  arr <- array(unlist(my.list),c(dim(my.list[[1]])[1],dim(my.list[[1]])[2],length(my.list)))
  rowMeans( arr , na.rm=TRUE, dims = 2 )
}

# print(str(adjdflist))
#full join to add NA's to df if missing Exam 1 or 5
afadf <- full_join(adjdflist$AFA1, adjdflist$AFA5, by = "sidno"); print("afa")
eurdf <- full_join(adjdflist$EUR1, adjdflist$EUR5, by = "sidno"); print("eur")
#chndf <- full_join(adjdflist$CHN1, adjdflist$CHN5, by = "sidno"); print("chn")
hisdf <- full_join(adjdflist$HIS1, adjdflist$HIS5, by = "sidno"); print("his")
alldf <- full_join(adjdflist$ALL1, adjdflist$ALL5, by = "sidno"); print("all")

#################################
# FINAL ADJUSTMENTS AND WRITE
#################################

for(pop in c('AFA','EUR','HIS','ALL')){
  df <- get(tolower(pop) %&% "df")
  df1na <- select(df,ends_with(".x"))
  df5na <- select(df,ends_with(".y"))
  meanmat <- elementMean(list(df1na, df5na)) ##take the mean of exam 1 and exam 5
  invmeanmat <- round(apply(meanmat, 2, rankinv),6) ##rank-inverse normalize and round
  finaldf <- cbind(df[,1],as.data.frame(invmeanmat)) ##add sidno to df
  colnames(finaldf) <- c("sidno", colnames(protdf[6:ncol(protdf)])) ##retrieve column names
  fwrite(finaldf, file = out.dir %&% "RNASeq_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse.txt",quote=F,sep="\t")
}
