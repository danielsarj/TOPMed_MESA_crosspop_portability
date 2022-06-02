# script based on Ryan Schubert's scripts, which can be found in here:
# https://github.com/RyanSchu/TOPMed_Proteome/tree/main/07Valid_CV

##################################################################################################

#SETUP ENVIRONMENT
cat("SETUP ENVIRONMENT\n")
##################################################################################################

library(dplyr)
library(qvalue)
library(data.table)
library(argparse) 


##################################################################################################

#DEFINE FUNCTIONS
cat("DEFINE FUNCTIONS\n")
##################################################################################################


"%&%" = function(a,b) paste(a,b,sep="")

match_orders<-function(rowmatrix,colmatrix, idcol=2){
  sampleOrderColMat<-colmatrix[,idcol] %>% unlist() %>% unname()
  col1RowMat<-colnames(rowmatrix)[1]
  rowmatrix<-rowmatrix %>% select(one_of(c(col1RowMat,sampleOrderColMat)))
  return(rowmatrix)
}

get_gene_corr<-function(gene_name,colmatrix,rowmatrix,skipcol=1)
{
  #gene name is the name of the gene
  #colmatrix is a matrix containing genes as columns
  #rowmatrix is a matrix containing genes as rows
  #skipcol = The first n columns of the row matrix do not contain values, skip these for correlation test  
  
  predicted_exp<-colmatrix %>% 
    select(gene_name) %>% ##select the gene
    unlist() %>% unname()
  
  measured_exp<-rowmatrix %>% 
    filter( gene_id == gene_name) %>% 
    select((skipcol+1):ncol(rowmatrix)) %>%
    unlist() %>% 
    unname()
  
  correlation<-cor.test(measured_exp,predicted_exp, method = "spearman")
  
  expression_corr<-list()
  expression_corr[["gene_id"]]<-gene_name
  expression_corr[["estimate"]]<-correlation$estimate
  expression_corr[["p.value"]]<-correlation$p.value
  
  return(expression_corr)
}


##################################################################################################

#SET GLOBAL VARIABLES
cat("SET GLOBAL VARIABLES\n")
##################################################################################################

parser <- ArgumentParser()
parser$add_argument('--tissue', help='tissue/cell type prediction models')
parser$add_argument('--model', help='if results are from EN, ENunf, mashr or matrixeqtl')
parser$add_argument('--model_pop', help='model population')
parser$add_argument('--predict_pop', help='predicted population')
args <- parser$parse_args()

predtiss <- args$tissue
model <- args$model
pop_list <- args$model_pop
obspop_list <- args$predict_pop

pip_list<-c("0")
clus_list<-c("0")
R2_filtering<-c(-1)

columns<-expand.grid(pop_list,clus_list,pip_list) %>% 
  arrange(Var1,Var2,Var3) %>%
  mutate(columns=Var1 %&% Var2 %&% Var3) %>% rename(pop=Var1,clus=Var2,pip=Var3)

col_list<-columns %>% select(columns) %>% unlist() %>% unname()

pi1_matrix<-matrix(NA, nrow = length(obspop_list), ncol = length(col_list))
rownames(pi1_matrix)<-obspop_list
colnames(pi1_matrix)<-col_list

##################################################################################################

#READ IN & PROCESS DATA
cat("READ IN & PROCESS DATA\n")
##################################################################################################

for (obspop in obspop_list)
{#per each set of observed expression data
  
  observed_expression<-fread("/home/daraujo1/data/Geuvadis/rna_counts/" %&% obspop %&% "_10PCAIR_PF_adj_rinv_TOPMED_expression10_peer_factor_adjusted_filt.txt.gz", header = T, sep = '\t',stringsAsFactors = F) %>%
    rename_at(vars(1),function(x){return("gene_id")}) %>% select("gene_id", sort(colnames(.)))
  
  for (i in 1:nrow(columns))
  {#How well does each model replicate
    
    predpop<-columns[i,1]
    predclus<-columns[i,2]
    predpip<-columns[i,3]
    pi1column<-columns[i,5]
    
    models<-fread("/home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis." %&% obspop %&% "_" %&% predtiss %&% "." %&% predpop %&% "_" %&% model %&% "_summary_predict.txt", header = F, sep = '\t',stringsAsFactors = F)
    #names<-fread("/home/daniel/Geuvadis/MESA_models_validation/geno_prediction_summary_header.txt",header=F,sep = '\t',stringsAsFactors = F)
    #names<-c(unlist(names))
    #colnames(models)<-c("gene_id","genename","gene_type","alpha","n.snps.in.window","n.snps.in.model","lambda_min_mse",
    #"test_R2_avg","test_R2_sd","cv_R2_avg","cv_R2_sd","in_sample_R2","nested_cv_fisher_pval","rho_avg","rho_se","rho_zscore",
    #"pred.perf.R2","pred.perf.pval","cv_rho_avg","cv_rho_se","cv_rho_avg_squared","cv_zscore_est","cv_zscore_pval","cv_pval_est","pred.perf.qval")
    colnames(models)<-c("gene_id","gene_name","n_snps_in_model","n_snps_used","test_R2_avg","pred_perf_pval")
    
    models<-models %>% select(gene_id, test_R2_avg)
    
    predicted_expression<-fread("/home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis." %&% obspop %&% "_" %&% predtiss %&% "." %&% predpop %&% "_" %&% model %&% "_predict.txt", header = T, sep = '\t',stringsAsFactors = F)
    predicted_expression <- predicted_expression[order(predicted_expression$FID),]

    obsGenes<-unlist(observed_expression$gene_id)
    predGenes<-colnames(predicted_expression)
    
    gene_list<-data.frame(gene_id=intersect(obsGenes,predGenes),stringsAsFactors = F)
    
    models<-inner_join(gene_list,models, by = "gene_id")
    cat("Existing Gene Models:\n")
    str(models)
    
    ##################################################################################################
    
    #ANALYZE DATA
    cat("ANALYZE DATA\n")
    ##################################################################################################
    
    for(j in R2_filtering)
    {#Calculate pi1 at different model R2 thresholds
      
      cat("testing replication rate at R2 threshold of ",j,"\n")
      filtered_gene_R2<-models
      
      if(dim(filtered_gene_R2)[1] == 0)
      {#check if there are any models that meet R2 threshold
        
        cat("WARNING: ",pi1column," had no models at R2 threshold ", j, " for genes in ",obspop,"\n")
        #        cat("WARNING: possible at higher R2 thresholds, but unlikely at low thresholds unless there is a gene name mismatch. Check that your genenames are the same between files")
        cat("WARNING: ASSIGNING 0 as pi1\n")
        print(pi1_matrix)
        pi1_matrix[obspop,pi1column] <- 0
        next
        
      }
      predictive_correlations<-sapply(X=gene_list$gene_id,FUN=get_gene_corr,colmatrix=predicted_expression,rowmatrix=observed_expression,simplify=T,USE.NAMES = T)
      predictive_correlations<-data.frame(gene_id=unlist(predictive_correlations[1,]),estimate=unlist(predictive_correlations[2,]),p.value=unlist(predictive_correlations[3,]))
      
      str(predictive_correlations)
      
      filtered_gene_R2<- filtered_gene_R2 %>% inner_join(predictive_correlations, by = "gene_id")
      
      ##################################################################################################
      
      #WRITE OUT CORRELATION DATA
      cat("WRITE OUT CORRELATION DATA\n")
      ##################################################################################################
      
      fwrite(filtered_gene_R2,"/home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis." %&% obspop %&% "_" %&% predtiss %&% "." %&% predpop %&% "_" %&% model %&% "_expression_spearman_correlation.txt",col.names = T,sep='\t')
      corr_pvals<-filtered_gene_R2$p.value
      qobjCorr <- tryCatch({qvalue(p = corr_pvals)},
                           error=function(cond){
                             cat("Error: ", pi1column, "and ",obspop, " caused error to qvalue, assigning 0 as pi1\n")
                             cond
                           })
      if(inherits(qobjCorr, "error")) {
        pi1_matrix[obspop,pi1column] <- 0
        next
      }
      pi1<- 1 - qobjCorr$pi0
      pi12<- signif(pi1,4)
      #      pi1_matrix[obspop,pi1column] <- pi12
    }#close
  }
}

##################################################################################################

#WRITE OUT SUMMARY DATA
cat("WRITE OUT SUMMARY DATA")
##################################################################################################


pi1_matrix<-data.table(pi1_matrix)
row.names(pi1_matrix)<-obspop_list
#fwrite(pi1_matrix, "/home/daraujo1/data/Geuvadis/WGS_predixcan/Geuvadis." %&% obspop %&% "_" %&% predtiss %&% "." %&% predpop %&% "_" %&% model %&% "_pi1_R2-1.txt", sep = '\t', col.names=T,row.names = T)


