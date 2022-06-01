library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
'%&%' = function(a,b) paste(a,b,sep="")

pops<-c("AFA")
tiss<-c("Mono")
ranks<-c(10,20,30)
for (rank in ranks) {
  cat(rank, '\n')
  for (tis in tiss){
    for ( pop in pops ){
      cat(pop,'\n')
      raw_expression <- as.data.frame(fread("/home/chris/topmed_expression_whole_genome/expression/" %&% tis %&% "_expression_" %&% pop %&% "age_sex_adj_rinv_PC10.txt",header=T))  
      #print(head(raw_expression))
      
      # Adjusting Matrix: remove columns with all zeros, transpose, 1st row as column names
      #raw_expression <- raw_expression[,colSums(raw_expression) > 0, drop=FALSE]
      expression<-as.data.frame(t(raw_expression))
      #print(head(expression))
      colnames(expression) <- expression[1, ]
      expression <- expression[-1,]
      print(head(expression))
      # Perform PCA
      expression.pca <- prcomp(expression,scale. = TRUE,rank.=rank)
      
      # Write PCs to file 
      PCs <- as.data.frame(expression.pca$rotation)
      PCs <- tibble::rownames_to_column(PCs, "sample_id")
      fwrite(PCs, "/home/chris/topmed_expression_whole_genome/expression/PCA/" %&% tis %&% "_" %&% pop %&% "age_sex_adj_PC10_expression_PCs" %&% rank %&% ".txt",quote=F,sep="\t")
      
      # Plot Percent of Variance Explaied  
      pov <- as.data.frame(summary(expression.pca)$importance[2,])
      pov <- tibble::rownames_to_column(pov, "PC")
      colnames(pov) <- c("PC","pov")
      pov$PC <- factor(pov$PC, levels = pov$PC)  
      print(pov)
      if(rank==10){
      pov_plot <- pov[1:10,]
      }
      if(rank==20){
      pov_plot <- pov[1:20,]
      }
      if(rank==30){
      pov_plot <- pov[1:30,]
      }
      print(pov_plot)
      plot_pov <-ggplot(pov_plot,aes(x=PC,y=pov)) + geom_bar(stat="identity") + scale_y_continuous(limits = c(0,0.8)) +
        ggtitle("Percent of Variance Explained PCs " %&% tis %&% " " %&% pop %&% "") + xlab("PCs") + ylab("Percent Variance")
      ggsave("/home/chris/topmed_expression_whole_genome/expression/PCA/" %&% tis %&% "_" %&% pop %&% "age_sex_adj_PC10_expression_PCs" %&% rank %&% "_percent_variance.png")
      
      # Plot PC1 v PC2
      cat("Loading in Sample Pop")
      sample_id <- as.data.frame(fread("/home/chris/topmed_expression_whole_genome/expression/sample_pop/" %&% tis %&% "_pop.txt"))
      colnames(sample_id) <- c("sample","pop")
      
      cat("PC1 v PC2")
      pc1_pc2 <- as.data.frame(expression.pca$rotation)[,1:2]
      pc1_pc2 <- tibble::rownames_to_column(pc1_pc2, "sample")
      sample_id$sample <- as.integer(sample_id$sample)
      pc1_pc2$sample <- as.integer(pc1_pc2$sample)
      pc1_pc2 <- inner_join(pc1_pc2,sample_id, by="sample")
      plot_pc1_pc2 <- ggplot(pc1_pc2,aes(x=PC1,y=PC2,color=pop)) + geom_point() +
        ggtitle("PC1 v PC2 " %&% tis %&% " " %&% pop %&% "") 
      ggsave("/home/chris/topmed_expression_whole_genome/expression/PCA/" %&% tis %&% "_" %&% pop %&% "age_sex_adj_PC10_expression_PCs" %&% rank %&% "_PC1_PC2.png")
    }
  }
}
