#Input parameters - change them accordingly:

#work.dir is where we want our results to be saved
work.dir <- "/home/user_name/project/"
analysis.dir <- "link to the directory where .txt files with PCs for each cluster are located; with columns (no names): ID, PC1, PC2,..."
mega.data <- "link to a .csv file with columns: phenotype, ID, Family, Sex, SNP-dosages"
SNPs <- "link to a .txt file with the list of SNPs to investigate"

system(paste0("mkdir ", work.dir, "full-summary"))
res.dir <- paste0(work.dir,"full-summary/")

PCsNumber <- 6

mega.data <- read.csv(mega.data,header=T)
mega.data <- subset(mega.data,select=-c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8))

snps <- read.csv(SNPs,header=F,as.is=T)
snps <-snps$V1

file.names <- dir(analysis.dir)

for (i in 1:length(snps)) {
  output <- data.frame(matrix(NA,nrow=length(file.names),ncol=16))
  names(output) <- c("Clusters","PC1.b","PC1.p","PC2.b","PC2.p","PC3.b","PC3.p","PC4.b","PC4.p","PC5.b","PC5.p","PC6.b","PC6.p","GLM.b","GLM.p","GLM.SE")
  output$Clusters <- gsub(".txt","",file.names)
  
  for (j in 1:length(file.names)) {
    
    cluster <- read.table(paste0(analysis.dir,file.names[j]))
    pca.names <- c("ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30") 
    colnames(cluster) <- pca.names
    mega.data.cluster <- merge(cluster,mega.data,by="ID")
    
    
    for (k in 1:PCsNumber) {
      mod <- lm(paste0("PC",k,"~Sex+",snps[i]), data=mega.data.cluster)
      output[j,2*k] <-  as.numeric(coef(summary(mod))[,1][3])  #betta
      output[j,2*k+1] <-  as.numeric(coef(summary(mod))[,4][3])  #p-value
    }
    
    mod <- glm(paste0("phenotype~Sex+PC1+PC2+PC3+PC4+PC5+PC6+",snps[i]), family=binomial, data=mega.data.cluster)
    output[j,14] <-  as.numeric(coef(summary(mod))[,1][9])  #betta
    output[j,15] <-  as.numeric(coef(summary(mod))[,4][9])  #p-value
    output[j,16] <-  as.numeric(coef(summary(mod))[,2][9])   #SE for the betta
    
  }
  
  write.table(output,paste0(res.dir,snps[i],".txt"),quote=F,sep=",",row.names=F)
}

