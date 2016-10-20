#Input parameters - change them accordingly:

#work.dir is where we want our results to be saved
work.dir <- "/home/user_name/project/"
PC.pvalue <- 0.05/6
SNP.pvalue <- 0.05

res.dir <- paste0(work.dir,"full-summary-yes-no/")
files.dir <- paste0(work.dir,"full-summary/")

system(paste0("mkdir ", work.dir, "full-summary-yes-no/"))

file.names <- dir(files.dir)

for (i in 1:length(file.names)) {
  input <- read.csv(paste0(files.dir,file.names[i]),colClasses=c("Clusters"="character"))
  output <- data.frame(matrix(NA,nrow=nrow(input),ncol=3))
  names(output) <- c("Clusters","PC1.PC6.signif","SNP.signif")
  output$Clusters <- input$Clusters
  output$PC1.PC6.signif <- ifelse(input$PC1.p <= PC.pvalue | input$PC2.p <= PC.pvalue | input$PC3.p <= PC.pvalue | input$PC4.p <= PC.pvalue | input$PC5.p <= PC.pvalue | input$PC6.p <= PC.pvalue,"yes","no")
  output$SNP.signif <- ifelse(input$GLM.p <= SNP.pvalue,"yes","no")
  write.table(output,paste0(res.dir,file.names[i]),quote=F,sep=",",row.names=F)
}