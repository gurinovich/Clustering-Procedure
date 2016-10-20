#Input parameters - change them accordingly:

#work.dir is where we want our results to be saved
work.dir <- "/home/user_name/project/"
files.YN.dir <- paste0(work.dir,"full-summary-yes-no/")
files.PV.dir <- paste0(work.dir,"full-summary/")

file.names <- dir(files.YN.dir)
snps <- gsub(".txt","",file.names)

output <- data.frame(matrix(NA,nrow=1,ncol=5))
names(output) <- c("SNP","Clusters","Betta","SE","P.value")

for (i in 1:length(file.names)) {
  input <- read.csv(paste0(files.YN.dir,file.names[i]),colClasses="character")
  temp <- input[which(input$PC1.PC6.signif == "no" & input$SNP.signif == "yes"),]
  if (nrow(temp) != 0) {
    input.PV <- read.csv(paste0(files.PV.dir,file.names[i]),colClasses=c("Clusters"="character"))
    betta <- input.PV$GLM.b[which(input.PV$Clusters %in% temp$Clusters)]
    se <- input.PV$GLM.SE[which(input.PV$Clusters %in% temp$Clusters)]
    p.values <- input.PV$GLM.p[which(input.PV$Clusters %in% temp$Clusters)]
    output.add <- cbind(rep(snps[i],nrow(temp)),temp$Clusters,betta,se,p.values)
    colnames(output.add) <- c("SNP","Clusters","Betta","SE","P.value")
    output <- rbind(output,output.add)
  }
}

output$Betta <- as.numeric(output$Betta)
output$Betta <- -output$Betta

temp.df <- output[-which(output$P.value == 0),]
if ( nrow(temp.df) != 0 )
  output <- output[-which(output$P.value == 0),]

#remove rows with NAs:
output <- output[complete.cases(output),]

#sort by p-values:
output <- output[order(output$P.value),]

output$SNP <- as.character(output$SNP)

write.table(output,paste0(work.dir,"summary_yes_no.csv"),quote=F,sep=",",row.names=F)
