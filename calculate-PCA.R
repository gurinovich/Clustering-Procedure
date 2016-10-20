#Input parameters - change them accordingly:

#work.dir is where we want our results to be saved
work.dir <- "/home/user_name/project/"
cluster <- "name of the cluster - save it in a file with this cluster name preceeded by cluster. and .txt in the end"
#cluster.cluster-name.txt is a tab-delimited file with two columns: Family and ID, txt format"
clust.dir <- "directory where cluster file is located"
geneticData <- "link to genetic data in PLINK bed/bim/fam format"
EIGENSTRATscript <- "link to a file for the first step in the EIGENSTRAT analysis - can be empty - example full name: par.PED.EIGENSTRAT"
IBD <- "link to a file with IBD estimations for the data"
smartpca <- "link to a file for the second step in the EIGENSTRAT analysis - can be empty - example full name: smartpca.pa"

pca.dir <- paste0(work.dir,"PCA/")
analysis.dir <- paste0(work.dir,"analysis/")

#QC:
system(paste0("mkdir ", pca.dir, cluster),wait=T)
system(paste0("mkdir ", pca.dir, cluster,"/QC"),wait=T)
setwd(paste0(pca.dir, cluster,"/QC"))
system(paste0("cp ", clust.dir, "cluster.", cluster, ".txt ."),wait=T)
system(paste0("plink --bfile ", geneticData, " --keep ", "cluster.", cluster, ".txt " ,"--make-bed --out all_QC_pruned"),wait=T)
fam <- read.table("all_QC_pruned.fam")
fam[,6]  <- 1
write.table(fam,"all_QC_pruned_edit.fam", quote = F, col.names = F, row.names = F) 
#Convert:
system(paste0("mkdir ", pca.dir, cluster,"/Convert"),wait=T)
setwd(paste0(pca.dir, cluster,"/Convert"))
system(paste0("cp ", EIGENSTRATscript, " ."),wait=T)
fileConn <- file("par.PED.EIGENSTRAT","w")
writeLines(c(paste0("genotypename: ",pca.dir, cluster,"/QC/all_QC_pruned.bed #.ped, .bed "),
             paste0("snpname: ",pca.dir, cluster,"/QC/all_QC_pruned.bim  # .map or .bim, either works "),
             paste0("indivname: ",pca.dir, cluster,"/QC/all_QC_pruned_edit.fam # .ped or .fam, either works"),
             "outputformat: EIGENSTRAT",
             "genotypeoutname: all.eigenstratgeno",
             "snpoutname: all.snp",
             "indivoutname: all.ind",
             "familynames: NO"), fileConn)
close(fileConn)

system("convertf -p par.PED.EIGENSTRAT > convert.log",wait=T)
#Eigen:
system(paste0("mkdir ", pca.dir, cluster,"/Eigen"),wait=T)
setwd(paste0(pca.dir, cluster,"/Eigen"))

ind <- read.table(paste0(pca.dir, cluster,"/Convert/all.ind"))

#########################
### make list of related individuals and then set to related
ibd <- read.table(IBD, header = T)

which.ibd.keep <- which(as.character(ibd[,2]) %in% ind[,1] & 
                          as.character(ibd[,4]) %in% ind[,1] & ibd$PI_HAT > 0.2)
length(which.ibd.keep)

if (length(which.ibd.keep) == 0) {
  ind[,3] <- "Phenotype"
  
} else {
  
  
  ibd.keep <- ibd[which.ibd.keep,]
  table(table(c( as.character(ibd.keep[,2]), as.character(ibd.keep[,4]))))	 
  
  related.ids <- c()
  for (i in 1:dim(ibd.keep)[1]){
    if ( !(as.character(ibd.keep[i,2]) %in% related.ids) & 
         !(as.character(ibd.keep[i,4]) %in% related.ids)){
      ran.index <- sample(c(2,4),1)
      related.ids <- c(related.ids, as.character(ibd.keep[i,ran.index]))
    }
  }
  write.csv(related.ids, "Related.IDS.csv", quote = F, row.names = F)
  which.rel <- which(ind[,1] %in% related.ids)
  length(which.rel)
  ind[,3] <- "Phenotype"
  ind[which.rel,3] <- "Related_InferredPC"
}

write.table(ind , "phenotype_PCA.ind", quote = F, col.names = F, row.names = F, sep = "    ")

## making a list of populations to be analyzed (don't put related in there b/c you want their eigenvectors inferred)
pop <- names(table(ind[,3]))
pop <- pop[!(pop %in% c("Ignore", "Related_InferredPC","BorderlineCR_InferredPC" ))]
write.table(pop, "pop_list.txt", quote = F, col.names = F, row.names = F)

system(paste0("cp ", smartpca, " ."),wait=T)

fileConn <- file("smartpca.par","w")
writeLines(c(paste0("genotypename: ",pca.dir, cluster,"/Convert/all.eigenstratgeno"),
             paste0("snpname: ",pca.dir, cluster,"/Convert/all.snp"),
             paste0("indivname: ",pca.dir, cluster,"/Eigen/phenotype_PCA.ind"),
             
             "evecoutname:     phenotype.evec",
             "evaloutname:     phenotype.eval",
             "numoutevec:      30",
             "numoutlierevec:  5",
             "outliersigmathresh:  9",
             paste0("poplistname: ",pca.dir, cluster,"/Eigen/pop_list.txt"),
             "snpweightoutname: phenotype_snp_wt.txt"), fileConn)
close(fileConn)

system("smartpca -p smartpca.par > pca.log",wait=T)

system(paste0("sed '1d' ", pca.dir, cluster, "/Eigen/phenotype.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | cut -f 1-31 -d \" \" > ", analysis.dir, cluster, ".txt"),wait=T)
