#Input parameters - change them accordingly:

#work.dir is where we want our results to be saved
work.dir <- "/home/user_name/project/"
clust.dir <- "link to the location of files with clusters info: Family & Sample Ids - derived from the previous step - see README.md"

clust.files <- dir(clust.dir)
clust.names <- gsub(".txt","",clust.files)


#For each cluster create a vector with cluster names that are included in the cluster:
for (i in 1:length(clust.files)) {
  assign(clust.names[i],vector())
  clust.curr <- read.csv(paste0(clust.dir,clust.files[i]),sep=" ",header=F)
  names(clust.curr) <- c("Family","ID")
  
  for (j in 1:length(clust.files)) {
    if (i == j) next
    clust.check <- read.csv(paste0(clust.dir,clust.files[j]),sep=" ",header=F)
    names(clust.check) <- c("Family","ID")
    check.var <- setdiff(clust.check$ID,clust.curr$ID)
    if (length(check.var) == 0) assign(clust.names[i],c(get(clust.names[i]),clust.names[j]))
  }
}

output <- data.frame(matrix(NA,nrow=length(clust.files),ncol=3))
names(output) <- c("Parent","Cluster1","Cluster2")
output$Parent <- clust.names

for (i in 1:length(clust.files)) {
  t1 <- sort(gsub("cluster.","",get(clust.names[i])),decreasing=T)
  if (length(t1) == 0) {
    next
  } else if (length(t1) == 1) {
    output$Cluster1[which(output$Parent == clust.names[i])] = t1[1]
  } else {
    clust.curr <- read.csv(paste0(clust.dir,clust.files[i]),sep=" ",header=F)
    names(clust.curr) <- c("Family","ID")
    for (j in 1:(length(t1)-1)) {
      check.var <- T
      clust.check <- read.csv(paste0(clust.dir,paste0("cluster.",t1,".txt")[j]),sep=" ",header=F)
      names(clust.check) <- c("Family","ID")
      for (k in (j+1):length(t1)) {
        clust.check2 <- read.csv(paste0(clust.dir,paste0("cluster.",t1,".txt")[k]),sep=" ",header=F)
        names(clust.check2) <- c("Family","ID")
        if (identical(sort(union(clust.check$ID,clust.check2$ID)),sort(as.vector(clust.curr$ID)))) {
          output$Cluster1[which(output$Parent == clust.names[i])] = t1[j]
          output$Cluster2[which(output$Parent == clust.names[i])] = t1[k]
          check.var <- F
          break
        }
      }
      if (!check.var) break
    }
    
    if (check.var) {
      check.var2 <- "CHECK."
      for (j in 1:length(t1)) {
        check.var2 <- paste0(check.var2,t1[j],"-")
      }
      check.var2 <- substr(check.var2,1,nchar(check.var2)-1)
      output$Cluster1[which(output$Parent == clust.names[i])] = check.var2
    }
    
    
  }
}    

children <- c(output$Cluster1,output$Cluster2)
children <- children[!is.na(children)]
children.used <- children[-grep("CHECK",children)]
all.clusters <- gsub("cluster.","",clust.names)
children.not.used <- setdiff(all.clusters,children.used)

index <- row.names(output[grep("CHECK",output$Cluster1),])

for (i in index) {
  j <- as.numeric(i)
  t1 <- unlist(strsplit(gsub("CHECK.","",output$Cluster1[j]),"-"))
  
  if (length(intersect(t1,children.not.used)) == 1) {
    output$Cluster1[j] <- intersect(t1,children.not.used)[1]
    children.not.used <- children.not.used[children.not.used != intersect(t1,children.not.used)[1]]
  } else {
    output$Cluster1[j] <- paste(intersect(t1,children.not.used),collapse="-")
  } 
}

index <- row.names(output[grep("-",output$Cluster1),])

for (i in index) {
  j <- as.numeric(i)
  t1 <- unlist(strsplit(output$Cluster1[j],"-"))
  
  if (length(intersect(t1,children.not.used)) == 1) {
    output$Cluster1[j] <- intersect(t1,children.not.used)[1]
    children.not.used <- children.not.used[children.not.used != intersect(t1,children.not.used)[1]]
  } else {
    output$Cluster1[j] <- paste(intersect(t1,children.not.used),collapse="-")
  } 
}

output$Parent <- gsub("cluster.","",clust.names)

write.table(output,paste0(work.dir,"cluster_membership.txt"),quote=F,sep=",",na=".",row.names=F)
