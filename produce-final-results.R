#Input parameters - change them accordingly:

#work.dir is where we want our results to be saved
work.dir <- "/home/user_name/project/"
input.result <- "summary_yes_no.csv from the previous step"
cluster.member <- "cluster_membership.txt from one of the previous steps"

# ----------------------- START FUNCTIONS ----------------------------

#  check if they clust1 and clust2 are siblings
check.siblings <- function(clust1, clust2, cluster.member) {
  siblings <- FALSE  
  if (!length(cluster.member$Cluster2[which(cluster.member$Cluster1 == clust1)]) == 0) {
    ifelse(cluster.member$Cluster2[which(cluster.member$Cluster1 == clust1)] == clust2, siblings <- TRUE, siblings <- FALSE) 
  } else {
    ifelse(cluster.member$Cluster1[which(cluster.member$Cluster2 == clust1)] == clust2, siblings <- TRUE, siblings <- FALSE)
  }
  return(siblings)
}

who.is.parent <- function(sib1, cluster.member) {
  ifelse(sib1 %in% cluster.member$Cluster1,parent <- cluster.member$Parent[which(cluster.member$Cluster1 == sib1)],parent <- cluster.member$Parent[which(cluster.member$Cluster2 == sib1)])
  return(parent)
}

# function check.parent checks if any of the clusters that are hierarchically parents of it are in the results vector: returns vector of clusters that needs to be removed
check.parent <- function(clust, cluster.member, results) {
  max.cluster <- as.character(max(as.numeric(cluster.member$Parent)))
  to.remove <- vector()
  parent <- clust
  while (!(parent == max.cluster)) {
    parent <- who.is.parent(clust,cluster.member)
    clust <- parent
    ifelse(parent %in% results,to.remove <- c(to.remove,parent),next)
  }
  return(to.remove)
}

#TRUE if bettas are the same, FALSE if bettas are different:
z.score.check <- function(betta1, se1, betta2, se2) {
  z.score <- abs(betta1-betta2) / sqrt(se1^2 + se2^2)
  return( ifelse(z.score < 2, TRUE, FALSE) )
}

# ----------------------- END FUNCTIONS ----------------------------

output.result <- input.result

snps <- levels(droplevels(unique(output.result$SNP)))
for (snp in snps) {
  temp.df <- output.result[which(output.result$SNP == snp),]
  clusters <- temp.df$Clusters
  to.remove <- vector()
  
  if (nrow(temp.df) > 1) {
    combs <- combn(temp.df$Clusters,2)
    siblings <- vector()
    for (i in 1:ncol(combs)) {
      clust <- combs[,i]
      siblings.test <- check.siblings(clust[1],clust[2],cluster.member)
      if (siblings.test) {
        siblings <- c(siblings,c(clust[1],clust[2]))
      } else next
    }
    
    if (!length(siblings)==0) {
      for (i in 1:(length(siblings)/2)) {
        parent <- who.is.parent(siblings[i*2-1], cluster.member)
        z.score <- z.score.check(output.result$Betta[which(output.result$SNP == snp & output.result$Clusters == siblings[i*2-1])], output.result$SE[which(output.result$SNP == snp & output.result$Clusters == siblings[i*2-1])], output.result$Betta[which(output.result$SNP == snp & output.result$Clusters == siblings[i*2])], output.result$SE[which(output.result$SNP == snp & output.result$Clusters == siblings[i*2])])
        if (z.score) {          
          to.remove <- c(to.remove,c(siblings[i*2-1],siblings[i*2]))
        } else if ((parent %in% combs) & !z.score) {
          to.remove <- c(to.remove,parent)
        }
      }
    }
    
    ## remove clusters from clusters vector the ones that are in to.remove
    if (!(length(to.remove) == 0)) {
      to.remove <- sort(unique(to.remove))
      clusters <- clusters[-which(clusters %in% to.remove)]
    }
    #clusters.t <- clusters
    #to.remove.t <- to.remove
    for (i in clusters) {
      to.remove <- c(to.remove,check.parent(i,cluster.member,clusters))
    }
    
    
    if (!(length(to.remove) == 0)) {
      to.remove <- unique(to.remove)
      output.result <- output.result[-which(output.result$SNP == snp & output.result$Clusters %in% to.remove),]
    }
    
  }
}

#sort by p-values:
output <- output.result[order(output.result$P.value),]

write.table(output,paste0(work.dir,"summary_filtered.csv"),quote=F,sep=",",row.names=F)

