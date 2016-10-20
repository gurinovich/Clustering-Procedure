# Clustering-Procedure

Clustering procedure to identify significant SNPs associated with a binary phenotype while accounting for population structure

1. Principal component analysis: calculate-PCA.R
  * Input:
    1. file with samples: tab-delimited; 2 columns: Family-ID & ID
    2. genetic data in PLINK bed/bim/fam format
    3. IBD and Kinship estimations, e.g. from [MERLIN][1]
2. Hierarchical clustering of the samples using PCs: clustering.R
  * Input:
    1. file (.csv) with at least the following columns: Family, Id, Principal components
    2. number of principal components used for hierarchical clustering of samples; suggested: 6
    3. size of the smallest allowed cluster in hierarchical clustering; suggested: 100
3. For each cluster derived from the previous step, extract the information who are a parent and children of the cluster: cluster-membership.R
   * Input: output from the previous step
   * Output: cluster_membership.txt -> file with 3 columns: Parent, Cluster1, Cluster2
4. Run respective linear and binary logistic regression models on the genetics data: create-summary.R
    * Input:
      1. .txt files for each cluster with columns (no names): "ID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30" 
      2. .csv file with columns: phenotype, ID, Family, Sex, SNP-dosages
      3. .txt file with the list of SNPs to investigate
    * Output: work.dir/full-summary/ -> .txt files for each SNP with columns: Clusters,PC1.b,PC1.p,PC2.b,PC2.p,PC3.b,PC3.p,PC4.b,PC4.p,PC5.b,PC5.p,PC6.b,PC6.p,GLM.b,GLM.p,GLM.SE
5. Assign yes/no labels to the model parameters that pass respective thresholds: convert-summary-yes-no.R
     * Input:
      1. output from the previous step (work.dir/full-summary/)
      2. p-values for the significances of associations between PCs&SNPs and phenotype&SNPs; suggested: 0.05/6 and 0.05 respectively
     * Output: work.dir/full-summary-yes-no/ -> .txt files for each SNP with columns: Clusters,PC1.PC6.signif,SNP.signif: yes/no independence/significance association
6. Create the result file with the following columns: SNP, Clusters, Betta, SE, P.value: create-summary-yes-no.R
     * Input: output from the previous steps: work.dir/full-summary-yes-no/ and work.dir/full-summary/
     * Output: summary_yes_no.csv -> a file with all the SNPs - clusters combinations that are independent of PCs & significantly associated with the phenotype.
7. Filtering out to produce final results: produce-final-results.R
     * Input: output from previous steps: summary_yes_no.csv and cluster_membership.txt
     * Output: summary_filtered.csv
      
      
      
[1]: http://csg.sph.umich.edu/abecasis/Merlin/tour/ibd.html "Title"
