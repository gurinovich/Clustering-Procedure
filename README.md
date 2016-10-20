# Clustering-Procedure

Clustering procedure to identify significant SNPs associated with a binary phenotype while accounting for population structure

1. Principal component analysis: calculate-PCA.R
  * Required input files:
    1. file with samples: tab-delimited; 2 columns: Family-ID & ID
    2. genetic data in PLINK bed/bim/fam format
    3. IBD and Kinship estimations, e.g. from [MERLIN][1]
2. Hierarchical clustering of the samples using PCs: clustering.R
  * Required:
    1. file (.csv) with at least the following columns: Family, Id, Principal components
    2. number of principal components used for hierarchical clustering of samples; suggested: 6
    3. size of the smallest allowed cluster in hierarchical clustering; suggested: 100


[1]: http://csg.sph.umich.edu/abecasis/Merlin/tour/ibd.html "Title"
