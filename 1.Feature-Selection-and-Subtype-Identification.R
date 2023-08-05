library(CancerSubtypes)
library(NMF)

#Stage 1 Feature Selection
#For the microbial data matrix, the rows represent the features, and the columns represents the samples.
#1 Data normalization
microbial.data.rela <-apply(microbial.data,2, function(x){x/sum(x)})
#2 Feature selection
#2a Feature selection by Variance
microbial.data.select<-FSbyVar(microbial.data.rela,cut.type = "topk", value=100)
#2b Feature selection by Variant median absolute deviation
microbial.data.select<-FSbyMAD(microbial.data.rela,cut.type = "topk", value=100)

#Stage 2 Subtype Identification
#3 Clustering
#3a Consensus clustering
result<-ExecuteCC(d=microbial.data.select,clusterNum=3,maxK=10,
                 clusterAlg="pam",distance="pearson",title="microbial subtypes")
#3b Non-negative matrix factorization
result<-nmf(microbial.data.select, ranks=1:10, nrun=10, seed=0)
#3c Similarity network fusion
result<-ExecuteSNF(microbial.data.select, clusterNum=3, K=10, alpha=0.5, t=20)
#3d Integrative clustering (for multi-omics data)
multi.omics.dataset=list(microbiome=microbial.data.select,metabolism=metabolism.data)
result=ExecuteiCluster(datasets=multi.omics.dataset, k=3, lambda=list(0.44,0.33,0.28))   
#4 Validation
#4a Silhouette width
sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
#4b Proportion of ambiguous clustering
Kvec = 2:4
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = result$originalResult[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
optK = Kvec[which.min(PAC)]
#4c Consensus matrix
consensusMatrix<-result$originalResult[[i]]$consensusMatrix
#4d Cophenetic coefficient (for nmf)
coph <- result$measures$cophenetic

