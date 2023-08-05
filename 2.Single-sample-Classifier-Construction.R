library(pROC)

#Stage 3 Single-sample Classifier Construction
#5 Marker identification
#5a AUC
subtype.info$S1<-ifelse(subtype.info$subtype=='S1','S1','others')
subtype.info$S2<-ifelse(subtype.info$subtype=='S2','S2','others')
subtype.info$S3<-ifelse(subtype.info$subtype=='S3','S3','others')
auc4S1<-c();auc4S2<-c();auc4S3<-c()
for (i in rownames(microbial.data.rela)){
  s1<-auc(subtype.info[,'S1'],t(microbial.data.rela)[,i])
  s2<-auc(subtype.info[,'S2'],t(microbial.data.rela)[,i])
  s3<-auc(subtype.info[,'S3'],t(microbial.data.rela)[,i])
  auc4S1<-c(auc4S1,s1);auc4S2<-c(auc4S2,s2);auc4S3<-c(auc4S3,s3)
}
auc.res<-data.frame(S1.auc=auc4S1,S2.auc=auc4S2,S3.auc=auc4S3)
rownames(auc.res)<-rownames(microbial.data.rela)
cutoff=0.7
marker4S1<-rownames(auc.res[auc.res$S1.auc>cutoff,])
marker4S2<-rownames(auc.res[auc.res$S2.auc>cutoff,])
marker4S3<-rownames(auc.res[auc.res$S3.auc>cutoff,])
marker<-unique(c(marker4S1,marker4S2,marker4S3))
centroids.table <- matrix(NA, nrow=length(marker), ncol=3,
                          dimnames=list(marker,c('S1','S2','S3')))
marker.table<-as.data.frame(t(microbial.data.rela[marker,]))
marker.table$group<-subtype.info$subtype
for (i in colnames(marker.table)[1:(ncol(marker.table)-1)]){
  s1<-mean(marker.table[marker.table$group=='S1',i])
  s2<-mean(marker.table[marker.table$group=='S2',i])
  s3<-mean(marker.table[marker.table$group=='S3',i])
  centroids.table[i,]<-c(S1,S2,S3)
}  
#5b RandomForest
#5c Wilcoxon rank sum test
#5d Kruskal-Wallis test

#6 Classifier construction
classifyMicrobial <- function(x, minCor = .2,centroids){
  
  classes<-colnames(centroids)
  
  gkeep <- intersect(rownames(centroids), rownames(x))
  if (length(gkeep) == 0) stop("Empty intersection between profiled features and the features used for classification.\n Make sure that gene names correspond to the type of identifiers specified by the gene_id argument")
  if (length(gkeep) < 0.6 * nrow(centroids)) warning("Input feature profiles include less than 60% of the features used for classification. Results may not be relevant")
  
  cor.dat <- as.data.frame(cor(x[gkeep, ], 
                               centroids[match(gkeep, rownames(centroids)), classes],
                               use = "complete.obs"), row.names = colnames(x))
  ##########
  cor.dat[is.na(cor.dat)]<-0
  # Best correlated centroid
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y){classes[which.max(y)]})
  cor.dat$corToNearest <- apply(cor.dat[, classes], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp){
    cor.test(x[gkeep, smp], centroids[match(gkeep, rownames(centroids)), as.character(cor.dat[smp, "nearestCentroid"])])$p.value
  })
  
  # Separation level metrics
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - cor.dat[, classes], 1, function(x){sort(x)[2]})
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, classes], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed
  cor.dat[is.na(cor.dat$cor_pval),'nearestCentroid']<-NA
  
  cor.dat$NMIBC_Class <- cor.dat$nearestCentroid
  
  # Set to NA if best correlation < minCor
  try(cor.dat[which(cor.dat$corToNearest < minCor), "NMIBC_Class"] <-  NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <-  NA)
  
  #cor.dat <- cor.dat[, c("NMIBC_class" , "cor_pval", "separationLevel", classes)]
  return(cor.dat)
}


