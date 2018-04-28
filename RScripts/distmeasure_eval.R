rm(list = ls())

require('stats')
require('proxy')
require('mclust')

Pred1 = read.csv("/Users/steliostheophanous/Documents/Bioinformatics/Project1/synergy/ActionPlan/200_lift_S.txt",sep="\t", header=FALSE )
Pred2 = read.csv("/Users/steliostheophanous/Documents/Bioinformatics/Project1/synergy/ActionPlan/200_support_S.txt",sep="\t", header=FALSE )

Pred1<-Pred1[,1:100] #Choose number of association mining rules
Pred2<-Pred2[,1:100] #Choose number of association mining rules
Pred<-cbind(Pred1, Pred2)

Groups = read.table("/Users/steliostheophanous/Documents/Bioinformatics/Project1/synergy/synergynew/groupsl.txt",sep="\t", header=FALSE )
Names = read.table("/Users/steliostheophanous/Documents/Bioinformatics/Project1/synergy/synergynew/names.txt",sep="\t", header=FALSE )
Names<-as.matrix(Names)
Names<-gsub("_IC50=Sensitive", "", Names)
Names<-paste(Names[,1],Groups[,1])

Features<-Pred[1,]
for (i in 2:nrow(Pred)){
  Features<-cbind(Features, Pred[i,])}
Features<-unique(Features)
Features<-as.matrix(Features)
Features<-as.character(Features)
Features<-unique(Features)

EmptyM<-matrix(, nrow = 251, ncol = length(Features))
colnames(EmptyM)<-Features
Drugs<-as.matrix(Names)
rownames(EmptyM)<-Drugs

L<-nrow(EmptyM)
for(b in 1:L)
{
  
  Ds<-Pred[b,]
  Ds<-Ds[,-1]
  Ds[Ds=="NA"] <- 0
  Ds <- Ds[, colSums(is.na(Ds)) == 0] 
  C<-ncol(Ds)
  C<-C
  A<-colnames(EmptyM)
  B<-c()
  
  for(i in 1:C)
  {
    
    X<-Ds[i]
    X<-as.matrix(X)
    X<-X[1,1]
    R<-pmatch(X,A)
    B<-c(B,R)
    EmptyM[b,R]<-1
    
  }
}
EmptyM[is.na(EmptyM)] <- 0
BinaryMatrix <- EmptyM

P<-dist(BinaryMatrix,method="jaccard")
hc<-hclust(P, method="average")
plot(hc, cex=0.3, main = "Hierarchical Clustering Dendrogram", cex.main=0.8)

##Mean Silhouette Width calculation
s<-NULL
silhlist <- NULL
silhlist <- list()
for(i in 2:200){
  k<-cutree(hclust(P,method="average"),i)
  s[i]<-summary(silhouette(k,P))$si.summary[3]
  silhlist[["newElement"]] <- s
}
lapply(silhlist, function(x) x[which.max(abs(x))])
silhlist2 <- as.data.frame(silhlist)
write.table(silhlist2, file = "/Users/steliostheophanous/Desktop/150jaccard2")
silhlist2[silhlist2 == 0] <- NA
mean(silhlist2[2:200,], na.rm = T)

#Adjusted Rand Index calculation
mycl <- cutree(hc, k = 55)
mycl2 <- as.data.frame(mycl)
clustering <- read.table("/Users/steliostheophanous/Desktop/PerfClust.txt", sep = "\t", header = T)
perfclust <- clustering[,2]
mycl3 <- mycl2[,1]
adjustedRandIndex(mycl3, perfclust)

#Purity Calculation
puritylist <- list()
puritylist <- NULL
for(k in 1:100) {
  randomclust <- list()
  randomclust <- NULL
  for (n in 1:251) {
    newnum <- sample(1:55, replace = TRUE)
    randomclust[n] <- newnum
  }
  ClusterPurity <- function(clusters, classes) {
    sum(apply(table(classes, clusters), 2, max)) / length(clusters)
  }
  newpurity <- ClusterPurity(mycl3,randomclust)
  puritylist[k] <- newpurity
}
mean(puritylist)





