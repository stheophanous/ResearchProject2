rm(list = ls())

require('stats')
require('proxy')
require('pvclust')

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
P<-dist(BinaryMatrix, method="jaccard")

##Significant clustering evaluation
P2<-as.matrix(P)
hc2<-pvclust(P2, method.hclust="average", nboot = 200)
clsig <- pvpick(hc2, alpha=0.95, pv="au", type="geq", max.only=TRUE)
x <- as.character(clsig$edges)
mean(hc2$edges$au)
length(x)
mean(hc2$edges[x,1])
clsig$clusters