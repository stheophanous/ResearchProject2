rm(list = ls())

require('vegan')
require('apcluster')
require('proxy')

load("/Users/steliostheophanous/Documents/Bioinformatics/Project1/synergy/ActionPlan/drugs_ic50_matrix.RData")
drugs.m[is.na(drugs.m)] <- 0
d <- dist(t(drugs.m), method = "jaccard")
d2<-negDistMat(d,r=2)
ap.clust<-apcluster(d2,p=-0.0001, details = TRUE)

#Plot AP clustering and find clusters
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(ap.clust,fit$points)
text(x, y, labels = row.names(drugs.m), cex=0.5)
ap.clust@clusters
ap.clust
plot(ap.clust)
preferenceRange(d2, exact = TRUE)

#Find best preference value
preferenceRange(d2, exact = TRUE)
netsimlist <- list()
netsimlist <- NULL
clusterlist <- list()
clusterlist <- NULL
for(n in 1:100) {
  ap.clust <- apcluster(d2,p=-0.0001, details = TRUE)
  newnetsim <- ap.clust@netsim
  newcluster <- length(ap.clust@clusters)
  netsimlist[n] <- newnetsim
  clusterlist[n] <- newcluster
}
mean(netsimlist)
mean(clusterlist)
