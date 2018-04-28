require('mclust')

#Adjusted Rand Index control
clustering <- read.table("/Users/steliostheophanous/Desktop/PerfClust.txt", sep = "\t", header = T)
perfclust <- clustering[,2]
control <- data.frame(perfclust)
indexlist <- list()
indexlist <- NULL
for(n in 1:1000) {
  control2 <-as.data.frame(lapply(control, function(cc) cc[sample(c(TRUE, NA), prob = c(0.41, 0.59), size = length(cc), replace = TRUE)]))
  control2[is.na(control2)] <- sample(1:55, 1, replace = TRUE)
  control3 <- as.vector(control2)
  newindex <- adjustedRandIndex(perfclust, control3$perfclust)
  indexlist[n] <- newindex
}
mean(indexlist)
sd(indexlist)
sd(indexlist)/sqrt(length(indexlist))

