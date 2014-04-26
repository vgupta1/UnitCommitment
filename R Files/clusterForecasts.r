##### 
# Attempts some cluster analaysis on the forecasts
#####

setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/Data/ISO-NE Load Data")
dat = read.csv("LoadTest.csv")
dts = as.Date(dat[, 2])
dat = dat[, 3:26]
dat = 1e-3 * dat

dat2 = read.csv("PredTest.csv")
dat2 = dat2[, 3:26]
dat2 = 1e-3 * dat2

dat = dat - dat2
rm(dat2)

### using k means
wss <- (nrow(dat)-1)*sum(apply(dat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(dat, 
                                     centers=i)$withinss)

library(ggplot2)
qplot(x=1:10, y=wss) + geom_point() + geom_line()

qplot(x=1:10, y=wss) + 
  geom_point() + geom_line() + 
  xlab("Clusters") + 
  ylab("Within Group SoS") + 
  geom_text(aes(label=1:10, x=1:10, y=wss),  vjust=-1, hjust=0) +
  ylim(0, 26000) + xlim(1, 10.5) +
  theme_bw(base_size=18)

#3 maybe 4 clusters?
fit3 <- kmeans(dat, 3)
fit3$means = aggregate(dat,by=list(fit3$cluster),FUN=mean)

fit4 <- kmeans(dat, 4)
fit4$means = aggregate(dat,by=list(fit4$cluster),FUN=mean)

# Cluster Plot against 1st 2 principal components
# vary parameters for most readable graph
library(cluster) 
clusplot(dat, fit4$cluster, color=TRUE, shade=TRUE, 
         labels=1, lines=0, col.txt="black",
         main="", 
         sub="", xlab="", ylab="")


#identify the cluster labels for following indices
INDXSET = c(116, 51, 239, 118, 73, 59, 218, 220, 99, 227)  #obtained by 10 k means
clusters = fit4$cluster[INDXSET]

fit4$cluster[c(72, 88, 221, 160)]
#72 -> 3
#88 -> 2
#221 -> 1
#160 -> 4

library(plyr)
clustermap = revalue(factor(fit4$cluster), c("3"="72", "2"="88", "1"="221", "4"="160"))
write.csv(file="testClusterMap.csv", clustermap)


## create the clusterMap for the vlaidation set
dat.validate = read.csv("PredValidate.csv")
dts.validate = as.Date(dat.validate[, 2])
dat.validate = dat.validate[, 3:26]
dat.validate = 1e-3 * dat.validate

closest.cluster <- function(x) {
  cluster.dist <- apply(fit4$centers, 1, function(y) sqrt(sum((x-y)^2)))
  return(which.min(cluster.dist)[1])
}
clusters2 <- apply(dat.validate, 1, closest.cluster)

#now need to remap them
#72 -> 3
#88 -> 2
#221 -> 1
#160 -> 4
library(plyr)
clustermap2 = revalue(factor(clusters2), c("3"="72", "2"="88", "1"="221", "4"="160"))


write.csv(file="../../results/All/ValidateSet/validateClusterMap.csv", clustermap2)
clustermap2[1:5]

