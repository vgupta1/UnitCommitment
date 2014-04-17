##### 
# Attempts some cluster analaysis on the forecasts
#####

setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/Data/ISO-NE Load Data")
dat = read.csv("LoadTest.csv")
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
         main="Loads")


