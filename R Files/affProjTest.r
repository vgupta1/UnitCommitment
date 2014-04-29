####
# Aff Proj and Clustering Results
###
library(ggplot2)
library(plyr)
library(reshape)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/Size_10/")

dat.proj = read.delim("./affprojresults.txt", header=TRUE, sep="\t")

dat.by.method = cast(data=dat.proj, formula=Ind + NumEigs ~ Set, value="TotCost")
ggplot(aes(x=Budget * 1e-3, y=UCS * 1e-3, color=factor(Indx)), 
       data=dat.by.method) + geom_point()

### UserLimit seem to kick in faster than expected and destroy results
ggplot(aes(x=NumEigs, y=1e-3*(TotCost-nomVal), color=Status, group=Set), 
       data=dat.proj) +
  geom_point() + 
  facet_grid(Set~., scales="free")


#Solve Time Plots
ggplot(aes(x=NumEigs, y=SolveTime, fill=Set, group=NumEigs), data=dat.proj)+
  geom_boxplot() + 
  facet_grid(.~Set)

##################
#  Analysis fo the Cluster Stuff
####
dat.clust = read.delim("./clusterresults.txt", header=FALSE, sep="\t")
names(dat.clust)<- c("mode", "indx", "status", "TotCost", "NomCost", "SolveTime")
timediffs = ddply(dat.clust, 
                  .fun=function(df){c(mean(df$SolveTime), sd(df$SolveTime))}, 
                  .variables=~mode)
ggplot(aes(x=mode, y=SolveTime, fill=mode), data=dat.clust) + 
  geom_boxplot()

ggplot(aes(x=indx, y=SolveTime, color=mode, group=mode), 
       data=subset(dat.clust, mode != "none")) + 
  geom_point() +
  geom_line()

ggplot(aes(x=SolveTime, fill=mode, group=mode), 
       data=subset(dat.clust, mode != "none")) + 
  geom_density(position="identity", alpha=.7) + 
  geom_histogram(position="identity", alpha=.7)
