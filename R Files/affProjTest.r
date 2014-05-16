####
# Aff Proj and Clustering Results
###
library(ggplot2)
library(plyr)
library(reshape)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/Size_10/TestSet/")

dat.proj = read.delim("./affprojresults2.txt", header=TRUE, sep="\t")

#boxplots
ggplot(aes(x=NumEigs, y=NomTotCost-TotCost, group=NumEigs), 
       data=subset(dat.proj, Gap < .9)) + 
  geom_boxplot() +
  facet_grid(Set~., scales="free")

#point values
#User limit is an issue...
#separate plots for UCS and Ubudget
#Ucs
mean.savings = ddply(subset(dat.proj, Set=="UCS" & Gap<.9), ~NumEigs, function(df)mean(df$NomTotCost - df$TotCost))

ggplot(aes(x=NumEigs, y=(NomTotCost - TotCost)*1e-3), 
       data=subset(dat.proj, Gap < .9 & Set=="UCS")) + 
  geom_point(color="red") + 
  geom_line(aes(x=NumEigs, y=V1*1e-3), linetype="dashed",
            data=mean.savings) + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank()) + 
  xlab("Proj. Dim") + ylab("Savings ($M)")

mean.savings = ddply(subset(dat.proj, Set=="Budget" & Gap<.9), ~NumEigs, function(df)mean(df$NomTotCost - df$TotCost))
  
ggplot(aes(x=NumEigs, y=(NomTotCost - TotCost)*1e-3), 
         data=subset(dat.proj, Gap < .9 & Set=="Budget")) + 
  geom_point(aes(color=Status, size=100*Gap)) +
  geom_text(aes(label=round(100*Gap, 0), x=NumEigs, 
                   y=1e-3*(NomTotCost-TotCost), 
                  hjust=-.35, position="dodge"), 
              data=subset(dat.proj, Set=="Budget" &
                      Gap < .9 & Status=="UserLimit" & Gap> .2)) + 
  geom_line(aes(x=NumEigs, y=V1*1e-3), linetype="dashed",
            data=mean.savings) + 
  guides(color=FALSE, size=FALSE) + 
  theme_bw(base_size=16) + 
  xlab("Proj. Dim") + ylab("Savings ($M)") + 
  xlim(0, 10.5)
  

dat.by.method = cast(data=dat.proj, formula=Ind + NumEigs ~ Set, value="TotCost")
ggplot(aes(x=Budget * 1e-3, y=UCS * 1e-3, color=factor(Ind)), 
       data=dat.by.method) + geom_point()

#Solve Time Plots
ggplot(aes(x=NumEigs, y=SolveTime / 60, fill=Set, group=NumEigs), data=dat.proj)+
  geom_boxplot() + 
  facet_grid(.~Set) + 
  theme_bw(base_size=18)+
  theme(legend.position="none")+
  xlab("Proj. Dim.") +
  ylab("Sol. Time (min)")


### Same plot, just for UCS for presentation
ggplot(aes(x=NumEigs, y=SolveTime / 60, group=NumEigs), 
       data=subset(dat.proj, Set=="Budget"))+
  geom_boxplot(fill="light blue") + 
  theme_bw(base_size=18)+
  theme(legend.position="none")+
  xlab("Proj. Dim.") +
  ylab("Sol. Time (min)")



##################
#  Analysis fo the Cluster Stuff
####
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/All/TestSet/")
dat.clust = read.delim("./lpclustertestresults.txt", header=TRUE, sep="\t")
timediffs = ddply(dat.clust, 
                  .fun=function(df){c(mean(df$SolveTime), sd(df$SolveTime))}, 
                  .variables=~Mode)

ggplot(aes(x=Mode, y=SolveTime, fill=Mode), data=dat.clust) + 
  geom_boxplot()

ggplot(aes(x=Indx, y=SolveTime, color=Mode, group=Mode), 
       data=dat.clust) + 
  geom_point() +
  geom_line()

ggplot(aes(x=SolveTime, fill=Mode, group=Mode), 
       data=subset(dat.clust, Mode != "none")) + 
  geom_density(position="identity", alpha=.7) + 
  geom_histogram(position="identity", alpha=.7)

ggplot(aes(x=as.factor(Indx), y=SolveTime, group=Mode, fill=Mode), 
       data=dat.clust) + 
  geom_bar(stat="identity", position="dodge")
  geom_point() + geom_line()

