## 
# Analyzing Jacknife Results
##

library(ggplot2)
library(plyr)
library(reshape)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/Size_10")
dat.ucs.split = read.delim("jacknife_ucs_results.txt", header=TRUE, sep="\t")

table(dat.ucs.split$Epsilon, dat.ucs.split$Gamma2, dat.ucs$Gamma1)
dat.ucs = ddply(dat.ucs.split, .variables=~ Epsilon + Gamma1 + Gamma2, 
      .fun=function(df){c(mean(df$NomVals-df$TotCost), sd(df$TotCost-df$NomVals)/sqrt(10), length(df$TotCost))})
names(dat.ucs)[4:6] = c("Mean", "SD", "Num")

dat.ucs$Gamma1 = round(dat.ucs$Gamma1, 2)
dat.ucs$Gamma2 = round(dat.ucs$Gamma2, 2)

#build the graph up by bits
ggplot(aes(x=Epsilon, y=Mean, 
           color=factor(Gamma1), factor(group=Gamma1), 
           position="dodge"), 
       data=subset(dat.ucs, Gamma2==.23)) + 
  geom_point() + geom_line(linetype="dashed") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD))
  theme(legend.title=element_blank())

#Seemingly no difference between parameters!!!!
#marginally better than the rest:
#Epsilon Gamma1 Gamma2     Mean       SD Num
#32     0.9   0.14   0.27 8855.703 3545.829  10



##### Repeat this exercise for the budget
dat.bud.split = read.delim("jacknife_budget_results.txt", header=TRUE, sep="\t")
dat.bud = ddply(dat.bud.split, .variables=~ Gamma1 + GammaBound, 
                .fun=function(df){c(mean(df$NomVals-df$TotCost), sd(df$TotCost-df$NomVals)/sqrt(10), length(df$TotCost))})
names(dat.bud)[3:5] = c("Mean", "SD", "Num")

#drop the last one for now
dat.bud = dat.bud[-c(10), ]
dat.bud$Gamma1 = round(dat.bud$Gamma1, 2)

dat.ucs$Gamma1 = round(dat.ucs$Gamma1, 2)
dat.ucs$Gamma2 = round(dat.ucs$Gamma2, 2)

#looks like these actually make a difference
pos = position_dodge(1)
ggplot(aes(x=Gamma1, y=Mean, 
           color=factor(GammaBound), factor(group=GammaBound)), 
       data=dat.bud) + 
  geom_point(position=pos) + geom_line(linetype="dashed") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=pos) +
theme(legend.title=element_blank())




