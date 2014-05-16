## 
# Analyzing Jacknife Results
##

library(ggplot2)
library(plyr)
library(reshape)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/Size_10/TestSet/")
dat.ucs.split = read.delim("jacknife_ucs_results.txt", header=TRUE, sep="\t")

table(dat.ucs.split$Epsilon, dat.ucs.split$Gamma2, dat.ucs.split$Gamma1)
dat.ucs = ddply(dat.ucs.split, .variables=~ Epsilon + Gamma1 + Gamma2, 
      .fun=function(df){c(mean(df$NomVals-df$TotCost), sd(df$TotCost-df$NomVals)/sqrt(10), length(df$TotCost))})
names(dat.ucs)[4:6] = c("Mean", "SD", "Num")

dat.ucs$Gamma1 = round(dat.ucs$Gamma1, 2)
dat.ucs$Gamma2 = round(dat.ucs$Gamma2, 2)

#build the graph up by bits
dat.ucs[ ,4:5] = dat.ucs[ ,4:5] * 1e-3

pos = position_dodge(.05)
ggplot(aes(x=Epsilon, y=Mean, 
           color=factor(Gamma1), group=factor(Gamma1)),
       position=pos,
       data=subset(dat.ucs, Gamma2 %in% c(.17, .27))) + 
  geom_point(position=pos) + 
  geom_line(linetype="dashed", position=pos) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=pos) +
  facet_grid(.~Gamma2) + 
  theme_bw(base_size=16) +
#  theme(legend.title=element_blank()) + 
  scale_color_discrete(name="Gamma1") +
  ylab("Savings ($M)") + 
  theme(legend.position="top")

#Seemingly no difference between parameters!!!!
#marginally better than the rest:
#Epsilon Gamma1 Gamma2     Mean       SD Num
#32     0.9   0.14   0.27 8855.703 3545.829  10



##### Repeat this exercise for the budget
dat.bud.split = read.delim("jacknife_budget_results.txt", header=TRUE, sep="\t")
dat.bud = ddply(dat.bud.split, .variables=~ Gamma1 + GammaBound, 
                .fun=function(df){c(mean(df$TotCost-df$NomVals) * 1e-3, sd(df$TotCost-df$NomVals)/sqrt(10) * 1e-3, length(df$TotCost))})
names(dat.bud)[3:5] = c("Mean", "SD", "Num")

dat.bud$Gamma1 = round(dat.bud$Gamma1, 2)
dat.bud.split$Gamma1 = round(dat.bud.split$Gamma1, 2)

#looks like these actually make a difference
pos = position_dodge(1)
ggplot(aes(x=Gamma1/sqrt(24), y=Mean, 
           color=factor(GammaBound),
           group=factor(GammaBound)),
       position=pos,
       data=dat.bud) + 
  geom_point(position=pos) + 
  geom_line(linetype="dashed", position=pos) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=pos) +
  theme_bw(base_size=16) + 
  scale_color_discrete(name="Gamma2" ) + 
  ylab("Savings ($M)") + theme(legend.position="top") + 
  xlab("Gamma1")
  

ggplot(aes(x=GammaBound, y=Mean, 
           color=factor(Gamma1),
           group=factor(Gamma1)),
       position=pos,
       data=dat.bud) + 
  geom_point(position=pos) + 
  geom_line(linetype="dashed", position=pos) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=pos) +
  theme(legend.title=element_blank())

ggplot(aes(x=factor(Gamma1), y=TotCost - NomVals, 
           group=factor(Gamma1), 
           ), 
       data=dat.bud.split) + 
  geom_boxplot() +
  facet_grid(.~GammaBound)

ggplot(aes(x=factor(round(Gamma1/sqrt(24), 2)), y=TotCost, 
           color=factor(GammaBound), 
           group=factor(GammaBound)), 
       data=dat.bud.split) + 
  geom_point() + 
  geom_smooth() +
  theme(legend.position="top")
1.282 * sqrt(24)


