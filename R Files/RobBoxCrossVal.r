#Cross Val Parameter Rob, Box

setwd("Documents/Research/UnitCommittment/UnitCommitment/results/All/")
dat = read.delim("cv_boxrob.txt", header=FALSE)
names(dat) <- c("ratio", "Mean", "Std")

#round the Gammas for readability
dat$ratio = round(dat$ratio, 3)

library(ggplot2)
ggplot(aes(x=100*(1 - 2*(1-ratio)), y=-Mean), data=dat) +
  geom_point() + 
  geom_errorbar(aes(ymin=-Mean + Std, ymax=-Mean-Std))+
  ylab("Improv.over Nominal ($K)") + xlab("% Pts Covered") +
  theme_bw(base_size=18)
