## Create plots of the residual variance

setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/Data/ISO-NE Load Data")
dat = read.csv("LoadTest.csv")
dts = as.Date(dat[, 2])
dat = dat[, 3:26]
dat = 1e-3 * dat

dat2 = read.csv("PredTest.csv")
dat2 = dat2[, 3:26]
dat2 = 1e-3 * dat2

dat = dat - dat2


## Eigen space
library(ggplot2)
Sigma = cov(dat)
vals = eigen(Sigma, TRUE, only.values=TRUE)
exp_var = cumsum(vals$values)/sum(vals$values)

qplot(x=c(1:24), y=1e2*(1- exp_var)) + 
  geom_point() +
  geom_line() + 
  ylab("(%)") + 
  xlab("k") + 
  theme_bw(base_size=18)


### Buget set
stds = sort(diag(Sigma), decreasing=TRUE)
exp_var2 = cumsum(stds)/sum(stds)
qplot(x=c(1:24), y=1e2*(1- exp_var2)) + 
  geom_point() +
  geom_line() + 
  ylab("(%)") + 
  xlab("k") + 
  theme_bw(base_size=18)
