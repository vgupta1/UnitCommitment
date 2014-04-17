library(boot)

setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/Data/ISO-NE Load Data")

dat.pred = read.csv("PredTest.csv")
dat.pred = dat.pred[, 3:26]
dat.pred = 1e-3 * dat.pred

dat.true = read.csv("LoadTest.csv")
dat.true = dat.true[, 3:26]
dat.true = 1e-3 * dat.true

resids = dat.true - dat.pred
gam1 = function(data, indices){ norm(colMeans(data[indices, ]) - colMeans(data), "2")}
gam2 = function(data, indices){ norm( cov(data[indices, ]) - cov(data), "F") }

boot.out = boot(resids, gam1, 5000)
quantile(boot.out$t, c(.95, .9, .85, .8, .75))

boot.out2 = boot(resids, gam2, 5000)
quantile(boot.out2$t, c(.95, .9, .85, .8, .75))
