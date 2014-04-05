## Analyzing the Hindsight Costs with all reserves and Predictions

require(ggplot2)
require(reshape)

dat_hind = read.csv("allResPred_costs.csv")

dat_hind = cbind(dat_hind, totNom=dat_hind$Fixed + dat_hind$RealVariable, 
                            totHind=dat_hind$HindFixed + dat_hind$HindVariable)

ggplot(aes(x=totHind, y=totNom), data=dat_hind) + 
  geom_point() + 
  geom_abline(slope=1, intercept=0)

#melt it down to do two different histograms
