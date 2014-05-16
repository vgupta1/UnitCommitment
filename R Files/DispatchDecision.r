library(ggplot2)
library(reshape)
library(plyr)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/")

dat = read.delim("DispatchScheds", "\t", header=TRUE)

ggplot(aes(x=FuelType, y=HR1, fill=Set, group=Set),
       data=dat) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw(base_size=16) +
  theme(legend.title=element_blank(),
        legend.position=c(.2, .8)) + 
  xlab("") + ylab("GWh") + ggtitle("12:00 AM")

ggplot(aes(x=FuelType, y=HR13, fill=Set, group=Set),
       data=dat) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw(base_size=16) +
  theme(legend.title=element_blank(), 
        legend.position=c(.2, .8)) + 
  xlab("") + ylab("GWh") + ggtitle("12:00 PM")
