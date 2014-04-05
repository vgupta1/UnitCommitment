#Does some analysis of the different generator types and capabilities

library(gdata)
library(ggplot2)
library(reshape)

setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment")
dat = read.csv(file="resource_details.csv")


#look at capacity sizes by fuel_type
ggplot(aes(x=Capacity, color = FuelType, fill=FuelType, group=FuelType), data=dat) + 
  geom_histogram() + 
  facet_grid(FuelType~.)

ggplot(aes(x=FuelType, color=FuelType, y=Capacity), data=dat) + 
  geom_boxplot()

ggplot(aes(x=StartCost + 1, y=Capacity + 1, color=FuelType), data=dat) + geom_point() + 
  scale_y_log10() + scale_x_log10()


#split-apply-combine strategy to find the avg capacity by fuel_type and similar
require(plyr)
#summarize is a function which allows you to define a dataframe by various functions
ddply(dat, "FuelType", .fun=summarize, 
      Capacity=mean(Capacity), 
      StartCost=mean(StartCost), 
      MinUp = mean(MinUp), 
      MinDown = mean(MinDown))

dat_offer = read.csv("offerCurves.csv")
ggplot(aes(x=Amount,  y=Price, color=FuelType, group=Name), 
       data=subset(dat_offer, FuelType %in% c("Diesel"))) + 
  geom_line() + geom_point() + 
  facet_grid(FuelType~.)

ggplot(aes(x=Amount,  y=Price, color=FuelType, group=Name), 
       data=dat_offer) + 
  geom_line() + geom_point() + 
  facet_grid(FuelType~.)

