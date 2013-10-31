## Build some advanced regression models on teh load data set

library(gdata)
library(ggplot2)
library(reshape)

setwd("/Users/VGupta/Documents/Research/UnitCommittment/ISO-NE/")
load_hist.orig = read.csv('load_hist_orig.csv')
load_hist.orig$dateTime <- as.POSIXct(load_hist.orig$dateTime)
load_hist.orig <- load_hist.orig[order(load_hist.orig$Date, load_hist.orig$Hour) , ]

#Limit Peak Days in June, July, August
load.summer <- subset(load_hist.orig, months(dateTime) %in% c("June", "July", "August") &
                        !isWkEnd )

#save the last 186 days for validation.  First 500 for training/testing:
# 2010-07-19

#drop 2013 for testing
load.validate <- subset(load.summer, as.Date(Date) > "2010-07-09"  )
load.train <- subset( load.summer, as.Date(Date) <= "2010-07-09")
rm(load_hist.orig, load.summer)

#reshape the training data to a form suitable for CART
#for now just use the hourly temperatures for the day.  
#shoudl be replaced with forecasts
#should add previous load profile and indicator for day of week.
load.train.hour = cast(Date~Hour, data = load.train, value="DEMAND")
drybulb.train.hour = cast(Date~Hour, data=load.train, value="DryBulb")
dewpnt.train.hour = cast(Date~Hour, data = load.train, value="DewPnt")
names(load.train.hour) = paste("Load", names(load.train.hour), sep="")
names(drybulb.train.hour) = paste("Dry", names(drybulb.train.hour), sep="")
names(dewpnt.train.hour) = paste("Dew", names(dewpnt.train.hour), sep="")
load.train.hour = cbind(load.train.hour, 
                        drybulb.train.hour[, 2:25], 
                        dewpnt.train.hour[, 2:25])

#just try to fit hour 13 first
library(rpart)
trial_tree = rpart(Load ~ ., method="anova", data=cbind(Load = load.train.hour$Load13, 
                              drybulb.train.hour[, 2:25], 
                              dewpnt.train.hour[, 2:25]), 
                    cp=0, xval=10)
fit_tree = prune(trial_tree, cp = .003)

