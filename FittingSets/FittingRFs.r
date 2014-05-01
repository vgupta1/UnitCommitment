###
#  Uses Lagged data to fit Random Forest Model
##

#read in the whole csv
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/Data/ISO-NE Load Data/")
dat = read.csv("laggedData.csv")
dat = dat[, -1 ]  #kill wonky column
dat$DateHour = as.POSIXlt(dat$Date) + (dat$HOUR-1) * 3600  #Fix the datehour
dat$Date = as.character(dat$Date)

#Filter down to only looking at summers, divide roughly in 3 (respecting days)
dat = subset(dat, Month %in% c("June", "July", "August") & !isWkEnd)
dat.train = dat[1:5784, ]
dat.test = dat[5785:11568, ]
dat.validate = dat[11569:17352, ]
rm(dat)

library(randomForest)
#library(Party)
library(ggplot2)  #for creating graphs

fitRF = function(df){ randomForest(LOAD ~ ., data=df[,-c(1,3)]) }
t = split(dat.train, dat.train$HOUR)  
fits = lapply(t ,fitRF)
sapply(fits, function(f){mean(f$rsq)})

#compute the out of sample fits and residuals
t = split(dat.test, dat.test$HOUR)
predRSE = function(df)
{ preds = predict(fits[[df$HOUR[1]]], df)
  c(mean((preds-df$LOAD)^2), var(df$LOAD)) }
t2 = sapply(t, predRSE)
rsq_resids = 1- t2[1, ]/t2[2, ]
overall = 1- sum(t2[1, ]) / sum(t2[2, ])

qplot(x=1:24, y=rsq_resids) + geom_line() + 
  xlab("") + ylab("R Sq.") + 
  theme_bw(base_size=18) + 
  geom_hline(aes(yintercept=overall), color="red", linetype="dashed")  

#compute a data set of fits, trues
#Not the most efficient way
dat.test$pred = 0
for(iHr in 1:24)
{
  filt = dat.test$HOUR == iHr
  dat.test$pred[filt] = predict(fits[[iHr]], dat.test[filt, ])
}

##Do a plot to check that things look rational...
#This isn't the best thing ever... but it is what it is...
ggplot(aes(x=DateHour, y=LOAD*1e-3), 
       data=subset(dat.test, as.Date(dat.test$Date) >= as.Date("2008-07-01") 
                            & as.Date(dat.test$Date) <= as.Date("2008-07-7"))) + 
  geom_point() + 
  geom_point(aes(y=pred * 1e-3), color="blue")

hourlyload = cast(Date~HOUR, value="LOAD", data=dat.test)
hourlypreds = cast(Date~HOUR, value="pred", data=dat.test)

write.csv(hourlyload, "LoadTest.csv")
write.csv(hourlypreds, "PredTest.csv")

## Write same csv for the validation set (oh my!)
dat.validate$pred = 0
for(iHr in 1:24)
{
  filt = dat.validate$HOUR == iHr
  dat.validate$pred[filt] = predict(fits[[iHr]], dat.validate[filt, ])
}

hourlyload = cast(Date~HOUR, value="LOAD", data=dat.validate)
hourlypreds = cast(Date~HOUR, value="pred", data=dat.validate)

write.csv(hourlyload, "LoadValidate.csv")
write.csv(hourlypreds, "PredValidate.csv")