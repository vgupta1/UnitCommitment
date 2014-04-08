###
#  Uses the lagged data to fit Random Forest Model
##

#read in the whole csv
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
library(Party)

#First try fitting one massive random forest
# fit <- randomForest(LOAD ~ ., data=dat.train[, -c(1,3)])

## compute out-of-sample predictions with error
# summer.test = subset(dat.test, Month %in% c("June", "July", "August"))
# summer.test = cbind(summer.test, 
#                     preds=predict(fit, summer.test))
# summer.test$resids = summer.test$LOAD - summer.test$preds

#create a plot RSq by hour
# t = split(summer.test, factor(summer.test$HOUR))
# sd_resids = sapply(t, function(df){sd(df$resids)})
# sd_loads = sapply(t, function(df){sd(df$LOAD)})
# mean_loads = sapply(t, function(df){mean(df$LOAD)})

# dat.sds = data.frame(HOUR=1:24, Residuals=sd_resids*1e-3, Load=1e-3*sd_loads)
# dat.sds = melt(dat.sds, id.vars="HOUR")
# ggplot(aes(x=HOUR, y=value, group=variable, color=variable), 
#        data=dat.sds) + 
#       geom_line() + geom_point() +
#       xlab("") + ylab("(GWh)") + 
#     theme_bw(base_size=18) + 
#   theme(legend.position=c(.2, .8), 
#         legend.title=element_blank() )

#Graph the Rsq separately
# ss_resids = sapply(t, function(df){mean(df$resids^2) * 1e-6})
# ss_loads = sapply(t, function(df){var(df$LOAD) * 1e-6})
# qplot(x=1:24, y=1-ss_resids/ss_loads) + geom_line() + 
#   xlab("") + ylab("R Sq.") + 
#   theme_bw(base_size=18) + 
#   geom_hline(aes(yintercept=overallR2), color="red", linetype="dashed")  


######
###That on the whole doesn't look great.  Try splitting by hour instead
######
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

#compute a data set of date fits, date trues
#Not the most efficient way
dat.test$pred = 0
for(iHr in 1:24)
{
  filt = dat.test$HOUR == iHr
  dat.test$pred[filt] = predict(fits[[iHr]], dat.test[filt, ])
}

##Do a plot to check that things look rational...
ggplot(aes(x=DateHour, y=(LOAD-pred)*1e-3), 
       data=subset(dat.test, Year==2008)) + 
  geom_point()

#This isn't the best thing ever... but it is what it is...
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
ggplot(aes(x=DateHour, y=(LOAD-pred)*1e-3), 
       data=subset(dat.validate, Year==2013)) + 
  geom_point()

hourlyload = cast(Date~HOUR, value="LOAD", data=dat.validate)
hourlypreds = cast(Date~HOUR, value="pred", data=dat.validate)

write.csv(hourlyload, "LoadValidate.csv")
write.csv(hourlypreds, "PredValidate.csv")

###
# Comparison for Dimitris, using the true temperatures
###
rm(list = ls())
dat = read.csv("HourlyLoadData.csv", stringsAsFactors=FALSE)
dat = dat[, -1] #Kill Wonky Column
hourlyDry = cast(dat, Date~HOUR, value="DRY.BULB")
hourlyDew = cast(dat, Date~HOUR, value="DEW.POINT")
rownames(hourlyDry) <- hourlyDry$Date
rownames(hourlyDew) <- hourlyDew$Date 

dts = as.character(dat$Date)
dat = cbind(dat, hourlyDry[dts, 2:25], hourlyDew[dts, 2:25])
names(dat)[6:29] <- paste("Dry", 1:24, sep="_")
names(dat)[30:53] <- paste("Dew", 1:24, sep="_")
dat = dat[, -c(4, 5)]

#flesh it out with a couple extra things
dat = cbind(dat, 
            isWkEnd = weekdays(as.Date(dat$Date)) %in% c("Saturday", "Sunday"), 
            Month = months(as.Date(dat$Date)), 
            Year = format(as.Date(dat$Date), "%Y"), 
            DayOfWeek= weekdays(as.Date(dat$Date)))
summer.train= subset(dat[1:32144, ], Month %in% c("June", "July", "August"))
summer.test= subset(dat[32145:64288, ], Month %in% c("June", "July", "August"))
rm(dat)

t=split(summer.train, summer.train$HOUR)
fits = lapply(t, function(df){  randomForest(LOAD~., data = df[, -c(1)] ) })
sapply(fits, function(f){mean(f$rsq)})

t = split(summer.test, summer.test$HOUR)
rsq_sum = function(df)
{ preds = predict(fits[[df$HOUR[1]]], df)
  c( mean((preds - df$LOAD)^2), var(df$LOAD)) }
t2 = sapply(t, rsq_sum)
rsq_resids= 1- t2[1, ] / t2[2, ]

qplot(x=1:24, y=rsq_resids) + geom_line() + 
  xlab("") + ylab("R Sq.") + 
  theme_bw(base_size=18) + 
  geom_hline(aes(yintercept= 1- sum(t2[1, ]) / sum(t2[2,])), color="red", linetype="dashed")  



###


