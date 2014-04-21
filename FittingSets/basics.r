###
#Basic analysis
##
# Basic plots for the paper

#read in the whole csv
dat.Loads = read.csv("HourlyLoadData.csv")
dat.Loads = dat.Loads[, -1]  #Kill fictional column
dat.Loads$Date = as.Date(dat.Loads$Date)

#Incorporate the hour information
#This is possibly not the smartest way
dat.Loads$DateHour = as.POSIXlt(dat.Loads$Date) + (dat.Loads$HOUR -1) * 3600


#Some simple plots
library(ggplot2)
library(reshape)

#1 Year of Data
ggplot(aes(x=DateHour, y=LOAD * 1e-3), 
       data=subset(dat.Loads, Date >= as.Date("2009-01-01") & Date <= as.Date("2010-01-01"))) + 
         geom_point() +
  xlab("") + ylab("(GWh)") +
  theme_bw(base_size=18)

#Just that Summer
isWkEnd = function(dt){ weekdays(dt) %in% c("Saturday", "Sunday")}
dat.Loads$isWkEnd = isWkEnd(dat.Loads$Date)
ggplot(aes(x=DateHour, y=LOAD * 1e-3), 
       data=subset(dat.Loads, Date >= as.Date("2009-07-01") & Date <= as.Date("2009-07-28"))) + 
  geom_point(aes(color=isWkEnd), size=3) + 
  geom_line() +
  xlab("") + ylab("(GWh)") +
  scale_color_discrete(labels=c("Weekday", "Weekend")) + 
  theme_bw(base_size=18) +
  theme(legend.title=element_blank(), 
        legend.position=c(.2, .8))

#Dependence on Temperature
ggplot(aes(x=DRY.BULB, y=LOAD * 1e-3), 
       data=subset(dat.Loads, HOUR == 13)) + 
  geom_point(aes(color=isWkEnd))+ 
  xlab("Temperaure") + 
  ylab("(GWh)") + 
  scale_color_discrete(labels=c("Weekday", "Weekend")) + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.2, .8))

#daily profiles
ggplot(aes(x=HOUR, y=LOAD * 1e-3, group=Date, color=Date), 
       data=subset(dat.Loads, Date >= as.Date("2009-07-01") & Date <= as.Date("2009-07-28"))
       ) +
  geom_point() + geom_line() + 
  ylab("(GWh)") + xlab("Hour") +
  theme_bw(base_size=18)

library(reshape)
hourlyLoads = cast(dat.Loads, Date~HOUR, value="LOAD")
hourlyDry = cast(dat.Loads, Date~HOUR, value="DRY.BULB")
hourlyDew = cast(dat.Loads, Date~HOUR, value="DEW.POINT")


#name the rows for ease
rownames(hourlyLoads) = hourlyLoads$Date
rownames(hourlyDry) = hourlyDry$Date
rownames(hourlyDew) = hourlyDew$Date

# Just give yourself one day of lags
# i.e. each laod observation has previous days Loads, Temps, it's own date and hour
newDat = dat.Loads[, c("Date", "HOUR", "DateHour", "LOAD")]
laggedDts = as.character(as.Date(newDat$DateHour - 3600 * 24))
newDat = cbind(newDat, hourlyLoads[laggedDts, 2:25])
names(newDat)[5:28] = paste("PLoad", 1:24, sep="_")
newDat = cbind(newDat, hourlyDry[laggedDts, 2:25])
names(newDat)[29:52] = paste("PDry", 1:24, sep="_")
newDat = cbind(newDat, hourlyDew[laggedDts, 2:25])
names(newDat)[53:76] = paste("PDew", 1:24, sep="_")

#Add a few more details that will be useful to the RF
#year, month, day of week, isWeekend
newDat = cbind(newDat, isWkEnd = isWkEnd(newDat$DateHour), 
                      Month = months(newDat$DateHour), 
                      Year = format(newDat$DateHour, "%Y"), 
                      DayOfWeek= weekdays(newDat$DateHour)) 


#write everything except first 24 entries, which are all NA
sum(is.na(newDat[-c(1:24), ]))
write.csv(newDat[-c(1:24), ], "laggedData.csv")

###It would be nice at this point to construct a periodogram of
# how many past hours you need

library(TSA)
t = spectrum(x=dat.Loads$LOAD)
t = spectrum(cos((1:1000)*20/2/pi))
