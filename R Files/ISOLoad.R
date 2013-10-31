###
# Exploring the Load Data set
###
library(gdata)
library(ggplot2)

setwd("/Users/VGupta/Documents/Research/UnitCommittment/ISO-NE/")
load_hist = read.xls("smd_hourly_2013b.xls", sheet=2)
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2012.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2011.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2010.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2009.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2008.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2007.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2006.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2005.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2004.xls", sheet=2)))
load_hist = data.frame(rbind(load_hist, read.xls("smd_hourly_2003.xls", sheet=2)))

#add a column of true dates
t = with(load_hist, paste(Date, Hour), sep=" ")
t = paste(t, ":00", sep = "")
load_hist = data.frame(load_hist, dateTime = as.POSIXct(t))

q10 <- function(v){ quantile(v, .05)}
q90 <- function(v){ quantile(v, .95)}

#label the weekday/weekend
#really should be done by holiday calendar
load_hist = data.frame(load_hist, 
                       isWkEnd=weekdays(as.POSIXct(load_hist$Date)) %in% c("Saturday", "Sunday"))

#identify the minimum and maximum loads
OFF_MAX = with(subset(load_hist, isWkEnd), max(DEMAND))
PK_MAX = with(subset(load_hist, !isWkEnd), max(DEMAND))

OFF_MIN = with(subset(load_hist, DEMAND > 0 & isWkEnd), min(DEMAND))
PK_MIN = with(subset(load_hist, DEMAND > 0 & !isWkEnd), min(DEMAND))

##limiting to just summer months
##Weekday = [8182, 27622]

sum(load_hist$DEMAND >= 27600)/90600
sum(0 < load_hist$DEMAND & load_hist$DEMAND<= 7900)/90600

#do a big plot of everything by hour
#mark the daylight savings
t = rep(FALSE, 8784)

#VG Find a more automatic way to do this
t[7394] = TRUE
t[1682] = TRUE
load_hist = data.frame(load_hist, dayLightSavings=t)

g <- ggplot(aes(x=dateTime, y=DEMAND), 
            data= load_hist.orig) +
  geom_point() + theme_bw() + theme(legend.position="None") + 
  xlab("") + ylab("ISO-NE System Load (MWh)")


#limit to a single month and zone
# smooth by hour
ggplot(aes(x=dateTime, y =DEMAND), 
       data=subset(load_hist.orig, 
                   months(dateTime)=="June")) +
  theme_bw() + 
  geom_point() + geom_line() + 
  stat_summary(fun.y="mean", geom="point") +
  stat_summary(fun.y="mean", fun.ymax=q90, fun.ymin=q10, 
               geom="errorbar", position="identity") +
  stat_smooth(aes(group=1))

#compute means by hour.  

#we will be a bit lazy and do this only for June
#VG Remove this
load_hist.orig <- load_hist
load_hist <- subset(load_hist.orig, months(dateTime) %in% c("June", "July", "August"))

#drop june 2013 for the testing
load_hist <- load_hist[722:22801, ]


#add a column to the data frame which is the mean
mean_df <- function(df){ mean(df$DEMAND)}
load_means = daply(load_hist,
                   .variables=~isWkEnd + Hour, .fun=mean_df)
mean_fun <- function(df){ load_means[as.character(df$isWkEnd), df$Hour]}
load_hist = adply(load_hist, .fun=mean_fun, .margins=1)
names(load_hist)[8] <- "meanLoad"

#do a plot by day + hour within june group by year
dayHourOfMonth <- function(dt)
{ 
  yr <- as.POSIXlt(dt)$year + 1900
  ref_date <- as.POSIXlt(paste(yr, "-06-01", sep=""))
  return( difftime(dt, ref_date, units="days"))
}
t <- dayHourOfMonth(load_hist$dateTime)
load_hist <- data.frame(load_hist, dayHour=t)
load_hist <- data.frame(load_hist, year=as.POSIXlt(load_hist$dateTime)$year + 1900)

#plots of the month by year scenarios
ggplot(aes(x=as.numeric(dayHour), y = DEMAND, 
           color=year), data=load_hist) + 
  geom_point() + 
  geom_line() + facet_grid(year~.)
  theme(legend.position = "none")

#isolate the weekdays
#compute the 24 x 24 covariance matrix
t = cast(load_hist, subset=!isWkEnd, 
     formula= Date ~ Hour, value="DEMAND" )

#drop the first row and its annoying na ness
filt = is.na(t[, 1])
t = t[!filt, ]

filt = is.na(t[, 2])
t = t[!filt, ]

filt = is.na(t[, 25])
t = t[!filt, ]

colMeans(t[, 2:25])
cov(t[ 2:25])
df_temp = data.frame(Hour=1:24, Means=colMeans(t[, 2:25]), 
                     Date="06-01-2000")

# Graph a single peak day with the means
g <- ggplot(aes(x=Hour, y = DEMAND, group=Date), data=load_hist) + 
  geom_line(aes(alpha=.1, linetype='dotted')) + geom_point(aes(alpha=.1)) +
  theme_bw() + theme(legend.position="none") 
g + geom_line(aes(linetype='solid', color='b',y=Means), data=df_temp)

##  
write.csv(df_temp, "RawMeans_JunePK.csv")
write.csv(cov(t[, 2:25]), "RawCoVar_JunePK.csv")
##

##perform a temperature regression
## for the peak hours
load.hourly <- cast(load_hist, subset=!isWkEnd, 
                    formula= Date ~ Hour, value="DEMAND" )
load.hourly <- data.frame(load.hourly, PkLoad=rowMeans(load.hourly[, 8:24]))

temp.hourly <- cast(load_hist, subset=!isWkEnd, 
                    formula= Date ~ Hour, value="DryBulb" )
temp.hourly <- data.frame(temp.hourly, PkTemp=rowMeans(temp.hourly[, 8:24]))

resid.hourly <- load.hourly
load.hourly <- data.frame(load.hourly, PkTemp = temp.hourly$PkTemp)

ggplot(aes(x=PkTemp, y=X24), data=load.hourly) + 
  geom_point() + theme_bw() + 
  stat_smooth()

t_fun = function(load_col){ lm(load_col ~ PkTemp + I(PkTemp^2), data=load.hourly)}
temp.regs = apply(load.hourly[, 2:25], MARGIN=2, FUN=t_fun)

#this is a little dirty because we're doing everything in sample,
#but whatever....

pred.hourly <- load.hourly[, 1:25]
for (i in 1:24)
{
  pred.hourly[, i+1] <- predict(temp.regs[[i]], load.hourly )  
}

i = 13
ggplot(aes(x=PkTemp, y=load.hourly[, i+1]), data=load.hourly) + 
  geom_point() + theme_bw() + 
  geom_line(aes(y=pred.hourly[, i+1]), colour='blue') + 
  theme(legend.position='none')

#form temperature normalized loads.
i=12
resid.hourly <- pred.hourly
resid.hourly[, 2:25] <- load.hourly[, 2:25] - pred.hourly[, 2:25]
ggplot(aes(x=load.hourly$PkTemp, y=resid.hourly[, i+1]), data=resid.hourly) + 
  geom_point() + theme_bw()

RESID_MAX = apply(resid.hourly[,2:25], 2, max, na.rm=TRUE)
RESID_MIN = apply(resid.hourly[, 2:25], 2, min, na.rm=TRUE)

#compute upper and lower bound on the temperature normalization function.
t = read.csv(file="2002_hourly_data.txt")
back_temp = data.frame(t[, 1:4], Temp=t$DRY.BULB)

for(yr in 2001:1983)
{
  file_name = paste(yr, "_hourly_data.txt", sep="")
  t = read.csv(file=file_name)
  back_temp = rbind(back_temp, data.frame(t[, 1:4], Temp=t$DRY.BULB))
}

# only use the summer months.
t <- back_temp[ back_temp$MON %in% c(6, 7, 8), ]

#recast this to be one observation per day
back_pk_temps <- cast(t, subset= (HOUR %in% 7:23),
            YEAR + MON + DAY ~., fun.aggregate=mean, value="Temp")

names(back_pk_temps)[4] <- "PkTemp"


#magic number for the threshold
#we divy up the var heavily favoring the temp 
#we also divy up delta heavily favoring the temp
#for the temps then up
indx1_temp = qbinom(size=1840, prob=1-.05, p=.2, lower.tail=FALSE)
TEMP_MIN = min(t)
TEMP_MAX = t[indx1_temp]

#now do the UCS split
indx2_temp = qbinom(size=1840, prob=1-.025, p = .1, lower.tail=FALSE)
TEMP_MIN2 = min(t)
TEMP_MAX2 = t[indx2_temp]

names(resid.hourly)
RESID_MEAN = colMeans(resid.hourly[, 2:25], na.rm=TRUE)
RESID_COV = cor(resid.hourly[, 2:25], use="pair")

#dump the releavnt shite to python
##  
write.csv(cov(resid.hourly[, 2:25], use="pair"), "TempCovar_SummerPK.csv")
##
