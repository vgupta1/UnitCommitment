###
# Extracts the load info from .xls sheets
#   Augments information as necessary
###
library(gdata)
library(ggplot2)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/Data/ISO-NE Load Data/")

#There is a set of data that is saved as .txt files
# Taken from: http://www.iso-ne.com/markets/hstdata/hourly/his_data_pre/index.html
# openData = function(year){ read.csv(paste(year, "_hourly_data.txt", sep="")) }
# #Up to 1998, have one column structure
# allData= lapply(1980:1998, openData)
# dat.Loads1 = do.call(rbind, allData)

#1998 to 2002, different structure
# allData = lapply(1999:2002, openData)
# dat.Loads2 = do.call(rbind, allData)

#combine them together....
#dat.Loads1 = rbind(dat.Loads1[, 1:7], dat.Loads2[1:7])
# dat.Loads1 = dat.Loads2[, 1:7]
# rm(list=c("allData", "dat.Loads2", "temp"))

#there is another set of data saved as .xls files
#Taken from: http://www.iso-ne.org/markets/hstdata/znl_info/hourly/index.html
openData = function(year){ read.xls(paste("smd_hourly_", year, ".xls", sep=""), sheet=2)}

## This call takes a long time!
allData = lapply(2003:2014, openData)

dat.Loads2 = do.call(rbind, allData)
rm(allData)

#unfortunately, these are very different columns.  Rearrange to match
dat.Loads2$Date = as.character(dat.Loads2$Date)  #should have done this while reading in.
dat.Loads2 = dat.Loads2[, c("Date", "Hour", "DEMAND", "DryBulb", "DewPnt")]
names(dat.Loads2) = c("Date", "HOUR", "LOAD", "DRY.BULB", "DEW.POINT")

# dat.Loads1 = cbind(Date= paste(dat.Loads1$YEAR, dat.Loads1$MON, dat.Loads1$DAY, sep="-"), dat.Loads1)
# dat.Loads1$Date = as.Date(dat.Loads1$Date, "%y-%m-%d")
# dat.Loads1$Date = as.character(dat.Loads1$Date)
# dat.Loads1 = dat.Loads1[, c("Date", "HOUR", "LOAD", "DRY.BULB", "DEW.POINT")]
# names(dat.Loads2) = names(dat.Loads1) 

dat.Loads = dat.Loads2
rm(list=c("dat.Loads1", "dat.Loads2"))

####
# Some sanity checks that the pre and post 2002 numbers are comparable
####
# dat.Loads$Date = as.Date(dat.Loads$Date)
# t = split(dat.Loads, list(getMonth(dat.Loads$Date), dat.Loads$Date < as.Date("2003-01-01")))
# t2 = lapply(t, function(df){mean(df$LOAD)})
# do.call(rbind, t2)
# t2 = data.frame(as.numeric(t2))
# t2 = cbind(Load=t2, Pre=c(rep("FALSE", 12), rep("TRUE", 12)), Month=c(1:12, 1:12))
# t2 = data.frame(t2)
# names(t2)[1] = "Load"
# ggplot(aes(x=Month, y=Load, group=Pre, color=Pre), data=t2) + 
#   geom_point() + geom_line()

#Unfortunately, pre-2002 Load is clearly different by a substantive amount
#Only save down the ata post 2002
write.csv(dat.Loads, "HourlyLoadData.csv")

