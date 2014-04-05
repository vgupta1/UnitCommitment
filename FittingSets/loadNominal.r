### Analyzing the load nominal simulation

library(gdata)
library(ggplot)
library(reshape)
load.python.temp = read.csv("load_validate_set_hourly_temps.csv")

setwd("/Users/VGupta/Documents/Research/UnitCommittment/ISO-NE/HistoricalBackTest_1")
costs = read.csv("nominal_backtest_costs.csv")
sched = read.csv("nominal_backtest.csv")s



#rescale the costs and look at some histograms
costs = cbind(costs, relative=( 100 * (costs[, 3] / costs[1, 3] - 1)))
summary(costs$relative)
ggplot(aes(x=Variable), data=costs) + geom_histogram() + theme_minimal()

#add a color coding of the costs by temperature for the day
#load python temp can be read off a .csv file somewhere
#drop the nas
load.python.temp = na.omit(load.python.temp)
#off by two for the columns
load.python.temp = cbind(load.python.temp, 
                         DRY_PK = rowMeans(load.python.temp[ , 10:25]), 
                         DRY_OFF = rowMeans(load.python.temp[, c(3:9, 26)]), 
                         DRY_BASE = rowMeans(load.python.temp[,3:26]))

costs = cbind(costs, load.python.temp[, 27:29])
ggplot(aes(x=DRY_OFF, y=Variable), data=costs) + geom_point()


#make load.python.temp look more like sched
temp.sched = load.python.temp
names(temp.sched)[3:26] <- names(sched)[3:26]
temp.sched = cbind(temp.sched, Type=rep("Temp", 187))


#Interestingly, on the whole it looks like costs are from the fact that it was a
#colder than usual summer
sched.mel = melt(sched, id.vars=c("Date", "Type"))
temps.mel = melt(temp.sched, id.vars=c("Date", "Type") )



sched.mel = rbind(sched.mel, temps.mel)
ggplot(aes(x=variable, y=value, group=Date), 
       data=subset(sched.mel, Date != "Planned" & Type=="Total")) + 
  geom_line(aes(color=)) + 
  geom_point() + 
  geom_line( data = subset(sched.mel, Date=="Planned" & Type=="Total"), color="red") + 
  theme_minimal()


