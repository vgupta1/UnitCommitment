#####
# Analyzing backtesting output
####
library(ggplot2)
library(reshape)
library(plyr)
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/All/TestSet/")

#VG Fix the robust output to name it as FixedCost, VarCost and Shed
loadData = function( path )
{
  dat = read.delim(path)
  dat$Date = as.Date(dat$Date)
  dat$TotCost = dat$FixedCost + dat$VarCost + dat$Shed
  return(dat)
}

nom.dat = loadData("nomSummary.txt")
rob.dat = loadData("robSummary.txt")
qplot(x=rob.dat$PredVarCost, y=nom.dat$VarCost) + geom_point()

#### The affine on the smaller sets
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/Size_10/")
dat.aff = read.delim(file="affUCSbacktest.txt", sep="\t")
dat.aff = unique(dat.aff)
dat.aff$TotCost = dat.aff$FixedCost + dat.aff$VarCost + 5e3 * dat.aff$Shed

dat.rob = read.delim(file="robbacktest_testset.txt")
dat.rob$TotCost = dat.rob$FixedCost + dat.rob$VarCost + 5e3 * dat.rob$Shed
dat.rob$Method = "Robust"

dat = rbind(dat.aff, dat.rob)
rm(dat.aff)
rm(dat.rob)
dat$Method = factor(dat$Method, 
                    levels=c("Nominal", "Robust", "UBudget", "UCS"))

tab = ddply(dat, .variables=~Method, 
      .fun=function(df){
        c(sd(df$FixedCost), sd(df$VarCost), sd(df$Shed), sd(df$TotCost))})
names(tab)<- c("Method", "FixedCost", "VarCost", "Shed", "TotCost")

#Generate the histograms
#Tot Cost
library(grid)
ggplot(aes(x=TotCost * 1e-3), data=dat) +
  geom_histogram(aes(y=..density..)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Total Cost ($M)") + 
  theme_bw(base_size=18) + 
  theme(panel.margin = unit(.5, "cm")) 

#Shed
ggplot(aes(x=Shed), data=dat) +
  geom_histogram(aes(y=..density..)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Mismatch (GWh)") + 
  theme_bw(base_size=18)

#VarCost
ggplot(aes(x=VarCost * 1e-3), data=dat) +
  geom_histogram(aes(y=..density..)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Production Cost ($M)") + 
  theme_bw(base_size=18)

#Fixed Cost
ggplot(aes(x=FixedCost * 1e-3), data=dat) +
  geom_histogram(aes(y=..density..)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Startup Cost ($M)") + 
  theme_bw(base_size=18) + 
  theme(panel.margin = unit(.25, "cm")) 

dat.by.method = cast(dat, Date ~ Method, value="TotCost")

#Savings to Nominal Plots
#highlighted by shed amounts
sheds.by.method = cast(dat, Date~Method, value="Shed")
t = cbind(dat.by.method, sheds.by.method[, 2:5])
names(t)[6:9]<- paste(names(t)[6:9], "Shed", sep="_")

ggplot(aes(x=Nominal * 1e-3), data=t) + 
  geom_point(aes(y=(Nominal-UCS)*1e-3, size=UCS_Shed, color=UCS_Shed) ) +
#  geom_point(aes(y=(Nominal-UBudget)*1e-3, size=UBudget_Shed, color=UBudget_Shed) ) +
  scale_size_continuous(name="Mismatch (GWh)", 
                        breaks=c(5, 10, 15), 
                        limits=c(0, 16)) + 
  scale_color_continuous(name="Mismatch (GWh)", 
                         breaks=c(5, 10, 15), 
                         limits=c(0, 16)) +
  guides(colour=guide_legend() ) +
  ylim(-1, 41) + 
  xlab("Nominal ($M)") + 
  ylab("Savings ($M)") + 
  theme_bw(base_size=18) + 
  theme(legend.position="top")

### Solution Times
ggplot(aes(x=Method, y=Time/60, fill=Method), 
       data =dat)+
  geom_boxplot() + 
  ylab("(Min)") + xlab("") +  
  theme_bw(base_size=18) + 
  theme(legend.position="none")


## Compute the savings to nominal averages

dat.savings = cast(dat, Date~Method, value="TotCost")
dat.savings[, 3:5] = dat.savings[, 3:5] - dat.savings[,2]
dat.savings = melt(dat.savings, id.vars=Date, value.name="TotCost", variable.name="TotCost")
names(dat.savings)[2] = "TotCost"
sum_fun = function(df)
{
  c(mean(df$TotCost, na.rm=TRUE), sd(df$TotCost, na.rm=TRUE)/length(df$TotCost))
}


tab = ddply(dat.savings, .variables=~Method, .fun=sum_fun)
names(tab)<- c("Method", "Mean", "SD")

