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
setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/results/Size_10/ValidateSet/")
dat = read.delim(file="backtest_full.txt", sep="\t")
dat$TotCost <- dat$FixedCost + dat$VarCost + 5e3 * dat$Shed
dat$OpCost <- dat$FixedCost + dat$VarCost

dat$Method=revalue(dat$Method, c(Robust="Box"))
dat$Method = factor(dat$Method, 
                       c("Nominal", "Box", "UBudget", "UCS"))

tab = ddply(dat, .variables=~Method, 
      .fun=colwise(mean, .cols=c("FixedCost", "VarCost", "Shed", "TotCost", "OpCost")))

tab2 = ddply(dat, .variables=~Method, 
                   .fun=colwise(function(col)sd(col)/sqrt(length(col)), 
                    .cols=c("FixedCost", "VarCost", "Shed", "TotCost", "OpCost")))
names(tab2)[2:6] = paste("Std", names(tab2[2:6]), sep="_")
tab = merge(tab, tab2)
tab[2:11] = 1e-3 * tab[2:11]
tab$Shed = tab$Shed * 1e3
tab$Std_Shed = tab$Std_Shed * 1e-3

rm(tab2) 

#Generate the Histograms
#Operating Cost Histogram
ggplot(aes(x=Method, y=OpCost*1e-3, fill=Method), data=tab) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=round(OpCost*1e-3, 2), vjust=-.5), y=0) + 
  ylab("Cost ($M)") + xlab("") + 
  geom_errorbar(aes(ymin = (OpCost - OpStd)*1e-3, 
                    ymax = (OpCost + OpStd)*1e-3), 
                width=.4, position="dodge") + 
  theme_bw(base_size=16) + 
  theme(legend.position="none")

#Operating Costs as histograms
ggplot(aes(x=OpCost*1e-3, fill=Method), 
       data = dat) +
  geom_histogram() + 
  facet_grid(Method~.) + 
  geom_vline(aes(xintercept=OpCost*1e-3), data=tab) +
  geom_text(aes(label=round(OpCost*1e-3, 2), 
                x=OpCost*1e-3 + .5), 
            data=tab, 
            y=75) +
  theme_bw(base_size=16) + 
  theme(legend.position="none") + 
  xlab("Operating Cost ($M)") + 
  ylab("")


#Mismatches as histograms
ggplot(aes(x=Shed, fill=Method), 
       data = dat) +
  geom_histogram() + 
  facet_grid(Method~.) + 
  geom_vline(aes(xintercept=Shed), data=tab) +
  geom_text(aes(label=round(Shed, 2), 
                x=Shed + .9), 
            data=tab, 
            y=75) +
  theme_bw(base_size=16) + 
  theme(legend.position="none") + 
  xlab("Mismatch (GWh)") + 
  ylab("")

#subdivide the barplots
dat$ShedCost = dat$Shed * 5e3
dat.melt = melt(dat, id.vars=c("Method"), 
                measure.vars=c("OpCost", "ShedCost"), 
                variable_name="CostType")

#num for each
num = table(dat$Method)[1] * 1e3
ggplot(aes(x=Method), 
       data=dat.melt) + 
  geom_bar(aes(weight=value /num, 
               fill=CostType)) + 
  xlab("") + ylab("Cost ($M)") + 
#   geom_text(aes(label=round(TotCost*1e-3, 2), 
#                 hjust=1.2, vjust = -.25, y=1e-3 * TotCost), 
#             data=tab) + 
  theme_bw(base_size=16) +
  theme(legend.position="top", 
        legend.title=element_blank()) + 
  scale_fill_discrete(label = c("Operating", "Mismatch")) +
  geom_errorbar(aes(y = TotCost * 1e-3, 
                    ymin = (TotCost - Std_TotCost), 
                    ymax = (TotCost + Std_TotCost)), 
                    data=tab,
                width=.2, position="dodge")
  



ggplot(aes(x=Shed, fill=Method),
       data = subset(dat, Method %in% c("Nominal", "UCS"))) +
  geom_histogram(aes(y=..density..), 
                 position="identity", alpha=.7) + 
  theme_bw(base_size=16) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .8)) + 
  xlab("Mismatch (GWh)") + ylab("")

#Breakeven calculation?

ggplot(aes(x=Method, y=TotCost*1e-3, fill=Method), data=tab) + 
  geom_bar(stat="identity") +
  ylab("Savings ($M)") + xlab("") + 
  theme_bw(base_size=16) + 
  theme(legend.position="none") 



#Tot Cost
library(grid)
ggplot(aes(x=TotCost * 1e-3), data=dat) +
  geom_histogram(aes(y=..density.., fill=Method)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Total Cost ($M)") + 
  geom_vline(aes(xintercept=TotCost), data=tab) +
  geom_text(aes(label=round(TotCost, 2), 
                x=TotCost + 7), 
            data=tab, 
            y=.3) +
    theme_bw(base_size=16) + 
  theme(legend.position="none")
  
  
  #  theme(panel.margin = unit(.5, "cm")) 

#Shed
ggplot(aes(x=Shed), data=dat) +
  geom_histogram(aes(y=..density.., fill=Method)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Mismatch (GWh)") + 
  geom_vline(aes(xintercept=Shed), data=tab) +
  geom_text(aes(label=round(Shed, 2), 
                x=Shed + 1.2), 
            data=tab, 
            y=1.5) +
  theme_bw(base_size=16) + 
  theme(legend.position="none")

#VarCost
ggplot(aes(x=VarCost * 1e-3), data=dat) +
  geom_histogram(aes(y=..density.., fill=Method)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Production Cost ($M)") + 
  geom_vline(aes(xintercept=VarCost), data=tab) +
  geom_text(aes(label=round(VarCost, 2), 
                x=VarCost + .5), 
            data=tab, 
            y=2) +
  theme_bw(base_size=16) + 
  theme(legend.position="none")

#Fixed Cost
ggplot(aes(x=FixedCost * 1e-3), data=dat) +
  geom_histogram(aes(y=..density.., fill=Method)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Startup Cost ($M)") + 
  geom_vline(aes(xintercept=FixedCost), data=tab) +
  geom_text(aes(label=round(FixedCost, 2), 
                x=FixedCost + .05), 
            data=tab, 
            y=30) +
  theme_bw(base_size=16) + 
  theme(legend.position="none")

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

