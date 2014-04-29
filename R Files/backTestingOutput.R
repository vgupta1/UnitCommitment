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
#  names(dat)[c(4, 6, 8)] <- c("FixedCost", "VarCost", "Shed")
  dat$Date = as.Date(dat$Date)
  dat$TotCost = dat$FixedCost + dat$VarCost + dat$Shed
  return(dat)
}

nom.dat = loadData("nomSummary.txt")
rob.dat = loadData("robSummary.txt")

qplot(x=rob.dat$PredVarCost, y=nom.dat$VarCost) + geom_point()

#Bug in the Var Cost computation for Robust!

