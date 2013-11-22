## Build some advanced regression models on the load data set
library(gdata)
library(ggplot2)
library(reshape)

setwd("/Users/VGupta/Documents/Research/UnitCommittment/UnitCommitment/ISO-NE/")
load_hist.orig = read.csv('load_hist_orig2.csv')
load_hist.orig$dateTime <- as.POSIXct(load_hist.orig$dateTime) - 3600
load_hist.orig <- load_hist.orig[order(load_hist.orig$Date, load_hist.orig$Hour) , ]

#Limit Peak Days in June, July, August
load.summer <- subset(load_hist.orig, months(dateTime) %in% c("June", "July", "August") &
                        !isWkEnd )

#save the last 186 days for validation.  First 500 for training/testing:
# 2010-07-19

#drop 2013 for testing
load.validate <- subset(load.summer, as.Date(Date) >= "2010-07-09"  )
load.train <- subset( load.summer, as.Date(Date) < "2010-07-09")
rm(load_hist.orig, load.summer)

#reshape the training data to a form suitable for CART
#for now just use the hourly temperatures for the day.  
#should be replaced with forecasts
#should add previous load profile and indicator for day of week.
makeHourly = function(load_data)
{
  load.hour = cast(Date~Hour, data = load_data, value="DEMAND")
  drybulb.hour = cast(Date~Hour, data=load_data, value="DryBulb")
  dewpnt.hour = cast(Date~Hour, data = load_data, value="DewPnt")
  names(load.hour) = paste("Load", names(load.hour), sep="")
  names(drybulb.hour) = paste("Dry", names(drybulb.hour), sep="")
  names(dewpnt.hour) = paste("Dew", names(dewpnt.hour), sep="")
  load.hour = cbind(load.hour, 
                          drybulb.hour[, 2:25], 
                          dewpnt.hour[, 2:25])
}
load.validate.hour = makeHourly(load.validate)
load.train.hour = makeHourly(load.train)

write.csv(file="../load_validate.csv", load.validate.hour)

mean_load = colMeans(load.train.hour[, 2:25])
sd_loads = apply(load.validate.hour[,2:25], MARGIN=2, FUN=sd)
write.csv(file="../load_means.csv", mean_load)
write.csv(file="../load_std.csv", sd_loads)

#predict the hourly loads for any object that supports the predict interface
predictHourly = function(fitted_regs, load.hour)
{
  pred.hour = load.hour[, 1:25]  
  for(i in 1:24)
  {
    pred.hour[, i+1] <- predict(fitted_regs[[i]], load.hour)
  }
  return(pred.hour)
}

##Calculate the MSE
mse= function(sim, obs)
{
  filt = is.na(sim) | is.na(obs)
  mean((sim[!filt] - obs[!filt])^2)
}

MSEHourly = function(load.hour, pred.hour)
{
  out = rep(0, 24)
  for(ix in 1:24)
  {
    out[ix] <- mse(load.hour[, ix], pred.hour[, ix])
  }
  return(out)
}


#Fit the Trees
library(rpart)
library(rpart.plot)
fitTreeFun = function( load_col )
{
  trial_tree = rpart(Load ~ ., method="anova", 
                     data=cbind(Load = load_col, 
                                load.train.hour[, -(1:25)]), 
                     cp=0, xval=10)
  #try to identify 
  min_error = min(trial_tree$cptable[, "xerror"])
  min_indx = which(abs(trial_tree$cptable[, "xerror"] - min_error) < 1e-6)
  min_indx = as.numeric(min_indx[1])
  min_error_p = min_error + trial_tree$cptable[min_indx, "xstd"]
  min_indx = which(trial_tree$cptable[, "xerror"] - min_error_p < 1e-6)
  min_indx = as.numeric(min_indx[1])
  return(prune(trial_tree, cp = trial_tree$cptable[min_indx, "CP"]))
}
fittedTrees = apply(load.train.hour[, 2:25], MARGIN=2, FUN=fitTreeFun)
tree_preds = predictHourly(fittedTrees, load.validate.hour)
treesMSE = MSEHourly(tree_preds[,-1], load.validate.hour[,2:25])
treesMSErel = treesMSE/sd_loads^2

#write the fits on the validation set and the residuals
tree_resids = data.frame(LoadDate = load.validate.hour$LoadDate, 
                           load.validate.hour[, 2:25] - tree_preds[, -1])
write.csv(file="../treesFit.csv", tree_preds)
write.csv(file="../treesResid.csv", tree_resids)

#Fit Random Forest
library(randomForest)
#VG Figure out why it is that we have missing values
fitForestFun = function(load_col)
{
  fit <- randomForest(Load ~ .,   
                      data=cbind(Load=load_col, 
                                 load.train.hour[, -(1:25)]))
  return(fit)
}
fittedForests = apply(load.train.hour[, 2:25], MARGIN=2, FUN=fitForestFun)
forest_preds = predictHourly(fittedForests, load.validate.hour)
forestMSE = MSEHourly(forest_preds[, -1], load.validate.hour[, 2:25])
forestMSErel = forestMSE/sd_loads^2
forest_resids = data.frame(LoadDate = load.validate.hour$LoadDate, 
                           load.validate.hour[, 2:25] - forest_preds[, -1])
write.csv(file="../forestFit.csv", forest_preds)
write.csv(file="../forestResid.csv", forest_resids)

#Fit GAM
#simplest attempt, fit only against the current hour wet and dry
fitSimpleGAM = function( hr)
{
  form = paste("Load", hr, "~s(Dry", hr, ") + s(Dew", hr, ")", sep="")  
  return( gam(as.formula(form), data=load.train.hour))
}
fittedSimpleGAMS = apply(as.matrix(1:24), MARGIN=1, FUN=fitSimpleGAM)
simpleGAM_preds = predictHourly(fittedSimpleGAMS, load.validate.hour)
simpleGAMMSE = MSEHourly(simpleGAM_preds[, -1], load.validate.hour[, 2:25])
simpleGAMMSErel = simpleGAMMSE/sd_loads^2
simpleGAM_resids = data.frame(LoadDate = load.validate.hour$LoadDate, 
                           load.validate.hour[, 2:25] - simpleGAM_preds[, -1])
write.csv(file="../simpleGAMFit.csv", simpleGAM_preds)
write.csv(file="../simpleGAMResid.csv", simpleGAM_resids)

####
####
####
fitGAM = function( load_col )
{
  dry_terms = paste(paste("s(Dry", 1:24, ")", sep=""), collapse="+")
  wet_terms = paste(paste("s(Dew", 1:24, ")", sep=""), collapse="+")
  form = paste("Load", "~", dry_terms, "+", wet_terms, sep="")
  return( gam(as.formula(form), data=data.frame(Load=load_col, 
                                                load.train.hour[, -(1:25)]), 
              select=TRUE, method="REML"))
}

fittedGAMS = apply(load.train.hour[, 2:25], MARGIN=2, FUN=fitGAM)
GAM_preds = predictHourly(fittedGAMS, load.validate.hour)
GAMMSE = MSEHourly(GAM_preds[, -1], load.validate.hour[, 2:25])
GAMMSErel = GAMMSE/sd_loads^2
GAM_resids = data.frame(LoadDate = load.validate.hour$LoadDate, 
                              load.validate.hour[, 2:25] - GAM_preds[, -1])


#Fit LOESS



#### Assess the fits
#plot the fitted values in some way
tree_preds.melt = melt(tree_preds, id.vars="LoadDate", variable_name="Hour")
tree_preds.melt$Hour = unlist(strsplit(x=as.character(tree_preds.melt$Hour), split="Load", fixed=TRUE))[2 * 1:5688]
tree_preds.melt$Hour = as.numeric(tree_preds.melt$Hour)
t = strptime(paste(tree_preds.melt$LoadDate, " ", tree_preds.melt$Hour, ":00"), 
             "%Y-%m-%d %H")
tree_preds.melt = data.frame(tree_preds.melt, dateTime=t)
names(tree_preds.melt)[3] = "DEMAND"

ggplot(aes(dateTime, y=DEMAND), data=tree_preds.melt) +
  geom_line(color='blue') + theme_minimal() + 
  geom_point(data=load.validate) + 
  xlim(strptime("2011-07-01", "%Y-%m-%d"), 
       strptime("2011-09-10", "%Y-%m-%d"))

qplot(1:24, simpleGAMMSErel) + geom_line(color='red') + 
  geom_line(aes(y=treesMSErel), color="blue") + 
  geom_line(aes(y=forestMSErel), color="yellow") + 
  geom_line(aes(y=GAMMSErel), color="black") + 
  geom_point(aes(y=GAMMSErel), color="black")

