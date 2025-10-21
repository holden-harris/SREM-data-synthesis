#-----------------------------------setup---------------------------------------
rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library('raster')
library('rasterVis')
library('gstat')
library(dplyr); library(ggplot2); library(car); library(lme4); library(viridis)

color = viridis(30)

#################################################################################
##
## Predict nutrients

nutr = read.csv("./in/Nutrients_withFlow.csv", stringsAsFactors = TRUE)
nutr$Date = as.Date(nutr$Date)
nutr$logFlow = log(nutr$ma30)

## Visualize
par(mfrow=c(1,2))
boxplot(nutr$TNP ~ nutr$Site_ID)
palette(hcl.colors(3, "viridis"))
salform = formula(TNP ~ Sal)
plot(salform, data = nutr, ylim=c(0,2200), col=nutr$Source)
sal_lm = lm(salform, data=nutr); summary(sal_lm)
fr_lm  = lm(salform, data=subset(nutr, nutr$Source=="Frazer")); summary(fr_lm)
lw_lm  = lm(salform, data=subset(nutr, nutr$Source=="Lakewatch")); summary(lw_lm)
abline(sal_lm, col="blue", lwd=2)
abline(fr_lm, col=1, lwd=2)
abline(lw_lm, col =2)
par(mfrow=c(1,1))

## MIXED EFFECTS MODELS MODEL 
nutri_lmm = lmer(TNP ~ Sal + Temp + logFlow + (1|Site_ID), data=nutr); summary(nutri_lmm)


#################################################################################
##
## Get Stacks for Salinity, Temperature, and Flow

## Read in sal and temp stacks
sal.krg.stack = stack("./out/KRG/Salinity_KRG_Jan1997-Dec2020.gri")
temp.krg.stack = stack("./out/KRG/Temperature_KRG_Jan1997-Dec2020.gri")

## Get water data and make df of monthly flow averages
library(waterData)
dailyflow = importDVs(staid="02323500",code="00060",sdate="1997-01-01",edate="2020-12-31")
dailyflow$staid = NULL; dailyflow$qualcode = NULL
names(dailyflow) = c("Flow", "Date"); head(dailyflow); tail(dailyflow)
dailyflow$Date = as.Date(dailyflow$Date)
dailyflow$ym = as.factor(paste(format(dailyflow$Date, "%Y"), format(dailyflow$Date, "%m"), sep = "-"))
monthlyflow_logged = dailyflow %>% group_by(ym) %>% summarise(Flow = log(mean(Flow, na.rm=T)))
ym = monthlyflow_logged$ym


## Loop to make stack of nutrient raster layers--------------------------------
nutri_stack = stack()
for(i in 1:length(ym)){
#  i = 1
  print(paste("Processesing", ym[i], "flow was", round(exp(monthlyflow_logged$Flow[i]))))

  ## Make subsets for raster layers
  Sal = sal.krg.stack[[i]]
  Temp = temp.krg.stack[[i]]
  logFlow = raster::reclassify(Sal, cbind(-Inf, Inf, monthlyflow_logged$Flow[i]))
  monthstack = stack(Sal, Temp, logFlow)
  names(monthstack) = c("Sal","Temp","logFlow")
  
  ## Predict
  pred_nutri = raster::predict(monthstack, nutri_lmm, na.rm=T, re.form=NA)
  names(pred_nutri) = as.character(ym[i])
  pred_nutri
  par(mfrow=c(2,2)); plot(stack(monthstack, pred_nutri))
  nutri_stack = addLayer(nutri_stack, pred_nutri)
}

## Output--------------------------------------------------------------
nutri_stack
outname_nutrients =  paste0("./out/KRG/PDF Maps/Nutrients_Predicted_Jan1997-Dec2020")
writeRaster(nutri_stack, outname_nutrients, overwrite=TRUE)

## Make suffix names
sufx = paste(formatC(1:length(names(nutri_stack)), 
                     width=3, flag="0"), 
             names(sal.krg.stack), sep="_"); str(sufx)

## Write out ASCII files
dir.out = "./out/KRG/ASCII Maps/Nutrients/"
writeRaster(nutri_stack, paste0(dir.out, "Nutrients"), bylayer=T, suffix=sufx, format='ascii',overwrite=T)



#plot maps maps----------------------------------------------------------
#plt.stack = stack("./out/KRG/Nutrients_Predicted_Jan1997-Dec2020")
plt.stack = nutri_stack
wq.yrs = unique(nutr$Year)
peak = 2200
brks   = seq(0,peak,50)

dev.off()
pdf(paste0(outname_nutrients,'.pdf'),onefile=T)
for(y in wq.yrs){
#  y=wq.yrs[1]
  print(paste("Plotting",y))
  yr.idx  = which(substr(names(plt.stack),2,5) == y)
  plt.yr  = plt.stack[[yr.idx]] 
  par(mfrow=c(4,3),mar=c(1,1,2,0),oma=c(2,2,0,6))
  for(i  in 1:nlayers(plt.yr)){
 #   i=1
    plot(plt.yr[[i]],legend=F,col=color,colNA='darkgray',zlim=c(0,peak),breaks=brks,
         main=paste(substr(names(plt.yr)[i],2,5), substr(names(plt.yr)[[i]],7,8), sep="-"))
  }
  ## Add legend
  par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,1))
  fields::image.plot(plt.yr,legend.only=T,zlim=c(0,peak),col=color,add=T,legend.width=1,legend.mar=4,legend.line=3,
             legend.lab = 'Total Nutrients')
}
dev.off()
