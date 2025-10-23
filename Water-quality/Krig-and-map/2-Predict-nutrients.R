#-----------------------------------setup---------------------------------------
rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library('raster')
library('rasterVis')
library('gstat')
library(dplyr); library(Matrix); library(car); library(lme4); library(viridis)
color = viridis(30)

#################################################################################
##
## Predict nutrients
nutr <- read.csv("./Data/water-quality/processed/Nutrients_withFlow.csv", stringsAsFactors = TRUE)
nutr$Date = as.Date(nutr$Date)
nutr$logFlow = log(nutr$ma30)

## Visualize
par(mfrow=c(1,2))
boxplot(nutr$TNP ~ nutr$Site_ID,
        xlab = "TNP", ylab = "Site")
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

## MIXED EFFECTS MODELS MODEL ==================================================

nutr$Site_ID <- as.factor(nutr$Site_ID) ## Ensure factor + drop NAs (keeps lmer happy)

## Fit model
nutri_lmm <- lme4::lmer(TNP ~ Sal + Temp + logFlow + (1|Site_ID), data = nutr)

## Model diagnostics, residuals vs fitted and qq
par(mfrow = c(1,2))
plot(fitted(nutri_lmm), resid(nutri_lmm), xlab="Fitted", ylab="Residuals"); abline(h=0, lty=2)
qqnorm(resid(nutri_lmm)); qqline(resid(nutri_lmm))
par(mfrow = c(1,1))
car::vif(lm(TNP ~ Sal + Temp + logFlow, data = nutr_clean)) ## Multicollinearity (VIF) on a plain LM proxy

## Summaries --------------------------------------------------------------------
summary(nutri_lmm)
as.data.frame(VarCorr(nutri_lmm)) ## Random-effects variance summary

## R2 (marginal/conditional) and ICC
r2_vals  <- performance::r2(nutri_lmm)
icc_vals <- performance::icc(nutri_lmm)
cat("\nModel R2 (marginal / conditional): ",
    round(r2_vals$R2_marginal, 3), " / ", round(r2_vals$R2_conditional, 3), "\n")
cat("ICC (Site_ID): ", round(icc_vals$ICC_adjusted, 3), "\n")


#################################################################################
##
## Get Stacks for Salinity, Temperature, and Flow

## Read in sal and temp stacks
sal.krg.stack  = stack("./Water-quality/Krig-and-map/out/KRG/Salinity_KRG_Jan1997-Dec2020.gri")
temp.krg.stack = stack("./Water-quality/Krig-and-map/out/KRG/Temperature_KRG_Jan1997-Dec2020.gri")

## Get water data and make df of monthly flow averages
phys <- read.csv("./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv")
monthlyflow <- phys %>% 
  group_by(YM) %>% 
  summarise(avg = mean(Flow, na.rm = TRUE)) %>% 
  as.data.frame()
monthlyflow$logged <- log(monthlyflow$avg)
head(monthlyflow)
ym <- monthlyflow$YM

## Loop to make stack of nutrient raster layers--------------------------------
nutri_stack = stack()
for(i in 1:length(ym)){
#  i = 1
  print(paste("Processesing", ym[i], "flow was", round(monthlyflow$avg[i])))

  ## Make subsets for raster layers
  Sal = sal.krg.stack[[i]]
  Temp = temp.krg.stack[[i]]
  logFlow = raster::reclassify(Sal, cbind(-Inf, Inf, monthlyflow$logged[i]))
  monthstack = stack(Sal, Temp, logFlow)
  names(monthstack) = c("Sal","Temp","logFlow")
  
  ## Predict
  pred_nutri = raster::predict(monthstack, nutri_lmm, na.rm=T, re.form=NA)
  names(pred_nutri) = as.character(ym[i])
  pred_nutri
  par(mfrow=c(2,2)); plot(stack(monthstack, pred_nutri))
  nutri_stack = addLayer(nutri_stack, pred_nutri)
}

## Check stack
nutri_stack
nlayers(nutri_stack) ## Should be 288
par(mfrow = c(4, 3), mar = c(1, 1, 2, 1)) ## Set up plotting window for 12 panels (4x3 grid)
for (i in 1:min(12, nlayers(nutri_stack))) { 
  plot(
    nutri_stack[[i]],
    main = names(nutri_stack)[i],
    col  = viridis::viridis(80),
    colNA = "darkgray"
  )
}
par(mfrow = c(1, 1))

## Output--------------------------------------------------------------
outname_nutrients =  paste0("./Water-quality/Krig-and-map/out/KRG/Nutrients_Predicted_Jan1997-Dec2020")
writeRaster(nutri_stack, outname_nutrients, overwrite=TRUE)