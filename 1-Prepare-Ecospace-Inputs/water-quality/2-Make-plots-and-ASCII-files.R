#-----------------------------------setup---------------------------------------
rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library(raster)
library(dplyr)
library(tidyr)

## Depth map with salt marsh corrections
depth = raster('./in/depth/crm 60x50 18s 485m - saltmarsh corrected.gri')
plot(depth, colNA='black')


#################################################################################
##
## MAKE ASCII MAPS

## Read in KRG stacks-----------------------------------------------------------
sal.krg.stack  = stack("./out/KRG/Salinity_KRG_Jan1997-Dec2020.gri")
temp.krg.stack = stack("./out/KRG/Temperature_KRG_Jan1997-Dec2020.gri")
do.krg.stack   = stack("./out/KRG/DO_KRG_Jan1997-Dec2020.gri")
fc.krg.stack   = stack("./out/KRG/FC_KRG_Jan1997-Dec2020.gri")
nutr.krg.stack = stack("./out/KRG/Nutrients_Predicted_Jan1997-Dec2020.gri")

## Make suffix names
sufx = paste(formatC(1:length(names(sal.krg.stack)), 
                     width=3, flag="0"), 
             names(sal.krg.stack), sep="_"); head(sufx)

## Write out ASCII files--------------------------------------------------------

## Salinity
dir.out = "./out/KRG/ASCII Maps/Salinity/"
writeRaster(sal.krg.stack, paste0(dir.out, "Salinity"), bylayer=T, suffix=sufx, format='ascii',overwrite=T)

## Temperature
dir.out = "./out/KRG/ASCII Maps/Temperature/"
writeRaster(temp.krg.stack, paste0(dir.out, "Temperature"), bylayer=T, suffix=sufx, format='ascii',overwrite=T)

## DO
dir.out = "./out/KRG/ASCII Maps/DO/"
writeRaster(do.krg.stack, paste0(dir.out, "DO"), bylayer=T, suffix=sufx, format='ascii',overwrite=T)

## FC
dir.out = "./out/KRG/ASCII Maps/FC/"
writeRaster(fc.krg.stack, paste0(dir.out, "FC"), bylayer=T, suffix=sufx, format='ascii',overwrite=T)

## Nutrients
dir.out = "./out/KRG/ASCII Maps/Nutrients/"
writeRaster(nutr.krg.stack, paste0(dir.out, "FC"), bylayer=T, suffix=sufx, format='ascii',overwrite=T)


#################################################################################
##
## FUNCTIONS TO MAKE PDFS

## Function to make file name---------------------------------------------------
outname = function(wq.idw.stack, dir, env_name, modtype){
  out = gsub("X","", paste0(dir, env_name, "_", modtype, "_", names(wq.idw.stack)[1],"-",names(wq.idw.stack)[nlayers(wq.idw.stack)]))
  return(out)
}


## Function to make PDF of maps by year-----------------------------------------
pdf_map = function(plt.stack, colscheme, dir, env_name, modtype, maxval){
  ## Determine color scheme
  brks = seq(0,maxval,0.5)
  if(colscheme == 'turbo'){
     color   = viridis(min(length(brks),100), option = "H")
     colid   = "col-turbo"
  } else if(colscheme == 'virid'){
     color   = viridis(min(length(brks),100), option = "D")
     colid   = "col-viridis"
  } else if (colscheme == 'brks') {
     colv    = c("purple4","purple", "blue", "darkblue", "cyan", "green","darkgreen", "yellow", "orange", "red", "darkred")
     funpal  = colorRampPalette(colv,bias=2)
     nbcols  = length(brks)-1
     color   = funpal(nbcols) 
     colid   = "col-brks"
  }
  
  ## Make PDF  
  pdf(paste0(outname(plt.stack, dir, env_name, modtype), "_", colid, ".pdf"), onefile = T)
  wq.yrs = seq(1997, 2020, 1)
  
  for(y in wq.yrs){
    print(paste("Plotting",y))
    yr.idx  = which(substr(names(plt.stack),4,7) == y)
    plt.yr  = plt.stack[[yr.idx]] 
    par(mfrow=c(4,3),mar=c(1,1,2,0),oma=c(2,2,0,6))
    ## Plot 12 months
    for(i  in 1:nlayers(plt.yr)){
      plot(plt.yr[[i]],legend=F,col=color,colNA='darkgray',zlim=c(0,40),breaks=brks,
           main=paste(substr(names(plt.yr)[i],1,3),substr(names(plt.yr)[[i]],4,7)))
    }
    ##Add legend
    par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,1))
    image.plot(plt.yr,legend.only=T,zlim=c(0,maxval),col=color,add=T,legend.width=1,
               legend.mar=4,legend.line=3,
               legend.lab = env_name)
  }
  dev.off()
}


#################################################################################
##
## Output maps

dir.out = "./out/KRG/PDF Maps/"
wq = read.csv("./in/Spatial-temp_Phys-Flow_xMonth.csv")
vec = c('virid', 'turbo', 'brks')

## Salinity
maxval = ceiling(max(wq$Salinity, na.rm = T))
for(i in length(vec)) {
  pdf_map(sal.krg.stack, colscheme = vec[i], dir = dir.out, "Salinity", "KRG", maxval)
}

## Temperature
maxval = ceiling(max(wq$Temperature, na.rm = T))
for(i in vec) {
  pdf_map(temp.krg.stack, colscheme = i, dir = dir.out, "Temperature", "KRG", maxval)
}  

## DO
maxval = ceiling(max(wq$DO, na.rm = T))
for(i in vec) {
  pdf_map(do.krg.stack, colscheme = i, dir = dir.out, "DO", "KRG", maxval)
}  

## FC
maxval = ceiling(max(wq$FC, na.rm = T))
for(i in vec) {
  pdf_map(fc.krg.stack, colscheme = i, dir = dir.out, "FC", "KRG", maxval)
}  



###### Parkinglot
pdf_map(nutr.krg.stack, colscheme = 'virid', dir = dir.out, "Nutrients", "KRG", maxval)
pdf_map(nutr.krg.stack, colscheme = 'turbo', dir = dir.out, "Nutrients", "KRG", maxval)
pdf_map(nutr.krg.stack, colscheme = 'brks', dir = dir.out,  " Nutrients", "KRG", maxval)

