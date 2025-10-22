#-----------------------------------setup---------------------------------------
rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library(raster)
library(dplyr)
library(viridis)
library(stringr)


#################################################################################
##
## FUNCTIONS TO MAKE PDFS

## Function to make file name---------------------------------------------------
outname = function(wq.idw.stack, dir, env_name, modtype){
  #wq.idw.stack = stack_sal
  strt = str_sub(names(wq.idw.stack)[1], -7)
  stop = str_sub(names(wq.idw.stack)[nlayers(wq.idw.stack)], -7)
  out = gsub("X","", paste0(dir, env_name, "_", modtype, "_", strt,"-", stop))
  return(out)
}

outname(stack_sal, dir = dir.out, env_name = "Salinity", modtype = "Linkage")

## Function to make PDF of maps by year-----------------------------------------
pdf_map = function(plt.stack, colscheme, dir, env_name, modtype, maxval){
  
#  plt.stack = stack_nutr
#  colscheme = 'virid'
#  dir = dir.out
#  env_name = "Nutrients"
#  modtype = "Linkage model"
#  maxval = ceiling(max(wq$Salinity, na.rm = T))
#  maxval = 2200
  
  ## Determine color scheme
  brks = seq(0, maxval, maxval/50)
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
    #y = 1997
    print(paste("Plotting",y))
    raster_years = str_sub(names(plt.stack), -4) ## Get the last 4 digits
    yr.idx  = which(raster_years == y)
    plt.yr  = plt.stack[[yr.idx]] 
    par(mfrow=c(4,3),mar=c(1,1,2,0),oma=c(2,2,0,6))
    ## Plot 12 months
    for(i  in 1:nlayers(plt.yr)){
      plot(plt.yr[[i]], legend=F, col=color, colNA='darkgray', zlim=c(0, maxval),breaks=brks,
           main = str_sub(names(plt.yr)[i], -7))
    }
    ##Add legend
    par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,1))
    fields::image.plot(plt.yr,legend.only=T,zlim=c(0,maxval),col=color,add=T,legend.width=1,
               legend.mar=4,legend.line=3,
               legend.lab = env_name)
  }
  dev.off()
}


#################################################################################
##
## Ready hydro outputs

dir_sal = "./in/Hydro model outputs/Salinity_new/"
f_sal = list.files(path = dir_sal, pattern = '.*\\.(tif|asc)$')
r_sal = lapply(paste0(dir_sal,f_sal), raster)
stack_sal = stack(r_sal)


dir_temp = "./in/Hydro model outputs/Temperature_new/"
f_temp = list.files(path = dir_temp, pattern = '.*\\.(tif|asc)$')
r_temp = lapply(paste0(dir_temp, f_temp), raster)
stack_temp = stack(r_temp)


dir_nutr = "./in/Hydro model outputs/Nutrients_new/"
f_nutr = list.files(path = dir_nutr, pattern = '.*\\.(tif|asc)$')
r_nutr = lapply(paste0(dir_nutr, f_nutr), raster)
stack_nutr = stack(r_nutr)




#################################################################################
##
## Output maps

dir.out = "./out/Linkage model/"
wq = read.csv("./in/Spatial-temp_Phys-Flow_xMonth.csv")
vec = c('virid', 'turbo', 'brks')

maxval_sal = ceiling(max(wq$Salinity, na.rm = T))
pdf_map(stack_sal, colscheme = 'virid', dir = dir.out, "Salinity", "Linkage", maxval_sal)

maxval_temp = ceiling(max(wq$Temperature, na.rm = T))
pdf_map(stack_temp, colscheme = 'turbo', dir = dir.out, "Temperature", "Linkage", maxval_temp)

maxval_nutr = 2200
pdf_map(stack_nutr, colscheme = 'virid', dir = dir.out, "Nutrients", "Linkage", maxval_nutr)








