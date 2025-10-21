#-----------------------------------setup---------------------------------------
  rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
#  .libPaths("C:\\R\\win-library")
  library(rgeos)
  library(raster)
  library(rasterVis)
  library(gstat)
  library(data.table)
  library(fields) 
  library(automap)
  library(viridis)
  library(dplyr)
  library(data.table)
  library(sp); library(ggmap)


color = viridis(80)

## Function to make file name
outname = function(wq.idw.stack, dir, env_name, modtype){
  out = gsub("X","", paste0(dir, env_name, "_", modtype, "_", names(wq.idw.stack)[1],"-",names(wq.idw.stack)[nlayers(wq.idw.stack)]))
  return(out)
}

###########################################################################################
##   Depth base map
depth.out = raster("./in/depth/crm 60x50 18s 485m - saltmarsh corrected.gri") ## High res depth map with salt marsh corrections  
depth.out; plot(depth.out, colNA='black')


###########################################################################################
##  WQ data
wq.base = read.csv("./in/Spatial-temp_Phys-Flow_xMonth.csv")

###########################################################################################
##
##Insert dummy data point NE corner that matches upstream point----------------------------

ne.df = data.frame()
for (ym in unique(wq.base$YM)){
#  ym = "1997-02"
  print(ym) 
  wq.sub = subset(wq.base, wq.base$YM == ym)
  dummy = subset(wq.sub, wq.sub$Salinity == min(wq.sub$Salinity, na.rm=T)) ## Replicate most upstream value recorded
  dummy$Source = "Dummy"
  ne.df = rbind(ne.df, dummy)
}

## Summarize data for one for YM
dumm = ne.df %>% group_by(YM) %>% summarize_all(mean)

## Replicate upstream data for three points
dumm1 = dumm; dumm1$Lat = 29.395; dumm1$Long = -83.029;  dumm1$Source="Dummy1"
dumm2 = dumm; dumm2$Lat = 29.375; dumm2$Long = -83.06;   dumm2$Source="Dummy2"
dumm3 = dumm; dumm3$Lat = 29.345;  dumm3$Long = -83.065; dumm3$Source="Dummy3"

## Append dummy data and check
wq.base2 = rbindlist(list(wq.base, dumm1, dumm2, dumm3)) ## Add dummy data to wq

## Basemap
lllon = -83.0; lllat = 29.1; urlon = -83.25; urlat = 29.4
bbox = c(lllon, lllat, urlon, urlat)
baseMap <- get_map(location=bbox, zoom = 10, source="google", maptype= "satellite", crop=FALSE) 

## Look at upstream points
lookmap = ggmap(baseMap) +
  geom_point(data = wq.base2, aes(x = Long, y = Lat, fill = Source, pch=Source), shape = 21); lookmap

## Check
n.added.rows = (nrow(wq.base2) - nrow(wq.base)); n.added.rows
n.added.rows / 288 ## Should = number of dummy data points (3)

###########################################################################################
##
## Make Water Quality Data Spatially Explicit-----------------------------------------------
wq = wq.base2
coordinates(wq) = ~Long+Lat
wq.yrs  = unique(wq.base$Year)
crs(wq) = "+proj=robin"
proj4string(wq) = CRS("+proj=longlat +datum=WGS84")
wq = spTransform(wq, CRS="+proj=utm +zone=16 +datum=WGS84 +units=km") 

## Make grid to write out spatial data output; Reproject data to UTM coordinates------------
wq.grid.ras = projectRaster(depth.out, method = 'ngb', crs= "+proj=utm +zone=16 +datum=WGS84 +units=km")
depth.out; wq.grid.ras
wq.grid = as(wq.grid.ras, "SpatialGridDataFrame")
summary(wq.grid)
plot(wq.grid) ## There's three 'island' water pixels that will be converted to land in ecosopace



##########################################################################################
##
##   KRIGING (SIMPLE ORDINARY)
##

###########################################################################################

graphics.off();rm(.SavedPlots);gc();windows(record=T)

## Salinity-----------------------------------------------------------------
env_name = "Salinity"
maxval = 40
modtype = "KRG"
dir.out = "./out/KRG/"
outname_variogram = paste0(dir.out, "Monthly_variogram_", env_name, ".pdf")

wq.krg.stack = stack()
vgpars.out = data.frame()
#dev.off()
pdf(outname_variogram,onefile=T)
for(y in wq.yrs){
  #    y=2018;m=6  
  print(paste("Processing",y)); flush.console()
  months = sort(unique(wq$Month[wq$Year==y]))
  for(m in months){
    print(paste("  month",m)); flush.console()
    kdat = wq[wq$Year==y & wq$Month==m,]
    kdat = kdat[!is.na(kdat$Salinity),]
    krig = autoKrige(Salinity~1,input_data=kdat, new_data= wq.grid,
                     model=c('Sph','Exp','Ste'), debug.level=0)
    plot(krig, sub=paste(month.abb[m],y))
    krast = raster(krig$krige_output)
    krast[krast<0] = 0 #; summary(krast)
    names(krast) = paste(month.abb[m],y,sep="")
    krast2 = projectRaster(krast, depth.out, method='ngb', crs="+proj=longlat +datum=WGS84") ## Reproject back to long lat to rstore 60x50 dimensions
    #plot(krast2, colNA='black')
    wq.krg.stack = addLayer(wq.krg.stack,krast2)
    vgpars.out = rbind(vgpars.out, data.frame(Year=y,Month=m, ## Collect model info
                        model=krig$var_model$model[2],
                        nug=krig$var_model$psill[1],
                        psill=krig$var_model$psill[2],
                        range=krig$var_model$range[2],
                        kappa=krig$var_model$kappa[2]))
    rm(kdat, krig);gc()
  }
}
dev.off()

## Output---------------------------------
writeRaster(wq.krg.stack, outname(wq.krg.stack, dir.out, env_name, modtype), overwrite=T)
sal.krg.stack = wq.krg.stack; rm(wq.krg.stack)
#pdf_map(wq.krg.stack, dir.out, env_name, modtype, maxval)


## Temperature-------------------------------------------------------------
env_name = "Temperature"
maxval = 40
modtype = "KRG"
dir.out = "./out/KRG/"
outname_variogram = paste0(dir.out, "Monthly_variogram_", env_name, ".pdf")

wq.krg.stack = stack()
vgpars.out = data.frame()
dev.off()
pdf(outname_variogram,onefile=T)
for(y in wq.yrs){
  print(paste("Processing",y)); flush.console()
  months = sort(unique(wq$Month[wq$Year==y]))
  for(m in months){
    #    y=1997;m=1  
    print(paste("  month",m)); flush.console()
    kdat = wq[wq$Year==y & wq$Month==m,]
    kdat = kdat[!is.na(kdat$Temperature),] ## Temperature
    krig = autoKrige(Temperature~1,input_data=kdat,new_data=wq.grid,model=c('Sph','Exp','Ste'),debug.level=0) ## Temperature
    plot(krig,sub=paste(month.abb[m],y))
    krast = raster(krig$krige_output)
    vgpars.out = rbind(vgpars.out,
                       data.frame(Year=y,Month=m,
                                  model=krig$var_model$model[2],
                                  nug=krig$var_model$psill[1],
                                  psill=krig$var_model$psill[2],
                                  range=krig$var_model$range[2],
                                  kappa=krig$var_model$kappa[2]))
    krast[krast<0] = 0
    names(krast) = paste(month.abb[m],y,sep="")
    krast2 = projectRaster(krast, depth.out, method='ngb', crs="+proj=longlat +datum=WGS84") ## Reproject back to long lat to rstore 60x50 dimensions
    wq.krg.stack = addLayer(wq.krg.stack,krast2)
    rm(kdat,krig);gc()
  }
}
dev.off()

## Output---------------------------------
file = outname(wq.krg.stack, dir.out, env_name, modtype)
writeRaster(wq.krg.stack, file, overwrite=T)
write.csv(vgpars.out,paste0(file, '_VGpars.csv'), row.names=F)


## DO-------------------------------------------------------------
env_name = "DO"
maxval = 16
modtype = "KRG"
dir.out = "./out/KRG/"
outname_variogram = paste0(dir.out, "Monthly_variogram_", env_name, ".pdf")

wq.krg.stack = stack()
vgpars.out = data.frame()
dev.off()
pdf(outname_variogram,onefile=T)
for(y in wq.yrs){
  print(paste("Processing",y)); flush.console()
  months = sort(unique(wq$Month[wq$Year==y]))
  for(m in months){
    #    y=1997;m=1  
    print(paste("  month",m)); flush.console()
    kdat = wq[wq$Year==y & wq$Month==m,]
    kdat = kdat[!is.na(kdat$DO),] ## DO
    krig = autoKrige(DO~1,input_data=kdat,new_data=wq.grid,model=c('Sph','Exp','Ste'),debug.level=0) ## DO
    plot(krig,sub=paste(month.abb[m],y))
    krast = raster(krig$krige_output)
    vgpars.out = rbind(vgpars.out,
                       data.frame(Year=y,Month=m,
                                  model=krig$var_model$model[2],
                                  nug=krig$var_model$psill[1],
                                  psill=krig$var_model$psill[2],
                                  range=krig$var_model$range[2],
                                  kappa=krig$var_model$kappa[2]))
    krast[krast<0] = 0
    names(krast) = paste(month.abb[m],y,sep="")
    krast2 = projectRaster(krast, depth.out, method='ngb', crs="+proj=longlat +datum=WGS84") ## Reproject back to long lat to rstore 60x50 dimensions
    wq.krg.stack = addLayer(wq.krg.stack,krast2)
    rm(kdat,krig);gc()
  }
}
dev.off()

## Output---------------------------------
file = outname(wq.krg.stack, dir.out, env_name, modtype)
writeRaster(wq.krg.stack, file, overwrite=T)
pdf_map(wq.krg.stack, dir.out, env_name, modtype, maxval)
write.csv(vgpars.out,paste0(file, '_VGpars.csv'), row.names=F)


## FC-------------------------------------------------------------
env_name = "FC"
maxval = 500
modtype = "KRG"
dir.out = "./out/KRG/"
outname_variogram = paste0(dir.out, "Monthly_variogram_", env_name, ".pdf")

wq.krg.stack = stack()
vgpars.out = data.frame()
dev.off()
pdf(outname_variogram,onefile=T)
for(y in wq.yrs){
  print(paste("Processing",y)); flush.console()
  months = sort(unique(wq$Month[wq$Year==y]))
  for(m in months){
    #    y=1997;m=1  
    print(paste("  month",m)); flush.console()
    kdat = wq[wq$Year==y & wq$Month==m,]
    kdat = kdat[!is.na(kdat$FC),] ## FC
    krig = autoKrige(FC~1,input_data=kdat, new_data = wq.grid,
                     model=c('Sph','Exp','Ste'),debug.level=0) ## FC
    plot(krig,sub=paste(month.abb[m],y))
    krast = raster(krig$krige_output)
    vgpars.out = rbind(vgpars.out,
                       data.frame(Year=y,Month=m,
                                  model=krig$var_model$model[2],
                                  nug=krig$var_model$psill[1],
                                  psill=krig$var_model$psill[2],
                                  range=krig$var_model$range[2],
                                  kappa=krig$var_model$kappa[2]))
    krast[krast<0] = 0
    names(krast) = paste(month.abb[m],y,sep="")
    krast2 = projectRaster(krast, depth.out, method='ngb', crs="+proj=longlat +datum=WGS84") ## Reproject back to long lat to rstore 60x50 dimensions
    wq.krg.stack = addLayer(wq.krg.stack,krast2)
    rm(kdat,krig);gc()
  }
}
dev.off()

## Output---------------------------------
file = outname(wq.krg.stack, dir.out, env_name, modtype)
writeRaster(wq.krg.stack, file, overwrite=T)
write.csv(vgpars.out,paste0(file, '_VGpars.csv'), row.names=F)


## TNP-------------------------------------------------------------
env_name = "TNP"
maxval = 500
modtype = "KRG"
dir.out = "./out/KRG/"
outname_variogram = paste0(dir.out, "Monthly_variogram_", env_name, ".pdf")

wq.krg.stack = stack()
vgpars.out = data.frame()
dev.off()
pdf(outname_variogram,onefile=T)
for(y in wq.yrs){
  print(paste("Processing",y)); flush.console()
  months = sort(unique(wq$Month[wq$Year==y]))
  for(m in months){
    #    y=1997;m=1  
    print(paste("  month",m)); flush.console()
    kdat = wq[wq$Year==y & wq$Month==m,]
    kdat = kdat[!is.na(kdat$TNP),] ## TNP
    krig = autoKrige(TNP~1,input_data=kdat,new_data=wq.grid,model=c('Sph','Exp','Ste'),debug.level=0) ## TNP
    plot(krig,sub=paste(month.abb[m],y))
    krast = raster(krig$krige_output)
    vgpars.out = rbind(vgpars.out,
                       data.frame(Year=y,Month=m,
                                  model=krig$var_model$model[2],
                                  nug=krig$var_model$psill[1],
                                  psill=krig$var_model$psill[2],
                                  range=krig$var_model$range[2],
                                  kappa=krig$var_model$kappa[2]))
    krast[krast<0] = 0
    names(krast) = paste(month.abb[m],y,sep="")
    krast2 = projectRaster(krast, depth.out, method='ngb', crs="+proj=longlat +datum=WGS84")
    wq.krg.stack = addLayer(wq.krg.stack,krast2)
    rm(kdat,krig);gc()
  }
}
dev.off()

## Output---------------------------------
file = outname(wq.krg.stack, dir.out, env_name, modtype)
writeRaster(wq.krg.stack, file, overwrite=T)
write.csv(vgpars.out,paste0(file, '_VGpars.csv'), row.names=F)





###########################################################################################
##
##  SIMPLE INVERSE DISTANCE WEIGHTING
##
###########################################################################################
idw.wt = 2  #higher values result in stronger influence of local points (less smooth)

## Make empty stacks
sal.idw.stack  = stack()
temp.idw.stack = stack()
DO.idw.stack   = stack()
TNP.idw.stack  = stack()
FC.idw.stack   = stack()

#IDW loop-------------------------------------------------------
for(y in wq.yrs){
  #  y=1997
  months = sort(unique(wq$Month[wq$Year==y]))
  print(paste("Processing",y)); flush.console()
  for(m in months){
    #    m=1  
    print(paste("  month",m)); flush.console()
    wq.sub = wq[wq$Year==y & wq$Month==m,]
    par(mfrow=c(2,3))
    
    ## Salinity
    sal.sub = wq.sub[!is.na(wq.sub$Salinity),]
    sal.idw  = idw(Salinity~1, sal.sub, wq.grid, idp=idw.wt, debug.level=0)
    sal.rast = raster(sal.idw)
    names(sal.rast) = paste(month.abb[m],y,sep="")
    plot(sal.rast, main=paste("Salinity", y, m))
    sal.idw.stack   = addLayer(sal.idw.stack, sal.rast)
    rm(sal.idw, sal.rast); gc()
    
    ## Temperature
    temp.sub = wq.sub[!is.na(wq.sub$Temperature),]
    temp.idw  = idw(Temperature~1, temp.sub, wq.grid, idp=idw.wt, debug.level=0)
    temp.rast = raster(temp.idw)
    names(temp.rast) = paste(month.abb[m],y,sep="")
    plot(temp.rast, main="Temperature")
    temp.idw.stack    = addLayer(temp.idw.stack, temp.rast)
    rm(temp.idw, temp.rast); gc()
    
    ## DO
    DO.sub = wq.sub[!is.na(wq.sub$DO),]
    DO.idw  = idw(DO~1, DO.sub, wq.grid, idp=idw.wt, debug.level=0)
    DO.rast = raster(DO.idw)
    names(DO.rast) = paste(month.abb[m],y,sep="")
    plot(DO.rast, main="DO")
    DO.idw.stack    = addLayer(DO.idw.stack, DO.rast)
    rm(DO.idw, DO.rast); gc()
    
    ## TNP
    TNP.sub = wq.sub[!is.na(wq.sub$TNP),]
    if(nrow(coordinates(TNP.sub))>0) { 
      TNP.idw  = idw(TNP~1, TNP.sub, wq.grid, idp=idw.wt, debug.level=0)
      TNP.rast = raster(TNP.idw)
      names(TNP.rast) = paste(month.abb[m],y,sep="")
      plot(TNP.rast, main="TNP")
      TNP.idw.stack    = addLayer(TNP.idw.stack, TNP.rast)
      rm(TNP.idw, TNP.rast); gc()
    }
    
    ## FC
    FC.sub = wq.sub[!is.na(wq.sub$FC),]
    if(nrow(coordinates(FC.sub))>0) { 
      FC.idw  = idw(FC~1, FC.sub, wq.grid, idp=idw.wt, debug.level=0)
      FC.rast = raster(FC.idw)
      names(FC.rast) = paste(month.abb[m],y,sep="")
      plot(FC.rast, main="FC")
      FC.idw.stack    = addLayer(FC.idw.stack, FC.rast)
      rm(FC.idw, FC.rast); gc()
    }
  }
}
par(mfrow=c(1,1))

## Output-----------------------------------------------------------------

## Write out the IDW stacks
dir.out = "./out/IDW/"
writeRaster(sal.idw.stack, outname(sal.idw.stack, dir.out, "Salinity", "IDW"), overwrite=T)
writeRaster(temp.idw.stack, outname(temp.idw.stack, dir.out, "Temperature", "IDW"), overwrite=T)
writeRaster(DO.idw.stack, outname(DO.idw.stack, dir.out, "DO", "IDW"), overwrite=T)
writeRaster(TNP.idw.stack, outname(TNP.idw.stack, dir.out, "TNP", "IDW"), overwrite=T)
writeRaster(FC.idw.stack, outname(FC.idw.stack, dir.out, "FC", "IDW"), overwrite=T)



## Write out PDF Plots
pdf_map(sal.idw.stack, dir.out, "Salinity", "IDW", 40)
pdf_map(temp.idw.stack, dir.out, "Temperature", "IDW", 40)
pdf_map(DO.idw.stack, dir.out, "DO", "IDW", 16)
pdf_map(FC.idw.stack, dir.out, "FC", "IDW", 50)
pdf_map(TNP.idw.stack, dir.out, "TNP", "IDW", 2200)














