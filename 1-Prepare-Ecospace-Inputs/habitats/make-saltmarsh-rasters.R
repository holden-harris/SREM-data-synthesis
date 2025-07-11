rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library('raster')

dir.maps <- "./analyses/habitats/data/"

dir.subs = list.dirs(dir.maps,recursive=F)

par(mfrow=c(4,3),mar=c(1,1,1,1))
depth = raster(paste0(dir.maps,"crm 60x50 18s 485m"))
depth.hires = raster(paste0(dir.maps,"crm 360x300 3s 81m"))

shp = shapefile(paste0(dir.maps, "Suwanee_River_Water_Management_District_(SRWMD)_2016-2017_Land_Use.shp"))
sort(unique(shp$LEVEL_3_DE))

plot(depth,colNA='black')
plot(depth.hires,colNA='black')


#salt marsh------------------------------------------------------------------------------------------------------
shp.sltmsh = shp[shp$LEVEL_3_DE=='Saltwater Marshes',]
plot(crop(shp.sltmsh,extent(depth)))

ras.sltmsh = rasterize(shp.sltmsh,depth,field='SHAPEAREA')
ras.sltmsh[ras.sltmsh>0]=1  #convert shape area to presence
ras.sltmsh[is.na(ras.sltmsh[])] = 0
ras.sltmsh[is.na(depth)] = NA
plot(ras.sltmsh,colNA='black')

ras.sltmsh.hires = rasterize(shp.sltmsh,depth.hires,field='SHAPEAREA')
ras.sltmsh.hires[ras.sltmsh.hires>0]=1  #convert shape area to presence
ras.sltmsh.hires[is.na(ras.sltmsh.hires[])] = 0
ras.sltmsh.hires[is.na(depth.hires)] = NA
plot(ras.sltmsh.hires,colNA='black')

#depth layer is masking some salt marsh cells, correct that
depth.cor = depth
depth.cor[ras.sltmsh==1] = min(depth.cor[depth.cor>0])  #if salt marsh is present and depth is NA, set to depth to small value (0.1)

depth.hires.cor = depth.hires
depth.hires.cor[ras.sltmsh.hires==1] = min(depth.hires.cor[depth.hires.cor>0]) 

writeRaster(depth.hires.cor, "crm 360x300 3s 81m - salt marsh corrected")

#create new rasters
ras.sltmsh.cor = rasterize(shp.sltmsh,depth.cor,field='SHAPEAREA')
ras.sltmsh.cor[ras.sltmsh.cor>0]=1  #convert shape area to presence
ras.sltmsh.cor[is.na(ras.sltmsh.cor[])] = 0
ras.sltmsh.cor[is.na(depth.cor)] = NA
plot(ras.sltmsh.cor,colNA='black')

ras.sltmsh.hires.cor = rasterize(shp.sltmsh,depth.hires.cor,field='SHAPEAREA')
ras.sltmsh.hires.cor[ras.sltmsh.hires.cor>0]=1  #convert shape area to presence
ras.sltmsh.hires.cor[is.na(ras.sltmsh.hires.cor[])] = 0
ras.sltmsh.hires.cor[is.na(depth.hires.cor)] = NA
plot(ras.sltmsh.hires.cor,colNA='black')


#resample hires.cor
ras.sltmsh.cor2 = resample(ras.sltmsh.hires.cor,depth.cor)
plot(ras.sltmsh.cor2,colNA='black')

writeRaster(ras.sltmsh.cor,"saltmarsh SRWMD 2016-2017 60x50 - presence.asc",overwrite=T)
writeRaster(ras.sltmsh.cor2,"saltmarsh SRWMD 2016-2017 60x50 - bilinear resampled.asc",overwrite=T)
