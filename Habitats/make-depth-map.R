rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)

library(raster)
library(sp)
library(rgdal)

#read and summarize bathymetry data from coastal relief model
dir.in = "./Data/habitats/inputs/"
crm    = raster(paste0(dir.in, "coastal-relief-map.tiff"))
dim(crm)

#replace land with NA and make depth positive
crm[crm>0] = NA
crm = crm*-1
#0s will be treated as land in Ecospace, so set 0 to min value
crm[crm==0] = min(crm[crm>0])

#aggregate to coarser resolution
crm.agg2 = aggregate(crm,fact=2)
#0s will be treated as land in Ecospace, so set 0 to min value
crm.agg2[crm.agg2==0] = min(crm.agg2[crm.agg2>0])

crm.agg3 = aggregate(crm,fact=3)
crm.agg3[crm.agg3==0] = min(crm.agg3[crm.agg3>0])

crm.agg4 = aggregate(crm,fact=4)
crm.agg4[crm.agg4==0] = min(crm.agg4[crm.agg4>0])

crm.agg5 = aggregate(crm,fact=5)
crm.agg5[crm.agg5==0] = min(crm.agg5[crm.agg5>0])

crm.agg6 = aggregate(crm,fact=6)
crm.agg6[crm.agg6==0] = min(crm.agg6[crm.agg6>0])

#get cell size in meters
dist1 = pointDistance(coordinates(crm)[1:2,],longlat=T)
dist2 = pointDistance(coordinates(crm.agg2)[1:2,],longlat=T)
dist3 = pointDistance(coordinates(crm.agg3)[1:2,],longlat=T)
dist4 = pointDistance(coordinates(crm.agg4)[1:2,],longlat=T)
dist5 = pointDistance(coordinates(crm.agg5)[1:2,],longlat=T)
dist6 = pointDistance(coordinates(crm.agg6)[1:2,],longlat=T)

#plots
png("./Data/habitats/plots/crm-resolutions.png", height = 9, width = 6.5, units = "in", res = 300)
par(mfrow=c(3,2))
plot(crm,colNA='black',main=paste0('3 sec / ',round(dist1[2,1]),' m / ',dim(crm)[1],'x',dim(crm)[2]))
plot(crm.agg2,colNA='black',main=paste0('6 sec / ',round(dist2[2,1]),' m / ',dim(crm.agg2)[1],'x',dim(crm.agg2)[2]))
plot(crm.agg3,colNA='black',main=paste0('9 sec / ',round(dist3[2,1]),' m / ',dim(crm.agg3)[1],'x',dim(crm.agg3)[2]))
plot(crm.agg4,colNA='black',main=paste0('12 sec / ',round(dist4[2,1]),' m / ',dim(crm.agg4)[1],'x',dim(crm.agg4)[2]))
plot(crm.agg5,colNA='black',main=paste0('15 sec / ',round(dist5[2,1]),' m / ',dim(crm.agg5)[1],'x',dim(crm.agg5)[2]))
plot(crm.agg6,colNA='black',main=paste0('18 sec / ',round(dist6[2,1]),' m / ',dim(crm.agg6)[1],'x',dim(crm.agg6)[2]))
dev.off()

par(mfrow=c(1,1))
plot(crm.agg6,colNA='black',main=paste0('18 sec / ',round(dist6[2,1]),' m / ',dim(crm.agg6)[1],'x',dim(crm.agg6)[2]))

#export maps as ascii files for ecospace
dir.out <- "./Data/habitats/processed/crm-resolutions/"
writeRaster(crm,     paste0(dir.out, 'crm ',dim(crm)[1], 'x',dim(crm)[2],' 3s ',round(dist1[2,1]),'m.asc'),format='ascii',NAflag=0,overwrite=T)
writeRaster(crm.agg2,paste0(dir.out, 'crm ',dim(crm.agg2)[1],'x',dim(crm.agg2)[2],' 6s ',round(dist2[2,1]),'m.asc'),format='ascii',NAflag=0,overwrite=T)
writeRaster(crm.agg3,paste0(dir.out, 'crm ',dim(crm.agg3)[1],'x',dim(crm.agg3)[2],' 9s ',round(dist3[2,1]),'m.asc'),format='ascii',NAflag=0,overwrite=T)
writeRaster(crm.agg4,paste0(dir.out, 'crm ',dim(crm.agg4)[1],'x',dim(crm.agg4)[2],' 12s ',round(dist4[2,1]),'m.asc'),format='ascii',NAflag=0,overwrite=T)
writeRaster(crm.agg5,paste0(dir.out, 'crm ',dim(crm.agg5)[1],'x',dim(crm.agg5)[2],' 15s ',round(dist5[2,1]),'m.asc'),format='ascii',NAflag=0,overwrite=T)
writeRaster(crm.agg6,paste0(dir.out, 'crm ',dim(crm.agg6)[1],'x',dim(crm.agg6)[2],' 18s ',round(dist6[2,1]),'m.asc'),format='ascii',NAflag=0,overwrite=T)

writeRaster(crm,paste0(dir.out, 'crm ',dim(crm)[1],'x',dim(crm)[2],' 3s ',round(dist1[2,1]),'m'),NAflag=0,overwrite=T)
writeRaster(crm.agg2,paste0(dir.out, 'crm ',dim(crm.agg2)[1],'x',dim(crm.agg2)[2],' 6s ',round(dist2[2,1]),'m'),NAflag=0,overwrite=T)
writeRaster(crm.agg3,paste0(dir.out, 'crm ',dim(crm.agg3)[1],'x',dim(crm.agg3)[2],' 9s ',round(dist3[2,1]),'m'),NAflag=0,overwrite=T)
writeRaster(crm.agg4,paste0(dir.out, 'crm ',dim(crm.agg4)[1],'x',dim(crm.agg4)[2],' 12s ',round(dist4[2,1]),'m'),NAflag=0,overwrite=T)
writeRaster(crm.agg5,paste0(dir.out, 'crm ',dim(crm.agg5)[1],'x',dim(crm.agg5)[2],' 15s ',round(dist5[2,1]),'m'),NAflag=0,overwrite=T)
writeRaster(crm.agg6,paste0(dir.out, 'crm ',dim(crm.agg6)[1],'x',dim(crm.agg6)[2],' 18s ',round(dist6[2,1]),'m'),NAflag=0,overwrite=T)
