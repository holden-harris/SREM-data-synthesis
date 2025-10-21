rm(list=ls());graphics.off();rm(.SavedPlots);windows(record=T)
library('raster')
library('gridExtra')
library('colorRamps')
library('RColorBrewer')

#####################################################################################
# SETUP
#####################################################################################

#set directories
dir.depth = paste0("./Data/habitats/processed/crm-resolutions/")

#read depth maps, one high res for rasterizing, and one for low res output
depth.3s = raster(paste0(dir.depth,"crm 360x300 3s 81m"))

#depth map with output resolution 
depth.out = raster(paste0(dir.depth,"\\crm 60x50 18s 485m"))

#####################################################################################
# LOAD AND EXPLORE SHAPEFILE
#####################################################################################
#From website: All data are provided in the Florida State Plane Projection, North Zone, Units US Feet, Datum HPGN (83/90). In some cases data may be Florida 
#State Plane Projection, North Zone, Units US Feet, Datum NAD83. Currently, most data are available as district wide layers. No customized clipping 
#or joining will be performed.

#seagrass 2001----------------------------------------------------------------------
sg01.shp <- shapefile("./Data/habitats/inputs/seagrass2001_SRWMD/seagrass01.shp")
sg01.shp; dim(sg01.shp); summary(sg01.shp); head(sg01.shp); 

#get and set proj4string based on description above; FL state plan projection is in the USAboundaries package data
#state_proj = data.frame(state_proj)
#state_proj$proj4_string[state_proj$state=='FL' & state_proj$zone=='north']

#change units to ft
FLproj = "+proj=lcc +lat_1=30.75 +lat_2=29.58333333333333 +lat_0=29 +lon_0=-84.5 +x_0=600000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=ft +no_defs"

#define projection of shapefile
proj4string(sg01.shp) = CRS(FLproj)

names(sg01.shp)
table(sg01.shp$FLUCSDESC,sg01.shp$FLUCSCODE)
#sg01.shp$AREA = as.numeric(as.factor(sg01.shp$FLUCSDESC))

spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Not Classified'),],colorkey=F,
       col='black',fill='tan',zcol='AREA',main='Seagrass, Patchy & Continuous',
       sp.layout=list(list(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Continuous'),],fill='darkgreen',col='darkgreen',first=F),
                      list(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Patchy'),],fill='green',col='green',first=F)))

#explore data
table(sg01.shp$FLUCSDESC)
grid.arrange(
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Bays and Estuaries','Not Classified'),],zcol='AREA',main='Bays and Estuaries'),
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Marine'),],zcol='AREA',main='Marine'),
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Not Classified'),],zcol='AREA',main='Not Classified'),
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Continuous','Not Classified'),],zcol='AREA',main='Seagrass, Continuous'),
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Patchy','Not Classified'),],zcol='AREA',main='Seagrass, Patchy'),
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Patchy','Seagrass, Continuous','Not Classified'),],zcol='AREA',main='Seagrass, All'),
  spplot(sg01.shp[sg01.shp$FLUCSDESC%in%c('Tidal Flats'),],zcol='AREA',main='Tidal Flats')
)


#####################################################################################
# RASTERIZE 
#####################################################################################
#get depth and seagrass on same CRS
crs(sg01.shp);crs(depth.3s)
sg01.shp = spTransform(sg01.shp,crs(depth.3s))

#convert shape area to presence
sg01.shp$AREA = 1

#rasterize
sg01.ras = rasterize(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Continuous','Seagrass, Patchy'),],depth.3s,field='AREA',fun='max')
sg01c.ras = rasterize(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Continuous'),],depth.3s,field='AREA',fun='max')
sg01p.ras = rasterize(sg01.shp[sg01.shp$FLUCSDESC%in%c('Seagrass, Patchy'),],depth.3s,field='AREA',fun='max')

#make depth cells NA and non-seagrass=0
sg01.ras[is.na(sg01.ras[])] = 0
sg01c.ras[is.na(sg01c.ras[])] = 0
sg01p.ras[is.na(sg01p.ras[])] = 0

sg01.ras[is.na(depth.3s[])] = NA
sg01c.ras[is.na(depth.3s[])] = NA
sg01p.ras[is.na(depth.3s[])] = NA

par(mfrow=c(1,3))
plot(sg01.ras,colNA='black')
plot(sg01c.ras,colNA='black')
plot(sg01p.ras,colNA='black')


#####################################################################################
# output
#####################################################################################
#match resolution of output depth map
sg01.ras.out = resample(sg01.ras,depth.out,method='bilinear')
sg01c.ras.out = resample(sg01c.ras,depth.out,method='bilinear')
sg01p.ras.out = resample(sg01p.ras,depth.out,method='bilinear')

rowColFromCell(sg01.ras.out,which.max(sg01.ras.out))

png("./Data/habitats/plots/seagrasses-continous-patch.png", 6.5, 6.5, units = "in", res = 200)
par(mfrow=c(2,3),mar=c(2,2,2,1))
plot(sg01.ras,colNA='black',main='all')
plot(sg01c.ras,colNA='black',main='continuous')
plot(sg01p.ras,colNA='black',main='patchy')
plot(sg01.ras.out,colNA='black',main='all')
plot(sg01c.ras.out,colNA='black',main='continuous')
plot(sg01p.ras.out,colNA='black',main='patchy')
dev.off()

dir.out = "./Data/habitats/processed/saltmarsh-rasters/"
sufx = paste0(nrow(depth.out),'x',ncol(depth.out))

writeRaster(sg01.ras.out, paste0(dir.out, 'seagrass01 all ',sufx,'.asc'),overwrite=T)
writeRaster(sg01c.ras.out,paste0(dir.out, 'seagrass01 continuous ',sufx,'.asc'),overwrite=T)
writeRaster(sg01p.ras.out,paste0(dir.out, 'seagrass01 patchy ',sufx,'.asc'),overwrite=T)

writeRaster(sg01.ras.out,paste0(dir.out, 'seagrass01 all ',sufx),overwrite=T)
writeRaster(sg01c.ras.out,paste0(dir.out, 'seagrass01 continuous ',sufx),overwrite=T)
writeRaster(sg01p.ras.out,paste0(dir.out, 'seagrass01 patchy ',sufx),overwrite=T)

sg.half = sg01.ras.out
sg.half[sg.half<0.5] = 0
hist(sg.half)
plot(sg.half,colNA='black',main='all')
writeRaster(sg.half,paste0(dir.out, 'seagrass all halved ',sufx,'.asc'),overwrite=T)
