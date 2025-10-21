rm(list=ls())
windows()
library(raster)
library(reshape2)
library(rgeos)
library(sf)

## Import shape file
shp.clams = shapefile("./Data/habitats/inputs/FDACs-clams-data/Suwannee_River_AUZs.shp")
shp.clams
plot(shp.clams)

sf.clams = st_as_sf(shp.clams); sf.clams
unique(sf.clams$PARCEL_NAM) ## Look at parcels

## Note: From Leslie: Derricks is still on the books but is not currently leased, and neither is Pine Island. 
sf.clams = subset(sf.clams, sf.clams$PARCEL_NAM != "DERRICKS")
unique(sf.clams$PARCEL_NAM) ## Confirm DERRICKS is removed

## Rasterize
r1 = raster(extent(sf.clams), nrow=60, ncol=50); r1
r.clams.tic = rasterize(shp.clams, r1, field = as.numeric(shp.clams$TIC), fun='sum', 
                    background = 0, update = TRUE)

r.clams.area = rasterize(shp.clams, r1, field = 'ACRES', fun='sum', 
                    background = 0, update = TRUE)

par(mfrow=c(1,2))
plot(r.clams.area); plot(r.clams.tic)

## Re-project onto depth map
depth = raster("./final depth 60x50 18s 485m.asc")
depth
crs(depth) = "+proj=robin"
proj4string(depth) = CRS("+proj=longlat +datum=WGS84")

## Reproject rasterized clams
r.clams.tic[r.clams.tic>0.001] = 1 ## Make uniform

# sum all values that equal 1, and divide by the total number of pixels that are being aggregated (34*34)
r.clams.downscaled = aggregate(r.clams.tic, fact = 3, fun = function(x, ...) {(sum(x == 1, ...)/(20*20))*100})

par(mfrow=c(1,2))
plot(r.clams.tic); plot(r.clams.downscaled)

r.clams2 = projectRaster(r.clams.downscaled, depth, method = "bilinear")
r.clams2[r.clams2<0.01] = NA

par(mfrow=c(1,3))
plot(depth, colNA='black'); plot(r.clams.tic); plot(r.clams2)

writeRaster(r.clams2, "./Clam-leases-2022-09", format='ascii', overwrite=TRUE)