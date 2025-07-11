# Coastal Relief Model Bathymetry Maps at Multiple Resolutions

---

## Overview

This script processes bathymetry data from the Coastal Relief Model (CRM) and generates raster maps at multiple spatial resolutions. These outputs are used in Ecospace, the spatial module of Ecopath with Ecosim, where depth influences species distribution and movement.

---

## Load Libraries

```r
rm(list = ls()); graphics.off(); gc()
library(raster)
library(sp)
library(rgdal)
<<<<<<< HEAD
```

---

## Load and Preprocess Bathymetry

```r
# Replace land elevations (positive values) with NA
crm[crm > 0] <- NA

# Flip depths to positive values
crm <- crm * -1

# Replace 0s (interpreted as land in Ecospace) with the minimum depth value
crm[crm == 0] <- min(crm[crm > 0])
```

---

## Aggregate to Coarser Resolutions

```r
crm.agg2 <- aggregate(crm, fact = 2)
crm.agg2[crm.agg2 == 0] <- min(crm.agg2[crm.agg2 > 0])

crm.agg3 <- aggregate(crm, fact = 3)
crm.agg3[crm.agg3 == 0] <- min(crm.agg3[crm.agg3 > 0])

crm.agg4 <- aggregate(crm, fact = 4)
crm.agg4[crm.agg4 == 0] <- min(crm.agg4[crm.agg4 > 0])

crm.agg5 <- aggregate(crm, fact = 5)
crm.agg5[crm.agg5 == 0] <- min(crm.agg5[crm.agg5 > 0])

crm.agg6 <- aggregate(crm, fact = 6)
crm.agg6[crm.agg6 == 0] <- min(crm.agg6[crm.agg6 > 0])
```

---

## Estimate Cell Sizes in Meters

```r
dist1 <- pointDistance(coordinates(crm)[1:2,], longlat = TRUE)
dist2 <- pointDistance(coordinates(crm.agg2)[1:2,], longlat = TRUE)
dist3 <- pointDistance(coordinates(crm.agg3)[1:2,], longlat = TRUE)
dist4 <- pointDistance(coordinates(crm.agg4)[1:2,], longlat = TRUE)
dist5 <- pointDistance(coordinates(crm.agg5)[1:2,], longlat = TRUE)
dist6 <- pointDistance(coordinates(crm.agg6)[1:2,], longlat = TRUE)
```

---

## Plot Bathymetry at Each Resolution

```r
par(mfrow = c(3, 2))
plot(crm, colNA = 'black', main = paste0('3 sec / ', round(dist1[2,1]), ' m'))
plot(crm.agg2, colNA = 'black', main = paste0('6 sec / ', round(dist2[2,1]), ' m'))
plot(crm.agg3, colNA = 'black', main = paste0('9 sec / ', round(dist3[2,1]), ' m'))
plot(crm.agg4, colNA = 'black', main = paste0('12 sec / ', round(dist4[2,1]), ' m'))
plot(crm.agg5, colNA = 'black', main = paste0('15 sec / ', round(dist5[2,1]), ' m'))
plot(crm.agg6, colNA = 'black', main = paste0('18 sec / ', round(dist6[2,1]), ' m'))
par(mfrow = c(1,1))
```

---

## Export Raster Maps as ASCII for Ecospace

```r
dir.out <- "./analyses/habitats/outputs/crm-bathy-resolutions/"

writeRaster(crm,      paste0(dir.out, 'crm ', dim(crm)[1],'x',dim(crm)[2],' 3s ', round(dist1[2,1]),'m.asc'), format = 'ascii', NAflag = 0, overwrite = TRUE)
writeRaster(crm.agg2, paste0(dir.out, 'crm ', dim(crm.agg2)[1],'x',dim(crm.agg2)[2],' 6s ', round(dist2[2,1]),'m.asc'), format = 'ascii', NAflag = 0, overwrite = TRUE)
writeRaster(crm.agg3, paste0(dir.out, 'crm ', dim(crm.agg3)[1],'x',dim(crm.agg3)[2],' 9s ', round(dist3[2,1]),'m.asc'), format = 'ascii', NAflag = 0, overwrite = TRUE)
writeRaster(crm.agg4, paste0(dir.out, 'crm ', dim(crm.agg4)[1],'x',dim(crm.agg4)[2],' 12s ', round(dist4[2,1]),'m.asc'), format = 'ascii', NAflag = 0, overwrite = TRUE)
writeRaster(crm.agg5, paste0(dir.out, 'crm ', dim(crm.agg5)[1],'x',dim(crm.agg5)[2],' 15s ', round(dist5[2,1]),'m.asc'), format = 'ascii', NAflag = 0, overwrite = TRUE)
writeRaster(crm.agg6, paste0(dir.out, 'crm ', dim(crm.agg6)[1],'x',dim(crm.agg6)[2],' 18s ', round(dist6[2,1]),'m.asc'), format = 'ascii', NAflag = 0, overwrite = TRUE)
```

---

## Summary

- Bathymetry data from the Coastal Relief Model was cleaned and scaled to positive depths.
- Raster resolution was aggregated from 3 seconds (approx. 90m) up to 18 seconds (approx. 540m).
- These maps are formatted for use in Ecospace to test the impact of spatial resolution on ecosystem dynamics.
=======
