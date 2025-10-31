library(ggplot2)
library(ggmap)

#phys = read.csv("./out/Spatial-temp_Phys_xYM.csv")
phys = read.csv("./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv")

phys$Source= as.factor(phys$Source); summary(phys$Source)
phys$Year  = as.numeric(substr(phys$YM, 1,4))
phys$Month = as.factor(substr(phys$YM, 6,7))

## Basemap
#lllon = -82.95; lllat = 29.05; urlon = -83.3; urlat = 29.50
#bbox = c(lllon, lllat, urlon, urlat)


## Minimal working example: ESRI satellite, cropped to bbox
bbox_sf <- st_as_sfc(st_bbox(c(xmin = lllon, ymin = lllat, xmax = urlon, ymax = urlat), crs = 4326))
ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = 1) +
  ## tip: put coord_sf AFTER annotation_map_tile; datum=NA hides graticule mismatch
  coord_sf(xlim = c(lllon, urlon), ylim = c(lllat, urlat), expand = FALSE, datum = NA) +
  geom_sf(data = bbox_sf, fill = NA, color = NA) +
  labs(x = NULL, y = NULL)


library(maptiles)
library(terra)
library(sf)
library(tidyterra)
library(dplyr)
library(ggplot2)

## Study bbox (WGS84 lon/lat)
lllon <- -83.30; lllat <- 29.05
urlon <- -82.95; urlat <- 29.50
bbox_sf <- st_as_sfc(st_bbox(c(xmin = lllon, ymin = lllat, xmax = urlon, ymax = urlat), crs = 4326))

## Paths
dir.create("./data/habitats/basemaps", recursive = TRUE, showWarnings = FALSE)
dir.create("./data/habitats/basemaps/cache/tiles", recursive = TRUE, showWarnings = FALSE)
tif_out <- "./data/habitats/basemaps/suwannee_esri_z12.tif"

## Build once if missing; afterwards this block is skipped
if (!file.exists(tif_out)) {
  r_esri <- get_tiles(
    x        = bbox_sf,
#    provider = "Esri.WorldImagery",  ## satellite look
    provider = "Carto.LightNoLabels",   ## <---- changed provider here
    zoom     = 12,                   ## increase to 13–14 for crisper output
    crop     = TRUE,
    cachedir = "./data/cache/tiles"
  )
  writeRaster(r_esri, tif_out, overwrite = TRUE)
}

###################################################
## MAP SALINITY
## Page by year, faceted by month

## Load the baked basemap (fast)
r <- rast(tif_out)                 ## CRS: EPSG:3857
bbox_3857 <- st_transform(bbox_sf, 3857)
bb <- st_bbox(bbox_3857)

dev.off()
filename <- "./Water-quality/Compile-WQ-data/Salinity-maps-xMY_1997-2020"
pdf(paste0(filename, ".PDF"), onefile = TRUE)

## Plot years
for (g in 1997:2020) {
  g = 1997
  year <- g
  message(year)
  
  ## Convert points to sf and match basemap CRS (accurate + fast drawing)
  phys_sf <- phys %>%
    filter(Year == year, !is.na(Salinity)) %>%
    st_as_sf(coords = c("Long","Lat"), crs = 4326, remove = FALSE) %>%
    st_transform(3857)
  
  ## Basemap + points (faceted by Month)
  map_byYM <-
    ggplot() +
    geom_spatraster_rgb(data = r) +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      crs  = 3857,
      datum = NA,
      expand = FALSE
    ) +
    geom_sf(
      data = phys_sf,
      aes(
        shape = factor(Source, levels = c("FIM", "FDACS", "P-Coast", "LCR", "VBuoys")),
        color = Salinity
      ),
      size = 1, alpha = 0.8
    ) +
    xlab("") + ylab("") +
    scale_color_continuous(limits = c(0, 40), name = "Salinity") +
    ## 5 sources -> provide 5 shapes
    scale_shape_manual(values = c(19, 15, 17, 18, 8), name = "Source") +
    labs(title = year, x = "", y = "") +
    facet_wrap(~Month)
  
  print(map_byYM)
}

dev.off()


###################################################
## MAP SALINITY
## Page by year, faceted by month

dev.off()
filename = "./Water-quality/Compile-WQ-data/Salinity-maps-xMY_1997-2020"
#filename = "./maps/Salinity-maps-xMY_1997-1998"
pdf(paste0(filename, ".PDF"), onefile=T)

for(g in 1997:2020){
  #  test = 2010; for(g in test:test){
  year = g
  print(year)
  map_byYM <- 
    ggmap(baseMap) +
    geom_point(data = phys %>% filter (Year == year) %>% filter(!is.na(Sal_mean)),
               aes(shape = factor(Source, levels = c("FIM", "FDACS", "P-Coast", "LCR", "VBuoys")),
                   #      color = factor(Source, levels = c("FIM", "FDACS", "Frazer", "LCR")),
                   #      size = factor(log(n_days)+1),
                   x=Long, y=Lat, color = Sal_mean),
               size=1, alpha = 0.8) +
    xlab ("") + ylab("") +
    scale_color_continuous(limits = c(0, 40)) +
    scale_shape_manual(values = c(19, 15, 17, 18)) +
    labs(title = year, color = "Salinity", shape = "Source", x="", y="") +
    facet_wrap(~Month); print(map_byYM)
}  
dev.off()


###################################################
## MAP TEMPERATURE
## Page by year, faceted by month

dev.off()
filename = "./maps/Temperature-maps-xMY_1997-2020"
#filename = "./maps/Salinity-maps-xMY_1997-1998"
pdf(paste0(filename, ".PDF"), onefile=T)

for(g in 1997:2020){
  #  test = 2010; for(g in test:test){
  year = g
  print(year)
  map_byYM <- 
    ggmap(baseMap) +
    geom_point(data = phys %>% filter (Year == year) %>% filter(!is.na(Temp)),
               aes(shape = factor(Source, levels = c("FIM", "FDACS", "Frazer", "LCR", "VBuoys")),
                   x=Long, y=Lat, color = Temp),
               size=1, alpha = 0.8) +
    xlab ("") + ylab("") +
    scale_color_gradient(limits = c(2, 37), low="white", high="red2") +
    scale_shape_manual(values = c(19, 15, 17, 18, 3)) +
    labs(title = year, color = "Temperature (C)", shape = "Source", x="", y="") +
    facet_wrap(~Month); print(map_byYM)
}  
dev.off()
