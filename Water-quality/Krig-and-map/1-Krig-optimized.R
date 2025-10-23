##---------------------------------- setup ------------------------------------
rm(list = ls()); graphics.off(); gc(); if (interactive()) try(suppressWarnings(windows(record = TRUE)), silent = TRUE)

suppressPackageStartupMessages({
  library(raster)
  library(rasterVis)
  library(gstat)
  library(data.table)
  library(fields)
  library(automap)
  library(viridis)
  library(dplyr)
  library(sp)
  library(ggmap)
})

color <- viridis(80)

## ========================= FUNCTIONS =========================================

## Compose output filename from a RasterStack (unchanged logic)
outname <- function(wq.idw.stack, dir, env_name, modtype) {
  paste0(
    dir, env_name, "_", modtype, "_",
    gsub("X", "", names(wq.idw.stack)[1]), "-",
    gsub("X", "", names(wq.idw.stack)[nlayers(wq.idw.stack)])
  )
}

## Safe add layer (skip if empty)
add_if_has_data <- function(stk, r) {
  if (!is.null(r) && nlayers(r) == 1) addLayer(stk, r) else stk
}

## Year→months index once (avoids re-computing)
year_month_index <- function(wq_df) {
  split(unique(wq_df$Month), f = wq_df$Year) |>
    lapply(sort)
}

## ========================= BASEMAPS & DATA ===================================

## High-res depth map with salt marsh corrections
depth.out <- raster("./Data/habitats/processed/crm-salt-marsh-corrected/crm 60x50 18s 485m - saltmarsh corrected.gri"); depth.out

## Water quality
wq.base <- fread("./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv") |> as.data.frame(); head(wq.base)

##-------------------- dummy NE-corner replication (upstream proxy) ------------
## 2025 update: streamlined & vectorized
ne.byYM <- wq.base |>
  dplyr::group_by(YM) |>
  dplyr::filter(Salinity == min(Salinity, na.rm = TRUE)) |>
  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

make_dummy <- function(df, lat, lon, label) {
  d <- df
  d$Lat <- lat; d$Long <- lon; d$Source <- label
  d
}
dumm1 <- make_dummy(ne.byYM, 29.395, -83.029, "Dummy1")
dumm2 <- make_dummy(ne.byYM, 29.375, -83.060, "Dummy2")
dumm3 <- make_dummy(ne.byYM, 29.345, -83.065, "Dummy3")

wq.base2 <- rbindlist(list(wq.base, dumm1, dumm2, dumm3), use.names = TRUE, fill = TRUE)

## Quick checks
n.added.rows <- nrow(wq.base2) - nrow(wq.base)
print(n.added.rows)
print(n.added.rows / 288) ##^ Should equal number of dummy points (i.e., should print "3")

## Basemap bbox (unchanged)
lllon <- -83.0; lllat <- 29.1; urlon <- -83.25; urlat <- 29.4
bbox  <- c(lllon, lllat, urlon, urlat)

##---------------------- make WQ spatial & target grid -------------------------
wq <- wq.base2
coordinates(wq)    <- ~ Long + Lat
proj4string(wq)    <- CRS("+proj=longlat +datum=WGS84")
wq                 <- spTransform(wq, CRS("+proj=utm +zone=16 +datum=WGS84 +units=km"))
wq.yrs             <- sort(unique(wq.base$Year))           ## years defined from original data
ym_index           <- year_month_index(wq)                 ## map: Year -> Months

## Grid: project depth to UTM and use as prediction grid
wq.grid.ras <- projectRaster(depth.out, method = "ngb",
                             crs = "+proj=utm +zone=16 +datum=WGS84 +units=km")
wq.grid     <- as(wq.grid.ras, "SpatialGridDataFrame")
plot(wq.grid)


##============================== KRIGING =======================================
# --- ensure output dirs exist & close any stray PDF devices ---
dir.krg <- "./Water-quality/Krig-and-map/out/KRG"
dir.idw <- "./Water-quality/Krig-and-map/out/IDW"
if (!dir.exists(dir.krg)) dir.create(dir.krg, recursive = TRUE)
if (!dir.exists(dir.idw)) dir.create(dir.idw, recursive = TRUE)

## Close any already-open PDF device to avoid collisions
while ("pdf" %in% names(dev.list())) try(dev.off(), silent = TRUE)

## --- years across full range in your data ---
wq.yrs <- sort(unique(wq$Year))

## --- variables to krige (name + maxval) ---
krig_vars <- list(
  list(env_name = "Salinity",    maxval = 40),
  list(env_name = "Temperature", maxval = 40),
  list(env_name = "DO",          maxval = 16),
  list(env_name = "FC",          maxval = 500),
  list(env_name = "TNP",         maxval = 2200)
)

## --- helper function to extract variables ---
extract_vgpars <- function(ak, y, m) {
  vm <- ak$var_model
  data.frame(
    Year  = y, Month = m,
    model = vm$model[2],
    nug   = vm$psill[1],
    psill = vm$psill[2],
    range = vm$range[2],
    kappa = vm$kappa[2]
  )
}

## ============================== KRIGING LOOP ================================
for (cfg in krig_vars) {
  env_name <- cfg$env_name
  maxval   <- cfg$maxval
  modtype  <- "KRG"
  
  message("==> Starting kriging for ", env_name)
  outname_variogram <- file.path(dir.krg, paste0("Monthly_variogram_", env_name, ".pdf"))
  wq.krg.stack <- stack()
  vgpars.out   <- data.frame()
  
  ## Open PDF defensively
  pdf_opened <- FALSE
  try({
    pdf(outname_variogram, onefile = TRUE)
    pdf_opened <- TRUE
  }, silent = TRUE)
  if (!pdf_opened) stop("Failed to open PDF at: ", outname_variogram, 
                        "\nCheck write permissions or path.")
  
  ## Outer loop through years
  for (y in wq.yrs) {
    cat("Processing", y, env_name, "\n"); flush.console()
    months <- sort(unique(wq$Month[wq$Year == y]))
    
    ## Inner loop through months
    for (m in months) {
      cat("  month", m, "\n"); flush.console()
      
      ## Subset and drop NA for current variable
      kdat <- wq[wq$Year == y & wq$Month == m, ]
      kdat <- kdat[!is.na(kdat[[env_name]]), ]
      if (nrow(kdat) == 0) next
      
      fml  <- as.formula(paste0(env_name, " ~ 1"))
      krig <- autoKrige(
        fml, input_data = kdat, new_data = wq.grid,
        model = c("Sph", "Exp", "Ste"), debug.level = 0
      )
      
      plot(krig, sub = paste(month.abb[m], y))
      
      krast <- raster(krig$krige_output)
      krast[krast < 0] <- 0
      names(krast) <- paste(month.abb[m], y, sep = "")
      krast2 <- projectRaster(
        krast, depth.out, method = "ngb",
        crs = "+proj=longlat +datum=WGS84"
      )
      
      wq.krg.stack <- addLayer(wq.krg.stack, krast2)
      vgpars.out   <- rbind(vgpars.out, extract_vgpars(krig, y, m))
      
      rm(kdat, krig, krast, krast2); gc()
    }
  }
  
  ## Close the PDF we opened
  if (pdf_opened) try(dev.off(), silent = TRUE)
  
  ## Skip output if no layers (just in case)
  if (nlayers(wq.krg.stack) == 0) {
    warning("No layers produced for ", env_name, "; skipping write.")
    next
  }
  
  writeRaster(wq.krg.stack, paste0(dir.krg, "/"), overwrite = TRUE)
  write.csv(vgpars.out, paste0(file_base, "_VGpars.csv"), row.names = FALSE)

  assign(paste0(tolower(env_name), ".krg.stack"), wq.krg.stack, inherits = TRUE)
  rm(wq.krg.stack, vgpars.out); gc()
}


##====================== SIMPLE INVERSE DISTANCE WEIGHTING =====================

idw.wt = 2  #higher values result in stronger influence of local points (less smooth)

## Make empty stacks
sal.idw.stack  = stack()
temp.idw.stack = stack()
DO.idw.stack   = stack()
TNP.idw.stack  = stack()
FC.idw.stack   = stack()

## IDW loop-------------------------------------------------------
par(mfrow=c(1,1))

## Outer loop over years
for(y in wq.yrs){
  #  y=1997
  months = sort(unique(wq$Month[wq$Year==y]))
  print(paste("Processing",y)); flush.console()
  
  ## Inner loop over months
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

## Write outputs -----------------------------------------------------------------
## IDW stacks
dir.out = "./Water-quality/Krig-and-map/out/IDW/"
writeRaster(sal.idw.stack, outname(sal.idw.stack, dir.out, "Salinity", "IDW"), overwrite=T)
writeRaster(temp.idw.stack, outname(temp.idw.stack, dir.out, "Temperature", "IDW"), overwrite=T)
writeRaster(DO.idw.stack, outname(DO.idw.stack, dir.out, "DO", "IDW"), overwrite=T)
writeRaster(TNP.idw.stack, outname(TNP.idw.stack, dir.out, "TNP", "IDW"), overwrite=T)
writeRaster(FC.idw.stack, outname(FC.idw.stack, dir.out, "FC", "IDW"), overwrite=T)

## Write out PDF Plots
#pdf_map(sal.idw.stack, dir.out, "Salinity", "IDW", 40)
#pdf_map(temp.idw.stack, dir.out, "Temperature", "IDW", 40)
#pdf_map(DO.idw.stack, dir.out, "DO", "IDW", 16)
#pdf_map(FC.idw.stack, dir.out, "FC", "IDW", 50)
#pdf_map(TNP.idw.stack, dir.out, "TNP", "IDW", 2200)
