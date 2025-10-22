#----------------------------------- setup ------------------------------------
rm(list = ls()); graphics.off(); gc(); if (interactive()) try(suppressWarnings(windows(record = TRUE)), silent = TRUE)

suppressPackageStartupMessages({
  library(rgeos)
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

#---------------------------- helpers / utilities ------------------------------

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

## Extract variogram params from an autoKrige() result
extract_vgpars <- function(ak, y, m) {
  vm <- ak$var_model
  data.frame(
    Year  = y,
    Month = m,
    model = vm$model[2],
    nug   = vm$psill[1],
    psill = vm$psill[2],
    range = vm$range[2],
    kappa = vm$kappa[2]
  )
}

## Single kriging pass for one variable, one Year–Month
krige_one <- function(var, y, m, wq_sp, grid_sgdf, depth_ll) {
  # Filter and drop NA
  kdat <- wq_sp[wq_sp$Year == y & wq_sp$Month == m, ]
  kdat <- kdat[!is.na(kdat[[var]]), ]
  if (nrow(kdat) == 0) return(list(r = NULL, vg = NULL, ak = NULL))
  
  # autoKrige and rasterize
  fml  <- as.formula(paste0(var, "~ 1"))
  ak   <- autoKrige(fml, input_data = kdat, new_data = grid_sgdf,
                    model = c("Sph", "Exp", "Ste"), debug.level = 0)
  kr   <- raster(ak$krige_output)
  kr[kr < 0] <- 0
  names(kr) <- paste(month.abb[m], y, sep = "")
  
  # Reproject back to long/lat to restore 60x50 dims of depth map
  kr_ll <- projectRaster(kr, depth_ll, method = "ngb",
                         crs = "+proj=longlat +datum=WGS84")
  
  list(r = kr_ll, vg = extract_vgpars(ak, y, m), ak = ak)
}

## Batch kriging for a variable across all Year–Month, with plotting to a single PDF
krige_variable <- function(var, env_name, years_vec, ym_index,
                           wq_sp, grid_sgdf, depth_ll,
                           dir_out, variogram_pdf) {
  
  wq_krg_stack <- stack()
  vgpars_all   <- data.frame()
  
  pdf(variogram_pdf, onefile = TRUE)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  
  for (y in years_vec) {
    cat("Processing", env_name, y, "\n")
    months <- ym_index[[as.character(y)]]
    for (m in months) {
      cat("  month", m, "\n")
      res <- krige_one(var, y, m, wq_sp, grid_sgdf, depth_ll)
      if (!is.null(res$ak)) {
        # Same plot you were making before per Y–M
        plot(res$ak, sub = paste(month.abb[m], y))
      }
      if (!is.null(res$r)) wq_krg_stack <- addLayer(wq_krg_stack, res$r)
      if (!is.null(res$vg)) vgpars_all   <- rbind(vgpars_all, res$vg)
      rm(res); gc()
    }
  }
  
  list(stack = wq_krg_stack, vg = vgpars_all)
}

## Single IDW pass for one variable, one Year–Month
idw_one <- function(var, y, m, wq_sp, grid_sgdf) {
  wq_sub <- wq_sp[wq_sp$Year == y & wq_sp$Month == m, ]
  wq_sub <- wq_sub[!is.na(wq_sub[[var]]), ]
  if (nrow(wq_sub) == 0) return(NULL)
  fml    <- as.formula(paste0(var, "~ 1"))
  r      <- raster(idw(fml, wq_sub, grid_sgdf, idp = idw.wt, debug.level = 0))
  names(r) <- paste(month.abb[m], y, sep = "")
  r
}

## Batch IDW for a variable across all Year–Month
idw_variable <- function(var, years_vec, ym_index, wq_sp, grid_sgdf) {
  stk <- stack()
  for (y in years_vec) {
    cat("Processing", var, y, "\n")
    months <- ym_index[[as.character(y)]]
    par(mfrow = c(2, 3))
    for (m in months) {
      cat("  month", m, "\n")
      r <- idw_one(var, y, m, wq_sp, grid_sgdf)
      if (!is.null(r)) {
        plot(r, main = paste(var, y, m))
        stk <- addLayer(stk, r)
      }
      rm(r); gc()
    }
    par(mfrow = c(1, 1))
  }
  stk
}

#------------------------------- basemaps & data -------------------------------

## High-res depth map with salt marsh corrections
depth.out <- raster("./0-Data/habitats/processed/crm-salt-marsh-corrected/crm 60x50 18s 485m - saltmarsh corrected.gri")
plot(depth.out, colNA = "black")

## Water quality
wq.base <- fread("./0-Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv") |> as.data.frame()

#-------------------- dummy NE-corner replication (upstream proxy) -------------
## (logic preserved; streamlined & vectorized)
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
print(n.added.rows / 288)  # Should equal number of dummy points (3)

## Basemap bbox (unchanged)
lllon <- -83.0; lllat <- 29.1; urlon <- -83.25; urlat <- 29.4
bbox  <- c(lllon, lllat, urlon, urlat)

#----------------------- make WQ spatial & target grid -------------------------

wq <- wq.base2
coordinates(wq)    <- ~ Long + Lat
proj4string(wq)    <- CRS("+proj=longlat +datum=WGS84")
wq                 <- spTransform(wq, CRS("+proj=utm +zone=16 +datum=WGS84 +units=km"))
wq.yrs             <- sort(unique(wq.base$Year))           # years defined from original data
ym_index           <- year_month_index(wq)                 # map: Year -> Months

## Grid: project depth to UTM and use as prediction grid
wq.grid.ras <- projectRaster(depth.out, method = "ngb",
                             crs = "+proj=utm +zone=16 +datum=WGS84 +units=km")
wq.grid     <- as(wq.grid.ras, "SpatialGridDataFrame")
plot(wq.grid)

#=============================== KRIGING =======================================

# --- ensure output dirs exist & close any stray PDF devices ---
dir.krg <- "./1-Prepare-Ecospace-Inputs/water-quality/out/KRG/"
dir.idw <- "./1-Prepare-Ecospace-Inputs/water-quality/out/IDW/"
if (!dir.exists(dir.krg)) dir.create(dir.krg, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir.idw)) dir.create(dir.idw, recursive = TRUE, showWarnings = FALSE)

# Close any already-open PDF device to avoid collisions
while ("pdf" %in% names(dev.list())) try(dev.off(), silent = TRUE)

# --- years across full range in your data ---
wq.yrs <- sort(unique(wq$Year))

# --- variables to krige (name + maxval) ---
krig_vars <- list(
  list(env_name = "Salinity",    maxval = 40),
  list(env_name = "Temperature", maxval = 40),
  list(env_name = "DO",          maxval = 16),
  list(env_name = "FC",          maxval = 500),
  list(env_name = "TNP",         maxval = 2200)
)

# --- helper function to extract variables ---
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

# ============================== KRIGING LOOP ================================
for (cfg in krig_vars) {
  env_name <- cfg$env_name
  maxval   <- cfg$maxval
  modtype  <- "KRG"
  
  message("==> Starting kriging for ", env_name)
  
  outname_variogram <- file.path(dir.krg, paste0("Monthly_variogram_", env_name, ".pdf"))
  
  wq.krg.stack <- stack()
  vgpars.out   <- data.frame()
  
  # open PDF defensively
  pdf_opened <- FALSE
  try({
    pdf(outname_variogram, onefile = TRUE)
    pdf_opened <- TRUE
  }, silent = TRUE)
  if (!pdf_opened) stop("Failed to open PDF at: ", outname_variogram, 
                        "\nCheck write permissions or path.")
  
  for (y in wq.yrs) {
    cat("Processing", y, "\n"); flush.console()
    months <- sort(unique(wq$Month[wq$Year == y]))
    for (m in months) {
      cat("  month", m, "\n"); flush.console()
      
      # subset and drop NA for current variable
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
  
  # close the PDF we opened
  if (pdf_opened) try(dev.off(), silent = TRUE)
  
  # skip output if no layers (just in case)
  if (nlayers(wq.krg.stack) == 0) {
    warning("No layers produced for ", env_name, "; skipping write.")
    next
  }
  
  file_base <- outname(wq.krg.stack, dir.krg, env_name, modtype)
  writeRaster(wq.krg.stack, file_base, overwrite = TRUE)
  write.csv(vgpars.out, paste0(file_base, "_VGpars.csv"), row.names = FALSE)
  
  # Optional map PDF (your function)
  # pdf_map(wq.krg.stack, 'virid', dir.krg, env_name, modtype, maxval)
  
  assign(paste0(tolower(env_name), ".krg.stack"), wq.krg.stack, inherits = TRUE)
  rm(wq.krg.stack, vgpars.out); gc()
}


#======================= SIMPLE INVERSE DISTANCE WEIGHTING =====================

idw.wt <- 2  # (unchanged) higher values => less smooth / more local influence
dir.idw  <- "./1-Prepare-Ecospace-Inputs/water-quality/out/IDW/"

## Build IDW stacks via the generic helper
sal.idw.stack  <- idw_variable("Salinity",    wq.yrs, ym_index, wq, wq.grid)
temp.idw.stack <- idw_variable("Temperature", wq.yrs, ym_index, wq, wq.grid)
DO.idw.stack   <- idw_variable("DO",          wq.yrs, ym_index, wq, wq.grid)
TNP.idw.stack  <- idw_variable("TNP",         wq.yrs, ym_index, wq, wq.grid)
FC.idw.stack   <- idw_variable("FC",          wq.yrs, ym_index, wq, wq.grid)

## Write out IDW stacks (kept function & naming)
writeRaster(sal.idw.stack,  outname(sal.idw.stack,  dir.idw, "Salinity",    "IDW"), overwrite = TRUE)
writeRaster(temp.idw.stack, outname(temp.idw.stack, dir.idw, "Temperature", "IDW"), overwrite = TRUE)
writeRaster(DO.idw.stack,   outname(DO.idw.stack,   dir.idw, "DO",          "IDW"), overwrite = TRUE)
writeRaster(TNP.idw.stack,  outname(TNP.idw.stack,  dir.idw, "TNP",         "IDW"), overwrite = TRUE)
writeRaster(FC.idw.stack,   outname(FC.idw.stack,   dir.idw, "FC",          "IDW"), overwrite = TRUE)

