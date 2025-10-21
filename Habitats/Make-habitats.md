# Make static habitat layers for Ecospace

The primary update from **SREMv1** to **SREMv2** is that SREMv2 simulates *spatial–temporal* dynamics to represent habitats, the movement of species and fishers, and spatially explicit water conditions and primary production, with a specific emphasis on nutrient, seagrass, and fish dynamics.

Static habitat layers were developed for **depth**, **seagrass**, **oyster reefs**, **salt marsh/mangrove**, and **clam lease areas** (Figure 1E). These layers constrain the spatial extent of key habitats in the Ecospace model domain and are used to parameterize species–habitat affinities, dispersal rates, and environmental forcing.

---

## Depth map

**Code:** [`make-depth-map.R`](https://github.com/holden-harris/SREM-project/blob/main/1-Prepare-Ecospace-Inputs/habitats/make-depth-map.R)

**Source & domain.** The Ecospace basemap was derived from the **NOAA Coastal Relief Model (CRM)** (*NOAA National Centers for Environmental Information, 2023*). Bathymetry was cropped to the SREM spatial domain extending from **29.05° N to 29.40° N** and **82.95° W to 83.30° W**.

**Input.** `./Data/habitats/inputs/coastal-relief-map.tiff`

**Pre-processing rules (as used in the R code):**
- CRM elevations **> 0 m** (land) are set to **NA**.
- Ocean depths are made **positive** by multiplying by **–1** (so deeper water has larger values).
- Because Ecospace treats **0** as land, any 0 values are set to the **minimum positive depth** in the raster (done **after each** aggregation to avoid accidental 0s).
  
**Aggregation & resolution.**
- The native CRM (3 arcsec) was aggregated by `fact = 2…6` to generate 6, 9, 12, 15, and **18 arcsec** candidate grids.
- Great-circle cell sizes were computed using `sp::pointDistance(..., longlat = TRUE)` to label plots/exports in meters.
- The **18 arcsec (0.005°)** grid (≈ **485 m × 485 m**) was selected and written as the final Ecospace base map.  
  *In our domain, the 18″ grid yields a 60 × 50 (rows × cols) basemap.*

**Visualization**
- A 6-panel figure comparing resolutions is written to  
  `./Data/habitats/plots/crm-resolutions.png`  
  with each panel labeled by arcseconds, cell size (m), and raster dimensions.
- A single-panel plot of the selected **18″** grid is also produced.

**Exports for Ecospace.**
- ASCII (`format = "ascii"`, `NAflag = 0`) and default raster formats are written for each resolution to  
  `./Data/habitats/processed/crm-resolutions/`,  
  with filenames that encode grid dimensions, arcseconds, and meters (e.g.,  
  `crm 60x50 18s 485m.asc`).

---

## Salt marsh / mangrove habitat

**Code:** [`make-saltmarsh-rasters.R`](https://github.com/holden-harris/SREM-data-synthesis/blob/main/Habitats/make-saltmarsh-rasters.R)

The marsh/mangrove layer constrains where estuarine fringing habitats occur in Ecospace and is used to set habitat affinities for resident fishes and invertebrates.

**Source & domain.**  
Polygons originate from the **SRWMD 2016–2017 Land Use** dataset (distributed by FDEP). Features are cropped to the SREM domain (**29.05–29.40° N**, **82.95–83.30° W**) and rasterized onto the same grid as the Ecospace depth layer.

**Inputs:**
- Depth rasters (for grid/NA mask):  
  `./Data/habitats/processed/crm-resolutions/crm 60x50 18s 485m`  
- Land-use shapefile:  
  `Data/habitats/inputs/SWMD-land-use/SRWMD_2016-2017_Land_Use.shp`

**Pre-processing rules (as used in the R code):**
- Select **salt marsh** polygons as `LEVEL_3_DE == "Saltwater Marshes"`.  
  *(Apply the same workflow to mangroves by selecting the corresponding `LEVEL_3_DE` class, e.g., “Mangrove Swamps”, if present.)*
- Rasterize polygons to both the **base grid (60×50; 18″)** and a **hi-res grid (360×300; 3″)** using `raster::rasterize(..., field = "SHAPEAREA")`.
- Convert any positive rasterized value to **presence (1)**; set non-intersecting cells to **0**; set cells that fall on **depth NA** to **NA** to keep land masked consistently.

**Grid alignment & depth-mask correction.**
- In some shoreline cells, valid marsh polygons overlapped locations where the depth raster was **NA** (land), which would incorrectly erase marsh presence when masking.  
- To correct this, copy the depth rasters and, **where marsh==1**, set depth to the **minimum positive depth** observed in that raster (i.e., a tiny positive value) so those shoreline cells are retained as valid (not NA) for masking/rasterization:
  - `depth.cor` for 60×50 (18″)  
  - `depth.hires.cor` for 360×300 (3″)  
  The corrected hi-res depth is saved as:  
  `crm 360x300 3s 81m - salt marsh corrected`

**Resolution & proportion logic.**
- Two marsh rasters are produced at **60×50 (18″)**:  
  1) **Presence/absence** (0/1) marsh on the corrected base grid.  
  2) A **bilinear-resampled** product created by resampling the corrected hi-res (3″) presence raster down to the 18″ grid.  
- The **bilinear-resampled** raster yields **fractional values (0–1)** that approximate the **proportion of each 18″ cell covered by marsh**—the quantity used in Ecospace habitat capacity calculations.

**Visualization & QA.**
- Quick-look maps are plotted at each step: base/hi-res depth, initial rasterized marsh, corrected depth, corrected marsh, and the resampled (proportion) product.  

**Exports for Ecospace.**
- Presence (0/1):  
  `saltmarsh SRWMD 2016-2017 60x50 - presence.asc`
- Proportion (0–1), via bilinear downsampling from 3″ to 18″:  
  `saltmarsh SRWMD 2016-2017 60x50 - bilinear resampled.asc`

---

## Seagrass habitat

**Code:** [`make-seagrass-maps.R`](https://github.com/holden-harris/SREM-data-synthesis/blob/main/Habitats/make-seagrass-maps.R)

The seagrass layers (all, continuous, patchy) define where seagrass habitat can occur in Ecospace and provide **proportional coverage** per cell for habitat-affinity and capacity calculations.

**Source & domain.**  
Seagrass polygons are from **SRWMD seagrass 2001** (`seagrass01.shp`). Features are cropped/used within the SREM domain (**29.05–29.40° N**, **82.95–83.30° W**) and rasterized onto the Ecospace grid.

**Inputs:**
- Depth rasters (for grid and NA mask):  
  `./Data/habitats/processed/crm-resolutions/crm 360x300 3s 81m` *(hi-res rasterization grid)*  
  `./Data/habitats/processed/crm-resolutions/crm 60x50 18s 485m` *(output grid for Ecospace)*
- Seagrass shapefile (SRWMD 2001):  
  `./Data/habitats/inputs/seagrass2001_SRWMD/seagrass01.shp`

**CRS handling.**
- The SRWMD data are provided in **Florida State Plane, North Zone, NAD83, US feet**; the script explicitly sets:  
  `+proj=lcc +lat_1=30.75 +lat_2=29.58333333333333 +lat_0=29 +lon_0=-84.5 +x_0=600000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=ft +no_defs`
- Shapefile is then **reprojected to the depth raster CRS** (`depth.3s`) using `spTransform`, ensuring alignment before rasterization.

**Class selection.**
- Uses `FLUCSDESC` to subset:  
  - **All seagrass:** `Seagrass, Continuous` ∪ `Seagrass, Patchy`  
  - **Continuous only:** `Seagrass, Continuous`  
  - **Patchy only:** `Seagrass, Patchy`

**Rasterization (presence at 3″).**
- Set `AREA = 1` on polygons to mark presence.  
- `rasterize(..., depth.3s, field = "AREA", fun = "max")` produces **0/1 presence rasters** at **3 arcsec (≈81 m)** for: all, continuous, and patchy.
- Enforce mask & zeros:  
  - Non-intersecting cells → `0`  
  - Cells where `depth.3s` is `NA` (land) → `NA` (keeps shoreline mask consistent)

**Downscaling to Ecospace grid (proportion at 18″).**
- Resample each 3″ **presence** raster to the **18″** grid (`depth.out`, 60×50) with **bilinear** interpolation:  
  `resample(x_3s, depth.out, method = "bilinear")`
- The bilinear product yields **fractional values 0–1**, which approximate the **proportion of each 18″ cell covered by seagrass** (the value used in Ecospace).

**Visualization & QA.**
- Exploratory maps: class distribution (via `spplot`), then **all/continuous/patchy** at both **3″** and **18″** resolutions.  
- A 2×3 comparison panel is saved to:  
  `./Data/habitats/plots/seagrasses-continous-patch.png`

**Exports for Ecospace.**
- **All seagrass (proportion, 18″):**  
`seagrass01 all 60x50.asc` and raster-native version  
- **Continuous (proportion, 18″):**  
`seagrass01 continuous 60x50.asc` and raster-native version  
- **Patchy (proportion, 18″):**  
`seagrass01 patchy 60x50.asc` and raster-native version  

---

## Clam aquaculture lease maps

**Code:** [make-clam-map.R](https://github.com/holden-harris/SREM-data-synthesis/blob/main/Habitats/make-clam-map.R)

The clam-lease layer defines where aquaculture planting, growth, and harvest can occur in Ecospace. The map is produced as **cell-wise coverage (proportion or %)**, allowing it to be combined with other habitat layers in capacity calculations and used to constrain fleet activity.

**Source & domain.**  
Polygons are from **FDACS Suwannee River Aquaculture Use Zones (AUZs)** (`Suwannee_River_AUZs.shp`). Features are filtered and rasterized over the SREM domain (≈ **29.05–29.40° N**, **82.95–83.30° W**) on the same grid used by the Ecospace depth layer.

Data were provided by the **Florida Department of Agriculture and Consumer Services (FDACS)**,  
available for viewing at  
<https://www.fdacs.gov/Agriculture-Industry/Aquaculture/Shellfish-Harvesting-Area-and-Aquaculture-Lease-Map>.

### **Processing steps**
1. **Load and explore data.**  The shapefile of clam leases is imported and plotted to inspect lease boundaries.

2. **Filter non-active leases.**  The *DERRICKS* parcel is removed since it is no longer an active lease.

3. **Rasterize leases.**  
   The lease polygons are rasterized onto a **60×50 grid** matching the Ecospace map.  
   - `TIC` (Tenants-in-Common code) → marks lease presence.  
   - `ACRES` → used for area verification.  
   - Output is a 0/1 **presence raster**.

4. **Aggregate and calculate coverage.**  The raster is aggregated to represent the **percentage of each cell covered by leases** (0–100%).

5. **Project and align with depth grid.**  The raster is reprojected to match the depth grid’s coordinate system and resolution. Very small coverage values (<1%) are set to `NA`.

**Outputs: ASCII raster**
`./Clam-leases-2022-09` Cell values represent **percent coverage (0–100)** of clam leases per Ecospace grid cell.

---

## References

- NOAA National Centers for Environmental Information (2023). *Coastal Relief Model.*  
- Hijmans, R. J. (2025). **terra**: Spatial Data Analysis. R package.  
- Florida Department of Environmental Protection (FDEP) Geospatial Open Data Portal.  
  <https://geodata.dep.state.fl.us/>  
- Florida Department of Agriculture and Consumer Services (FDACS). *Shellfish Harvesting Areas and Aquaculture Lease Map.*  
  <https://www.fdacs.gov/Agriculture-Industry/Aquaculture/Shellfish-Harvesting-Area-and-Aquaculture-Lease-Map>  
