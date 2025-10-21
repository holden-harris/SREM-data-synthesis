# Make static habitat layers for Ecospace

The primary update from **SREMv1** to **SREMv2** is that SREMv2 simulates *spatial–temporal* dynamics to represent habitats, the movement of species and fishers, and spatially explicit water conditions and primary production, with a specific emphasis on nutrient, seagrass, and fish dynamics.

Static habitat layers were developed for **depth**, **seagrass**, **oyster reefs**, **salt marsh/mangrove**, and **clam lease areas** (Figure 1E). These layers constrain the spatial extent of key habitats in the Ecospace model domain and are used to parameterize species–habitat affinities, dispersal rates, and environmental forcing.

---

## 1. Depth map

The Ecospace basemap was derived from the **NOAA Coastal Relief Model (CRM)**  
(*NOAA National Centers for Environmental Information, 2023*).  
Bathymetry was cropped to the SREM spatial domain extending from  
**29.05° N to 29.40° N** and **82.95° W to 83.30° W**.

The CRM raster (3 arc-second ≈ 90 m) was aggregated to produce a **60 × 50-cell grid** with **18 arc-second (0.005°)** square cells (≈ 485 m × 485 m).  
Depths were averaged within each grid cell using bilinear resampling and written to the `depth-map.csv` and `depth-map.tif` files used by Ecospace.

**Code:** [`make-depth-map.R`](https://github.com/holden-harris/SREM-project/blob/main/1-Prepare-Ecospace-Inputs/habitats/make-depth-map.R)

---

## 2. Salt marsh / mangrove habitat

Salt-marsh and mangrove polygons were obtained from the **Suwannee River Water Management District (SRWMD)** 2016–2017 land-use shapefiles distributed via the  
**Florida Department of Environmental Protection (FDEP) Geospatial Open Data API**  
<https://geodata.dep.state.fl.us/>.

These data were cropped to the SREM domain, re-projected to WGS 84, and downscaled to the 60 × 50 Ecospace grid.  
Each grid cell was assigned a value equal to the **proportion of cell area** covered by marsh or mangrove polygons using the `terra::rasterize()` function (Hijmans 2025).

**Code:** `make-saltmarsh-rasters.R`

---

## 3. Seagrass habitat

Seagrass polygons were obtained from the same SRWMD/FDEP 2016–2017 dataset and processed in an identical manner.  
Because seagrass distribution is highly dynamic, these static maps represent the **maximum observed extent** during the mapping period.  
Dynamic seagrass feedbacks are represented separately in Ecospace through mediation functions (see *Chagaris et al., SREM draft*).

**Code:** `make-seagrass-maps.R`

---

## 4. Clam aquaculture leases

Clam-lease polygons were provided by the **Florida Department of Agriculture and Consumer Services (FDACS)**,  
available for viewing at  
<https://www.fdacs.gov/Agriculture-Industry/Aquaculture/Shellfish-Harvesting-Area-and-Aquaculture-Lease-Map>.

Polygons were cropped to the Suwannee Estuary, merged by lease type, and rasterized to the Ecospace grid.  
These areas define where clam aquaculture biomass can occur and where the clam-harvest fleet operates in Ecospace.

**Code:** `make-clam-map.R`

---

## References

- NOAA National Centers for Environmental Information (2023). *Coastal Relief Model.*  
- Hijmans, R. J. (2025). **terra**: Spatial Data Analysis. R package.  
- Florida Department of Environmental Protection (FDEP) Geospatial Open Data Portal.  
  <https://geodata.dep.state.fl.us/>  
- Florida Department of Agriculture and Consumer Services (FDACS). *Shellfish Harvesting Areas and Aquaculture Lease Map.*  
  <https://www.fdacs.gov/Agriculture-Industry/Aquaculture/Shellfish-Harvesting-Area-and-Aquaculture-Lease-Map>  



