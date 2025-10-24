# Make krige maps
R workflow to **build monthly kriged maps** (and optional IDW maps) for multiple water-quality variables from synthesized station data, and writes the outputs (stacks + variogram PDFs) in the repo’s `Water-quality/Krig-and-map/out/` tree. Makes kriged maps for:
- Salinity
- Temperature
- Nutrients (TNP)
- Dissolved oxygen (DO)
- Fecal coliform (FC)

Salinity and temperature maps are used in the Ecospace model to make the predicted nutrient maps. Note that the kriged maps for TNP, DO, and FC are not currently used by the Suwannee River Ecospace Model, despite the workflow to produce them are included here. 

## Background: Spatial interpolation

Spatial interpolation techniques are essential when working with environmental monitoring data, as sampling is always conducted at discrete locations while management and modeling often require predictions at unsampled locations. **Inverse Distance Weighting (IDW)** offers a deterministic interpolation method that computes a predicted value at an unsampled location as a weighted average of nearby observations, with weights inversely proportional to distance raised to a power parameter (commonly 2). While IDW is computationally simpler and more intuitive than kriging, it does not model spatial covariance explicitly and thus lacks a formal error estimation framework. 

**Kriging** offers a statistical approach grounded in geostatistics that treats the observed phenomenon as a realization of a random function with spatial correlation. The key steps involve computing an empirical variogram, fitting a theoretical model (spherical, exponential, or stable), and then using the fitted model to derive interpolation weights that minimize prediction variance under the unbiasedness constraint (often termed “ordinary kriging”). In R, packages such as `gstat` and `automap` simplify the process by offering variogram estimation and automatic kriging in a few lines of code. 

This script performs **ordinary kriging**, implemented through the `automap::autoKrige()` function in R.  
Ordinary kriging assumes that the spatial process has an **unknown but constant mean** within the local neighborhood of interpolation and that spatial autocorrelation can be modeled through a **variogram function**.  

In each monthly iteration, `autoKrige()` automatically fits and selects among three common **variogram model types**:
- **Spherical (`"Sph"`)** – assumes correlation decreases linearly with distance until reaching a defined range, beyond which values are uncorrelated.  
- **Exponential (`"Exp"`)** – models correlation that decays exponentially with distance, producing smoother surfaces.  
- **Stable (`"Ste"`)** – a more flexible model that allows intermediate behaviors between spherical and exponential depending on its shape parameter (`kappa`).  

For each variable (Salinity, Temperature, DO, FC, TNP), the function:
1. Estimates the **empirical variogram** from observed data. These are saved in corresponding PDF documents. 
2. Fits candidate theoretical models (Sph, Exp, Ste) by weighted least squares.  
3. Selects the best-fitting model to predict values at all grid cells using ordinary kriging.  

This method provides **unbiased spatial predictions** with **minimum estimation variance**, making it particularly suitable for environmental datasets with moderate spatial autocorrelation and irregularly spaced sampling locations.

## Code Overview
1. **Load data & libraries**; sets a reproducible plotting palette.  
2. **Augment upstream coverage** by adding three small “dummy” stations near the NE corner (replicates the most upstream salinity each month; improves kriging near the boundary).  
3. **Project** observations to UTM and **set the spatial grid** by reprojecting a high-res depth raster.  
4. **Krige monthly fields** of each variable (Salinity, Temperature, DO, FC, TNP) using `automap::autoKrige`, saves:  
   - a **stack** (all months across all years) and  
   - a **Monthly_variogram_<var>.pdf** with fit/variograms.  
5. **Write outputs**, e.g.,
   `./Water-quality/Krig-and-map/out/KRG/Salinity_KRG_Jan1997-Dec2020.grd/gri`  
   and a companion `<…>_VGpars.csv` of variogram parameters by Year–Month.  
6. (Additional) **Build IDW maps** for the same variables and writes stacks to `out/IDW/`.

## Inputs
- **Station data:** `./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv`  
  (includes fields such as `Year`, `Month`, `YM`, `Lat`, `Long`, `Salinity`, `Temperature`, `DO`, `FC`, `TNP`, and flow summaries)
- **Depth grid:** `./Data/habitats/processed/crm-salt-marsh-corrected/crm 60x50 18s 485m - saltmarsh corrected.gri`  
  (used to define the kriging prediction grid after reprojection)

## Kriging loop

The loop automates the creation of monthly kriged water-quality maps for each environmental variable (Salinity, Temperature, DO, FC, TNP).
For each variable, it performs the following steps:

1. **Setup and output creation:**  
   Initializes an empty raster stack (`wq.krg.stack`) and a data frame (`vgpars.out`) to store variogram parameters.  
   Opens a multi-page PDF named `Monthly_variogram_<variable>.pdf` for diagnostic plots.

2. **Iterate through time:**  
   Loops over each year and month in the dataset, subsetting observations (`kdat`) for the current variable and removing missing values.  
   Skips months with no data.

3. **Perform kriging:**  
   Fits an ordinary kriging model using `automap::autoKrige()` with Spherical, Exponential, and Stable variogram models.  
   Produces a predicted surface on the UTM grid and appends the variogram fit and map to the open PDF.

4. **Post-process predictions:**  
   Converts the kriging output to a raster, clamps any negative values to zero, reprojects the layer to WGS84 (longitude–latitude),  
   and ensures it matches the spatial resolution and extent of the depth basemap (`depth.out`).

5. **Save outputs:**  
   Adds each month’s raster to the cumulative stack and records model parameters (model type, nugget, sill, range, kappa).  
   After all months are processed, the full raster stack and variogram parameter table are saved to:
   - `<Variable>_KRG_Jan1997-Dec2020.gri/.grd` (monthly kriged rasters)  
   - `<Variable>_KRG_Jan1997-Dec2020_VGpars.csv` (variogram metadata)  
   - `Monthly_variogram_<Variable>.pdf` (fit diagnostics)

## Inverse Distance Weighting (IDW) mapping workflow

This section generates monthly water-quality maps using **Inverse Distance Weighting (IDW)** interpolation as an alternative to kriging.  
The workflow computes spatially continuous rasters for Salinity, Temperature, Dissolved Oxygen (DO), Fecal Coliform (FC), and Total Nutrients (TNP).

1. **Define interpolation parameters:**  
   The IDW power (`idw.wt = 2`) controls how quickly influence decreases with distance.  
   Higher values make local observations more dominant (less smooth surfaces).
   
3. **Initialize empty raster stacks:**  
   One `RasterStack` is created for each variable (e.g., `sal.idw.stack`, `temp.idw.stack`, etc.) to store monthly predictions across all years.

4. **Iterate through time:**  
   The code loops through each **year** and **month**, subsetting the dataset (`wq.sub`) to include only valid measurements for the current variable.

5. **Run IDW interpolation:**  
   For each variable and month, the `gstat::idw()` function estimates values at unsampled grid cells based on the weighted average of nearby observations:

```text
           Σ [ Z(x_i) / d(x_i, x_0)^p ]
Z(x_0) =  --------------------------------
             Σ [ 1 / d(x_i, x_0)^p ]
```
   where \(d(x_i, x_0)\) is distance between points and \(p = 2\) is the inverse-distance power.

5. **Generate and plot monthly rasters:**  
   Each interpolated surface is converted to a raster layer, labeled with its month-year (e.g., `Jan1997`), and plotted for quick visualization. Each monthly raster is appended to the appropriate variable-specific stack (e.g., `sal.idw.stack` for salinity).  

---

# Nutrient Prediction Modeling and Mapping

### Overview
R code for predicting and mapping monthly nutrient concentrations (total nitrogen + phosphorus; *TNP*) across an estuarine ecosystem.  

Predictions are generated using a **linear mixed-effects model** that relates nutrient concentrations to salinity, temperature, and freshwater discharge. 

The fitted model is then applied spatially to monthly kriged predictor rasters to produce a continuous raster time series of nutrient maps from 1997–2020.

---

### Inputs
| Input File | Description |
|-------------|-------------|
| `./Data/water-quality/processed/Nutrients_withFlow.csv` | Monthly water-quality dataset containing nutrient concentrations, salinity, temperature, and a 30-day moving average of flow. |
| `./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv` | Physical data used to calculate monthly mean river flow for the study period. |
| `./Water-quality/Krig-and-map/out/KRG/Salinity_KRG_Jan1997-Dec2020.gri` | Kriged monthly salinity raster stack used as a spatial predictor. |
| `./Water-quality/Krig-and-map/out/KRG/Temperature_KRG_Jan1997-Dec2020.gri` | Kriged monthly temperature raster stack used as a spatial predictor. |

---

## Model fitting (linear mixed-effects model)
A linear mixed-effects model (LMM) was fitted to quantify relationships between nutrient concentrations and key environmental drivers.  

The response variable was **total nitrogen + phosphorus (TNP)**, and the fixed effects were **salinity (Sal; PSU)**, **temperature (Temp; °C)**, and the **log-transformed 30-day moving average of flow (logFlow; m³ s⁻¹)**.  
A random intercept for **Site_ID** was included to account for repeated observations across 16 monitoring sites.

**Model formula:**
```r
TNP ~ Sal + Temp + logFlow + (1 | Site_ID)
```

**R code:**
```r
nutri_lmm <- lme4::lmer(TNP ~ Sal + Temp + logFlow + (1|Site_ID), data = nutr)
```

### Results and explanation

The linear mixed-effects model revealed that nutrient concentrations (TNP) were strongly structured by salinity, temperature, and freshwater flow.  

| Fixed Effect | Estimate | Std. Error | t-value | Direction |
|---------------|-----------|-------------|----------|------------|
| (Intercept) | 137.49 | 89.04 | 1.54 | — |
| Sal | −18.21 | 0.85 | −21.52 | Negative |
| Temp | 7.99 | 0.73 | 10.91 | Positive |
| logFlow | 85.26 | 8.14 | 10.48 | Positive |

Across 2,520 monthly observations from 16 estuarine sites, **salinity** had the strongest negative effect on nutrient concentrations, while both **temperature** and **log-transformed river flow** had positive effects.

The model explained a large proportion of variance in nutrient concentrations (**marginal R² = 0.63**, **conditional R² = 0.68**), with moderate random variation among sites (**SD = 122 µg L⁻¹**). Residual patterns were evenly distributed, indicating that the model captured both spatial and temporal variation effectively.  

### Model diagnostics and evaluation  
Model evaluation included visual checks of residuals, assessment of variance inflation factors (VIF) for multicollinearity, and calculation of random-effect variance components. Residuals showed no strong heteroscedasticity or nonlinearity, and normal Q-Q plots confirmed approximate normality.

## Importing kriged environmental predictor rasters

Monthly kriged raster stacks of **salinity** and **temperature** were imported to serve as spatial predictor inputs for the model.  
Each raster stack consists of 288 layers (12 months × 24 years) at the same spatial resolution and projection.  

Additionally, monthly mean river discharge values were computed from the 30-day moving average of flow data and log-transformed to match the predictor variable used in the model. The resulting datasets provide spatially explicit salinity and temperature predictors for each grid cell and a uniform log-flow predictor for each month.

These layers were combined into monthly stacks and supplied to the fitted mixed-effects model to generate spatial predictions of nutrient concentrations.

## Applying the model spatially to generate predicted nutrient maps
For each month, the kriged **salinity** and **temperature** rasters and the uniform **log-flow** raster were stacked together to form a three-layer predictor dataset. This computation was performed across all raster cells using the `raster::predict()` function in R.  

Each grid cell’s predicted nutrient concentration was calculated as the linear combination of these three predictor values and their estimated coefficients from the fitted linear mixed-effects model (LMM):

**Predicted nutrient concentration for each grid cell:**
`TNP_cell = β0 + βSal * Sal_cell + βTemp * Temp_cell + βlogFlow * logFlow_month`

Where:
- `Sal_{cell}` and `Temp_{cell}` are the cell-specific values from the kriged salinity and temperature rasters,  
- `logFlow_{month}` is the spatially uniform log-transformed discharge for that month, and  
- `β₀`, `β_Sal`, `β_Temp`, and `β_logFlow` are the fixed-effect coefficients estimated from the LMM.

The resulting raster layer provided the modeled nutrient concentration for each grid cell in that month.  

These monthly rasters were sequentially appended into a multi-layer stack (`nutri_stack`), yielding a continuous time series of predicted nutrient distributions that reflect spatial gradients in salinity and temperature, modulated by monthly freshwater flow.

![Predicted monthly nutrient maps for 1997](./figures/nutri_stack_1997.png)

