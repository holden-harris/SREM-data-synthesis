# Make krige maps

This R workflow builds **monthly kriged maps** (and optional **IDW maps**) for multiple water-quality variables from synthesized station data. Outputs include raster stacks and variogram diagnostic PDFs saved to the repository’s `Water-quality/Krig-and-map/out/` directory.  
The workflow produces kriged maps for:
- Salinity  
- Temperature  
- Nutrients (Total Nitrogen + Phosphorus; *TNP*)  
- Dissolved Oxygen (DO)  
- Fecal Coliform (FC)  

The kriged salinity and temperature maps are used (1) directly by the Suwannee River Ecosystem and (2) to generate predicted nutrient (TNP) maps below. While this workflow produces kriged and IDW maps for TNP, DO, and FC, those products aren't currently used. 

---

## Background: Spatial interpolation

Spatial interpolation techniques are fundamental for environmental modeling, as field observations are collected at discrete stations, while management and modeling require spatially continuous predictions.  

**Inverse Distance Weighting (IDW)** is a deterministic interpolation technique that estimates values at unsampled locations as weighted averages of nearby observations, with weights inversely proportional to distance raised to a power parameter (commonly 2).  
IDW is computationally efficient and easy to interpret, but it does not account for spatial autocorrelation explicitly and lacks a formal error estimation framework. **Kriging**, by contrast, is a geostatistical method that models spatial dependence using a variogram and produces statistically optimal (minimum-variance, unbiased) estimates.  

This workflow performs **ordinary kriging**, which assumes the process has an **unknown but constant mean** within each local neighborhood, and models spatial covariance using theoretical variogram functions.  

For each monthly time step, the script uses `automap::autoKrige()` to:
1. Compute the **empirical variogram** from observed data (saved as monthly diagnostic PDFs).  
2. Fit three common variogram models—**Spherical (`"Sph"`)**, **Exponential (`"Exp"`)**, and **Stable (`"Ste"`)**—using weighted least squares.  
3. Select the best-fitting model and generate interpolated predictions across the spatial grid.

This approach provides **unbiased spatial predictions** with **quantified variance**, making it well-suited for environmental datasets with moderate spatial autocorrelation and irregularly spaced sampling locations.

## Code Overview
1. **Load data & libraries:** Imports required packages and sets a consistent color palette for plots.  
2. **Augment upstream coverage:** Adds three “dummy” upstream stations that replicate the lowest observed salinity each month to stabilize kriging near the northern boundary.  
3. **Project data and define grid:** Transforms station coordinates to UTM (Zone 16N, km units) and projects a depth raster (`depth.out`) to the same coordinate system to serve as the kriging prediction grid.  
4. **Kriging interpolation:**  
   Loops over all years and months for each variable, using `autoKrige()` to fit the variogram, perform kriging, and save:  
   - A **raster stack** containing all months and years.  
   - A **PDF file** with variogram and diagnostic plots (`Monthly_variogram_<Variable>.pdf`).  
5. **Write outputs**, e.g.,
   `./Water-quality/Krig-and-map/out/KRG/Salinity_KRG_Jan1997-Dec2020.grd/gri`  
   and a companion `<…>_VGpars.csv` of variogram parameters by Year–Month.  
6. (Additional) **Build IDW maps** for the same variables and writes stacks to `out/IDW/`.

### Inputs
| Input File | Description |
|-------------|-------------|
| `./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv` | Contains monthly water-quality observations including `Year`, `Month`, `YM`, `Lat`, `Long`, `Salinity`, `Temperature`, `DO`, `FC`, `TNP`, and flow metrics. |
| `./Data/habitats/processed/crm-salt-marsh-corrected/crm 60x50 18s 485m - saltmarsh corrected.gri` | High-resolution bathymetry and habitat raster used to define the kriging prediction grid and spatial extent after reprojection. |

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
   After all months are processed, the full raster stack and variogram parameter table are saved to the following. 

### Outputs
| Output File | Description |
|--------------|-------------|
| `<Variable>_KRG_Jan1997-Dec2020.gri/.grd` | Multi-layer raster stack containing monthly kriged predictions for each variable (Salinity, Temperature, DO, FC, TNP) from January 1997 to December 2020. |
| `<Variable>_KRG_Jan1997-Dec2020_VGpars.csv` | Table of fitted variogram parameters for each Year–Month, including model type, nugget, partial sill, range, and kappa values. |
| `Monthly_variogram_<Variable>.pdf` | Multi-page diagnostic PDF showing empirical and fitted variograms for each monthly kriging model, used to evaluate model fits and spatial dependence. |

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
   Like the krige maps, each interpolated surface is converted to a raster layer, labeled with its month-year (e.g., `Jan1997`), and plotted for quick visualization. Each monthly raster is appended to the appropriate variable-specific stack (e.g., `sal.idw.stack` for salinity).  

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

