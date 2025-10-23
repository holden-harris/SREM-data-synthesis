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

