# Coastal Relief Model Bathymetry Maps at Multiple Resolutions

**Author:** Holden Harris  
**Date:** `r Sys.Date()`

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
