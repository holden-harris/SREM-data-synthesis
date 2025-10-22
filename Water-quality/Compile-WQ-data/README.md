# Water Quality Data Synthesis (Suwannee Estuary)

**Purpose.**  
Compile, harmonize, and summarize water-quality observations from seven monitoring sources into daily and monthly datasets (with river‐flow covariates) for use in Ecospace spatial–temporal driver maps and model fitting.

---

## Data sources

1. **FDACS — Florida Department of Agriculture and Consumer Services (shellfish sanitation water quality; CK/SS/HB areas).**  
   The Division of Aquaculture within FDACS classifies and monitors shellfish-harvesting areas under Florida’s implementation of the National Shellfish Sanitation Program. Monitoring includes bacteriological (e.g., fecal coliform) sampling, area closures, and certification of processing plants. [FDACS Shellfish Aquaculture & Certification](https://www.fdacs.gov/Agriculture-Industry/Aquaculture/Shellfish) :contentReference[oaicite:0]{index=0}

2. **FIM — Fisheries Independent Monitoring (FWC/FWRI physical variables).**  
   The FWC Fish & Wildlife Research Institute’s FIM program uses a stratified-random sampling design (inshore and offshore) to monitor fish and invertebrate populations across Florida estuaries. At each sampling site physical/chemical data (temperature, salinity, dissolved oxygen, etc.) are also recorded, providing a long-term environmental time-series. [FWC FIM Program](https://myfwc.com/research/saltwater/fim/) :contentReference[oaicite:1]{index=1}

3. **Project COAST long-term monitoring (Prof. Tom Frazer et al.).**  
   Project COAST, led by Prof. Tom Frazer at UF, operated a long-term coastal water-quality monitoring network (late 1990s onward) in the Nature Coast / Big Bend region, tracking nutrient concentrations, chlorophyll, salinity and temperature to assess eutrophication risk and baseline conditions.

4. **Lakewatch — UF/IFAS Florida LAKEWATCH (field YSI + lab parameters).**  
   LAKEWATCH is a volunteer-based program coordinated by UF/IFAS that has collected monthly physical and chemical data from lakes, springs and estuarine sites for over 30 years. Parameters include temperature, Secchi depth, total nitrogen, total phosphorus and chlorophyll-a. [UF IFAS LAKEWATCH](https://lakewatch.ifas.ufl.edu/about-us/what-is-florida-lakewatch/)  

5. **LCR — UF Lone Cabbage Reef restoration monitoring (fixed sites).**  
   The Lone Cabbage Reef restoration project (UF/IFAS Nature Coast Biological Station) established fixed monitoring stations for oyster reef habitat restoration in the Suwannee Estuary. The dataset includes environmental measurements (salinity, temperature) at fixed sites through the restoration period.

6. **VBuoys — USF Optical Oceanography “Virtual Buoy” remote sensing stations (temperature, optical data).**  
   The Virtual Buoy System (VBS) from the Optical Oceanography Lab at USF (led by Dr. Chuanmin Hu) leverages satellite remote-sensing algorithms to produce near-real-time “virtual station” time series for coastal water temperature, optical water quality metrics and water-colour parameters. [USF Virtual Buoy System](https://optics.marine.usf.edu/projects/vbs.html) :contentReference[oaicite:2]{index=2}

7. **USGS Flow Rates — U.S. Geological Survey Suwannee River (Wilcox gage 02323500).**  
   Daily discharge data from the USGS gage near Wilcox provides long-term hydrologic context for the Suwannee River system. Data are accessible via USGS Water Data for the Nation.  
   [USGS Water Data – Station 02323500](https://waterdata.usgs.gov/nwis/uv?site_no=02323500)  


**Snapshot of source coverage (post-filters, 1997–2020):**

| Source     | n_samples | Years | Months | Days | Sites | Salinity | Temperature | FC  | TNP  |
|:------------|----------:|------:|-------:|-----:|------:|---------:|-------------:|----:|----:|
| **FDACS**   | 49,677    | 24    | 288    | 1,761| 142   | 402      | 302          | 103 | 0  |
| **FIM**     | 14,151    | 24    | 287    | 2,133| 14,151| 361      | 301          | 0   | 0  |
| **Frazer**  | 2,160     | 19    | 215    | 226  | 10    | 1,238    | 1,240        | 0   | 157|
| **Lakewatch** | 274     | 4     | 41     | 50   | 6     | 101      | 78           | 0   | 113|
| **LCR**     | 10,365    | 4     | 41     | 1,193| 10    | 10,364   | 9,704        | 0   | 0  |
| **VBuoys**  | 7,025     | 21    | 251    | 251  | 28    | 0        | 6,885        | 0   | 0  |

---

## Inputs (paths used in code)

- `./Data/water-quality/inputs/fdacs_wq_CK-SS-HB.csv`
- `./Data/water-quality/inputs/lcr_wq_total.csv`  
  `./Data/water-quality/inputs/lcr_sites_lat_long.csv`
- `./Data/water-quality/inputs/lab.csv` *(Lakewatch combined YSI + lab dataset)*
- `./Data/water-quality/inputs/frazer_suwannee_97-15.csv`
- `./Data/water-quality/inputs/FIM CK full physical dataset 2020.csv`
- `./Data/water-quality/inputs/virtual_buoy_temps_compiled.csv`

**River flow (USGS Wilcox gauge 02323500):** fetched via `waterData::importDVs()`.

---

## Outputs (written by the script)

- **Merged, point‐level (with daily flow MAs):**  
  `./Data/water-quality/processed/Sptatial-temp_Phys_all-measurements-with-Flow.csv`
- **Daily site averages (+ daily flow MAs):**  
  `./Data/water-quality/processed/Spatial-temp_Phys-Flow_xDay.csv`
- **Monthly site averages (+ monthly mean flow):**  
  `./Data/water-quality/processed/Spatial-temp_Phys-Flow_xMonth.csv`
- **Flow with moving averages & lags (1–45 days):**  
  `./Data/water-quality/processed/Flow-with-avgs-lags.csv`
- **Lakewatch joined with LCR (to fill salinity/temperature gaps):**  
  `./Data/water-quality/processed/lakewatch_sal-temp_joined.csv`

---

## 🧭 Processing Notes

### **1️⃣ Data Import and Cleaning**
- Each dataset is read from its respective CSV file in `./Data/water-quality/inputs/`.
- Columns are standardized to lowercase, and date/time fields are parsed and formatted (`Date`, `Year`, `Month`, `Day`, `YM`).
- Data are filtered to remove inactive or out-of-range records (e.g., `FDACS status == "ACTIVE"`, `LCR in_service == 1`, `year >= 1997`).
- Latitude and longitude columns are standardized for spatial joins.
- For **Lakewatch**, missing salinity and temperature values are supplemented by:
  - Merging *YSI sensor* readings when available.
  - Cross-referencing *LCR* observations at matching site and hourly timestamps.

---

### **2️⃣ Harmonization**
- Each data source is transformed into a unified table with the following columns:  
  `Source, Date, Time, Year, Month, YM, Lat, Long, Site_ID, Salinity, Temperature, DO, FC, Chl, TN, TP`.
- Derived fields are added:
  - `TNP = TN + TP`
  - `YM = "YYYY-MM"`
- Records are concatenated using `rbind()` to form a single harmonized dataset `physcomp`.

---

### **3️⃣ River Flow Data**
- Daily discharge data are imported from the **USGS Wilcox Gage (02323500)** using the `waterData` package.
- Flow variables are summarized using rolling statistics with `zoo::rollapply()`:
  - Moving averages: `ma01`, `ma05`, `ma15`, `ma30`
  - Rolling standard deviations: `sd05`, `sd15`, `sd30`
- Lagged flow variables (`lag0`–`lag45`) are generated to support time-lag analyses.

---

### **4️⃣ Data Integration**
- The harmonized water-quality dataset (`physcomp`) is merged with daily flow metrics by `Date` to create `physcomp2`.
- The merged table includes both **raw measurements** and **flow covariates**, enabling linked hydrologic–water quality analyses.

---

### **5️⃣ Aggregation**
- **Daily Site Means (`physcomp3`)**  
  Grouped by `Source`, `Site_ID`, and `Date` to calculate:
  - Mean, min, and max of `Salinity` and `Temperature`
  - Mean of `DO`, `FC`, `Chl`, `TN`, `TP`, and `TNP`
  - Sample counts (`n_samp`)
- **Monthly Site Means (`physcomp4`)**  
  Further aggregated by `YM` to produce average conditions per site per month.
  - Merged with monthly-mean flow.
  - Rounded values for reporting consistency.

---

### **6️⃣ Output Files**
| Output File | Description |
|--------------|-------------|
| `Spatial-temp_Phys_all-measurements-with-Flow.csv` | All individual water-quality measurements joined with daily flow stats |
| `Spatial-temp_Phys-Flow_xDay.csv` | Daily site means with daily flow data |
| `Spatial-temp_Phys-Flow_xMonth.csv` | Monthly site means with monthly mean flow |
| `Flow-with-avgs-lags.csv` | Daily flow dataset with moving averages and lag variables |
| `lakewatch_sal-temp_joined.csv` | Lakewatch dataset with missing salinity/temperature filled from LCR data |

---

## ⚙️ Key Notes

- **Temporal Scope:**  
  All data is restricted to **1997–2020** to align with Ecospace simulation years.

- **Variable Units:**  
  - Salinity (psu), Temperature (°C)
  - DO (mg/L)
  - FC (cfu/100 mL)  
  - Chlorophyll-a (µg/L) 
  - TNP = TN + TP, TN and TP (µg/L)  
  - Flow (cubic feet per second, cfs)

- **Gap-Filling:**  
  Lakewatch data enriched with YSI and LCR observations to improve salinity and temperature completeness.

- **Quality Assurance:**  
  - Filters applied for valid salinity ranges (`0–45 psu`) and active sites only.  
  - Spot-check plots (e.g., `Temperature ~ Date` by site) used to verify data integrity visually.  

- **Flow Features:**  
  Moving averages capture smooth hydrologic trends.  
  Lag variables (up to 45 days) enabled exploring delayed water-quality responses.
