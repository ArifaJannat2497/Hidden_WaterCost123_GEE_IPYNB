# The hidden water cost of rice intensification in Bangladesh

Code and data for the manuscript by Md. Monirul Islam, Arifa Jannat and Kentaka Aruga.

## Contents

| Folder / file | Description |
|---|---|
| `gee/01_Boro_intensity_panel.js` | Google Earth Engine: annual Boro maps, stable agricultural mask, district Boro intensity 2001–2025, multi-year intensity raster (Fig. 1) |
| `gee/02_productivity_water_trajectories.js` | Google Earth Engine: 2024 Boro footprint, ET- and LST-based productivity–water trajectory classes, class areas and diagnostics (Fig. 5; Supplementary Tables 10–11) |
| `gee/03_Boro_accuracy_assessment.js` | Google Earth Engine: stratified sampling, labelling tool and area-weighted accuracy estimators for the Boro map |
| `python/01_panel_regressions.py` | Level, first-difference and elasticity models with HC1, Driscoll–Kraay and wild-cluster bootstrap inference; leave-one-district-out and district-specific estimates (Table 1; Supplementary Tables 4–6) |
| `python/02_trend_changepoint_tests.py` | Mann–Kendall, Hamed–Rao, Sen's slope and Pettitt tests (Supplementary Table 3); run after `01` |
| `python/03_sensor_transition_robustness.py` | Robustness to the Landsat 7–8 sensor transition (Supplementary Table 12) |
| `data/Boro_district_season_panel_2001_2025.csv` | District-season panel (4 districts × 25 Boro seasons) used by the Python scripts |
| `notebooks/` | Original Google Colab notebooks used during the analysis |

## Data sources

All inputs are openly available in the Google Earth Engine catalogue: Landsat 5/7/8/9 Collection 2 Level-2, ERA5-Land monthly aggregates, CHIRPS daily rainfall and FAO GAUL level 2 boundaries.

## Running the Python analysis

```
pip install pandas numpy statsmodels linearmodels pymannkendall scipy
cd python
python 01_panel_regressions.py
python 02_trend_changepoint_tests.py
python 03_sensor_transition_robustness.py
```

Outputs (CSV tables) are written to the `python/` folder.

The Earth Engine scripts run in the Earth Engine Code Editor (https://code.earthengine.google.com); set the export folder or asset path at the top of each script.
