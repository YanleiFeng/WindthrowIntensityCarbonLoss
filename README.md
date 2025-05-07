# ESA_CCI_2020 Windthrow Biomass and Reference Biomass Analysis

**Author:** Yanlei Feng  
**Date Organized:** May 7, 2025  
**Last Updated:** October 16, 2024  

## Overview

This repository contains JavaScript code for Google Earth Engine (GEE) to prepare and analyze windthrow (tree-fall) biomass data and reference biomass in the Amazon region. The workflow includes:
- Selecting and masking forest and windthrow areas  
- Calculating biomass from ESA CCI and Avitabile datasets  
- Computing Non-Photosynthetic Vegetation (NPV) intensity and windthrow intensity  
- Extracting environmental and meteorological variables for regression analysis  

## Data Sources

- **ESA CCI Biomass Maps** (2019, 2020)  
- **Hansen Global Forest Change** (tree cover, loss, gain)  
- **LANDSAT 8 TOA** imagery for spectral mixture analysis  
- **NASA GPM IMERG** precipitation data  
- **MODIS LST** temperature data  
- **ERA5** wind shear  
- **MCS density** dataset  
- **SoilGrids-ISRIC** soil organic carbon and nitrogen  
- **Digital Elevation Model (SRTM)** for terrain variables  
- **Additional GEE assets** for distance to forest edge, main river, HAND, forest height, and soil fertility  

## Repository Structure

```
/
├── script.js            # Main Earth Engine script
├── README.md            # This documentation file
└── LICENSE              # Project license
```

## Key Functions

### maskL8Clouds(image)
Applies cloud, cloud shadow, cirrus, and dilated cloud masks to a Landsat image.

### smaima(image)
Performs Spectral Mixture Analysis (SMA) using predefined endmembers (NPV, Green Vegetation, Shade).

### selectBands(image)
Selects Landsat bands (B1–B7) for analysis.

### kernel_70(image), kernel_1070(image)
Applies focal maximum filters with 70m and 1070m radii to generate buffers.

### calculateReferenceMask(region, only_wt_mask, fm_amz)
Computes reference area mask based on elevation and soil 5th–95th percentiles within windthrow zones.

### getBiomass(mask, region), getAviBiomass(mask, region)
Extracts mean biomass from ESA CCI or Avitabile datasets for a masked region.

### reduceBiomass(biomassLayer, region)
Reduces a biomass layer to mean values within a region.

### NPV_wt(region, wt_mask, year, img, cm)
Calculates normalized NPV, determines threshold for windthrow intensity, and returns mask and mean NPV value.

### NPV_ref(region, ref_mask, year, img)
Calculates mean NPV for reference areas.

### getVariables(mask, region)
Extracts mean environmental variables (precipitation, temperature, wind shear, MCS density, topography, soil, etc.) for regression.

## Usage

1. Clone this repository to your local machine.  
2. Open `script.js` in the Google Earth Engine Code Editor.  
3. Update any asset IDs or geometry definitions as needed.  
4. Run the script to generate a `FeatureCollection` of windthrow metrics.  
5. Export the results as a CSV via `Export.table.toDrive()`.

## Update History

- **2024-10-16**: Added variables: Distance to main river, distance to forest edge, HAND, forest height, soil fertility.  
- **2024-08-07**: Extended variable list for download.  
- **2024-07-28**: Integrated David's datasets; reorganized 2014 events.  
- **2024-07-10**: Updated windthrow mask combining Hansen and custom NPV-based masks.  
- **2024-07-09**: Added NPV calculation for windthrow intensity.  
- **2025-03-30**: Modified forest mask threshold to 50%; filtered events 2015–2020.  
- **2025-03-27**: Changed buffer boundary to 1000m; switched detection to Hansen disturbance.

## License

This project is licensed under the MIT License. See `LICENSE` for details.
