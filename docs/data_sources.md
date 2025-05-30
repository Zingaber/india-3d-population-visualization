# Data Sources Documentation

## Primary Dataset (Not included due to size limits)

### India Population Density Grid 2020
- **File**: `ind_pd_2020_1km_ASCII_XYZ.csv`
- **Size**: 225.35 MB
- **Records**: ~2.3 million grid points
- **Format**: X (longitude), Y (latitude), Z (population density)
- **Resolution**: 1km grid cells
- **Coverage**: Complete India territory
- **Source**: [Provide original data source URL]

### Administrative Boundaries
- **Source**: GADM v4.1 (Global Administrative Areas)
- **Levels**: Country, State, District, Sub-district boundaries
- **Formats**: Shapefile, GeoPackage

### Population Centers
- **Source**: OpenStreetMap via Humanitarian OSM Team
- **Formats**: GeoJSON, KML, Shapefile, GeoPackage
- **Size**: 103-290 MB files

### Global Context Data
- **Source**: LandScan Global 2023
- **Format**: TIFF raster files
- **Purpose**: Population validation and global context

## How to Obtain Data

1. **Population Grid**: Download from [original source]
2. **GADM Boundaries**: Download from https://gadm.org/
3. **OSM Population Data**: Download from Humanitarian Data Exchange
4. **LandScan**: Download from Oak Ridge National Laboratory

## Data Processing Pipeline

1. Place downloaded data files in `Download/` directory
2. Run `source("R/My3DIndiaMap.R.r")` to process data
3. Outputs will be generated in `output/` directory

## Repository Contents

This repository includes:
- ✅ **Source code**: All R scripts for processing and visualization
- ✅ **Generated outputs**: 18 high-quality 3D visualizations
- ✅ **Documentation**: Complete project documentation
- ✅ **Project structure**: Professional organization
- ❌ **Raw data**: Excluded due to GitHub size limits (available separately)

## File Size Summary
- **Total project with data**: ~500MB
- **GitHub repository**: ~20MB (code + outputs only)
- **Data files**: ~480MB (download separately)
