# Priestley-Taylor Jet Propulsion Laboratory Soil Mositure (PT-JPL-SM) Evapotranspiration Model Python Implementation

[![CI](https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL-SM/actions/workflows/ci.yml/badge.svg)](https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL-SM/actions/workflows/ci.yml)

This software package is a Python implementation of the Priestley-Taylor Jet Propulsion Laboratory Soil Moisture (PT-JP-SM) model of evapotranspiration. It was re-implemented in Python by Gregory Halverson at Jet Propulsion Laboratory based on Python code developed by AJ Purdy. 

The original [PT-JPL](https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL) model was implemented in MATLAB by Joshua Fisher and re-implemented in Python by Gregory Halverson. The PT-JPL model was designed for processing remote sensing data. It has the ability to partition latent heat flux into canopy transpiration, interception, and soil evaporation. Purdy et al., 2018 incorporated additional constraints from soil water availability on soil evaporation. Additional controls on transpiration are driven by soil water availability and canopy height include a weighting scheme to shift control on transpiration rates from soil water availability to atmospheric demand based on aridity.
 
The software was developed as part of a research grant by the NASA Research Opportunities in Space and Earth Sciences (ROSES) program. It was designed for use by the Ecosystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) mission as a precursor for the Surface Biology and Geology (SBG) mission. However, it may also be useful for general remote sensing and GIS projects in Python. This package can be utilized for remote sensing research in Jupyter notebooks and deployed for operations in data processing pipelines.
 
The software is being released according to the SPD-41 open-science requirements of NASA-funded ROSES projects.

Gregory H. Halverson (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
Lead developer<br>
NASA Jet Propulsion Laboratory 329G

Adam J. Purdy (he/him)<br>
[adpurdy@csumb.edu](mailto:adpurdy@csumb.edu)<br>
Algorithm inventor<br>
California State University Monterey Bay

Joshua B. Fisher (he/him)<br>
[jbfisher@chapman.edu](mailto:jbfisher@chapman.edu)<br>
Algorithm inventor<br>
Chapman University

Margaret C. Johnson (she/her)<br>
[maggie.johnson@jpl.nasa.gov](mailto:maggie.johnson@jpl.nasa.gov)<br>
Sensitivity analysis<br>
NASA Jet Propulsion Laboratory 398L
 
Claire Villanueva-Weeks (she/her)<br>
[claire.s.villanueva-weeks@jpl.nasa.gov](mailto:claire.s.villanueva-weeks@jpl.nasa.gov)<br>
Code maintenance<br>
NASA Jet Propulsion Laboratory 329G
 
## Installation

### PyPi Deployment

Install the `PTJPLSM` package using pip:

```
pip install PTJPLSM
```

### GitHub Development

For development, clone this repository and install locally:

```
git clone https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL-SM.git
cd PT-JPL-SM
make environment
mamba activate PTJPLSM
make install
```

## Usage

### Processing Rasters (Gridded/Array Data)

You can process gridded (raster) or array-based environmental datasets using the PT-JPL-SM model. This workflow is suitable for remote sensing data such as ECOSTRESS or similar sources, and allows for partitioning evapotranspiration (ET) into its components.

#### Function Signature

```python
from PTJPLSM import PTJPLSM

results = PTJPLSM(
     NDVI,
     Rn_Wm2,
     geometry=None,
     time_UTC=None,
     hour_of_day=None,
     day_of_year=None,
     GEOS5FP_connection=None,
     ST_C=None,
     emissivity=None,
     albedo=None,
     G=None,
     Ta_C=None,
     RH=None,
     soil_moisture=None,
     field_capacity=None,
     wilting_point=None,
     Topt=None,
     fAPARmax=None,
     canopy_height_meters=None,
     delta_Pa=None,
     gamma_Pa=GAMMA_PA,
     epsilon=None,
     beta_Pa=BETA_PA,
     PT_alpha=PT_ALPHA,
     field_capacity_scale=FIELD_CAPACITY_SCALE,
     minimum_Topt=MINIMUM_TOPT,
     field_capacity_directory=SOIL_CAPACITY_DIRECTORY,
     wilting_point_directory=SOIL_CAPACITY_DIRECTORY,
     canopy_height_directory=GEDI_DOWNLOAD_DIRECTORY,
     floor_Topt=FLOOR_TOPT,
     resampling=RESAMPLING
)
```

#### Required and Optional Inputs

- **Required:**
  - `NDVI`: Normalized Difference Vegetation Index (Raster or np.ndarray)
  - `Rn_Wm2`: Net radiation (W/m²) (Raster or np.ndarray)

- **Optional:**
  - `geometry`: Geospatial metadata for rasters
  - `time_UTC`: Observation time
  - `hour_of_day`, `day_of_year`: Time information
  - `ST_C`: Surface temperature (°C)
  - `albedo`: Surface albedo
  - `G`: Soil heat flux
  - `Ta_C`: Air temperature (°C)
  - `RH`: Relative humidity (0-1)
  - `soil_moisture`: Soil moisture
  - `field_capacity`, `wilting_point`: Soil properties
  - `Topt`: Optimal plant temperature
  - `fAPARmax`: Maximum fAPAR
  - `canopy_height_meters`: Canopy height

If some inputs (e.g., `Ta_C`, `RH`, `soil_moisture`, `field_capacity`, `wilting_point`, `canopy_height_meters`, `Topt`, `fAPARmax`) are not provided, the function will attempt to load or compute them using the provided `geometry` and `time_UTC`.

All input rasters/arrays should be spatially aligned and have the same shape.

#### Example Workflow

Suppose you have loaded ECOSTRESS or similar remote sensing data as rasters or arrays:

```python
from PTJPLSM import PTJPLSM

# Example: Load your input data (replace with your actual data loading code)
geometry = ...      # RasterGeometry object
time_UTC = ...      # datetime object
NDVI = ...          # NDVI raster or array
Ta_C = ...          # Air temperature raster or array
RH = ...            # Relative humidity raster or array
Rn = ...            # Net radiation raster or array
ST_C = ...          # Surface temperature raster or array
albedo = ...        # Albedo raster or array

# Run the PTJPLSM model
results = PTJPLSM(
     geometry=geometry,
     time_UTC=time_UTC,
     NDVI=NDVI,
     Ta_C=Ta_C,
     RH=RH,
     Rn_Wm2=Rn,
     ST_C=ST_C,
     albedo=albedo
)

# Access the total latent heat flux (evapotranspiration)
LE = results["LE"]

# Optionally, set a colormap and export to GeoTIFF
from ECOv002_granules import ET_COLORMAP
LE.cmap = ET_COLORMAP
LE.to_geotiff("example_LE.tif")
```

#### Output Keys

The returned dictionary contains the following keys (each value is a Raster or np.ndarray):

- `G`: Soil heat flux
- `Rn_soil`: Net radiation of the soil
- `LE_soil`: Soil evaporation
- `Rn_canopy`: Net radiation of the canopy
- `PET`: Potential evapotranspiration
- `LE_canopy`: Canopy transpiration
- `LE_interception`: Interception evaporation
- `LE`: Total instantaneous evapotranspiration (constrained between 0 and PET)

#### Notes
- The function is robust to missing optional parameters, but key variables (e.g., `NDVI`, `Rn_Wm2`, `Ta_C`, `RH`, `soil_moisture`) must be provided or derivable.
- All input rasters/arrays should be spatially aligned and have the same shape.

---

### Processing Tables (Batch/Site-Level Data)

You can process tabular data (e.g., site-level or point measurements, or extracted pixel values) using the PT-JPL-SM model. This workflow is suitable for batch processing, sensitivity analysis, and is compatible with ECOSTRESS Cal-Val or similar datasets.

#### Required and Optional Columns

**Required columns:**
- `NDVI`: Normalized Difference Vegetation Index
- `ST_C`: Surface temperature (°C)
- `albedo`: Surface albedo
- `Ta_C` or `Ta`: Air temperature (°C)
- `RH`: Relative humidity (0-1)
- `SM`: Soil moisture
- `Rn`: Net radiation (W/m²) *(can be computed with `verma_net_radiation_table` if not present)*

**Optional columns (will be loaded if missing):**
- `Topt`: Optimal plant temperature
- `fAPARmax`: Maximum fAPAR
- `canopy_height_meters`: Canopy height
- `field_capacity`: Soil field capacity
- `wilting_point`: Soil wilting point
- `G`: Soil heat flux *(will be calculated if missing)*
- `geometry`: Geometry object
- `lat`, `lon`: Latitude and longitude *(used to construct geometry if needed)*

#### Typical Workflow

1. **Prepare your DataFrame**
    - Ensure your data includes the required columns listed above.
    - If you do not have `Rn`, you can compute it using the `verma_net_radiation_table` function.

2. **Process the Table**
    - Use `process_PTJPLSM_table` to process your DataFrame and run the PT-JPL-SM model.

3. **Analyze the Output**
    - The output DataFrame will include the original columns plus new columns for each model output.

#### Example

Suppose you have a CSV file with columns: NDVI, ST_C, albedo, Ta_C, RH, SM, Rn, lat, lon

```python
import pandas as pd
from PTJPLSM.process_PTJPLSM_table import process_PTJPLSM_table

# Load your data
df = pd.read_csv('my_input_data.csv')

# (Optional) Compute net radiation if not present
# from verma_net_radiation import verma_net_radiation_table
# df = verma_net_radiation_table(df)

# Process the table and run the PT-JPL-SM model
output_df = process_PTJPLSM_table(df)

# The output DataFrame will have new columns: 'G', 'Rn_soil', 'LE_soil', 'Rn_canopy', 'PET',
# 'LE_canopy', 'LE_interception', 'LE' in addition to the original columns.
print(output_df.head())
```

#### Output Columns

- `G`: Soil heat flux
- `Rn_soil`: Net radiation of the soil
- `LE_soil`: Soil evaporation
- `Rn_canopy`: Net radiation of the canopy
- `PET`: Potential evapotranspiration
- `LE_canopy`: Canopy transpiration
- `LE_interception`: Interception evaporation
- `LE`: Total instantaneous evapotranspiration (constrained between 0 and PET)

#### Notes
- If any required columns are missing, a KeyError will be raised.
- If geometry is not provided, latitude and longitude columns are required to construct spatial context.
- All input columns should be numeric and of compatible shape.
- This function is suitable for batch-processing site-level or point data tables for ET partitioning and for use in sensitivity analysis workflows (see the PTJPLSM Sensitivity notebook for an example).

## Testing

Run the unit tests using pytest:

```
make test
```

## Sensitivity

The PT-JPL-SM and [PT-JPL](https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL) models, when processing on top of pre-computed net radiation, do not take surface temperature directly into account on their own. When coupled with a net radiation model that is sensitive to surface temperature, such as [verma-net-radiation](https://github.com/JPL-Evapotranspiration-Algorithms/verma-net-radiation), PT-JPL-SM exhibits moderate sensitivity to surface temperature, with an average percent change in latent heat flux of 20%.

# ![PT-JPL-SM Latent Heat Flux Sensitivity Magnitude](notebooks/PT-JPL-SM%20Latent%20Heat%20Flux%20Sensitivity%20Magnitude.jpeg)

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file.

## References

Purdy, A. J., Fisher, J. B., Goulden, M. L., Colliander, A., Halverson, G. H., Tu, K., & Famiglietti, J. S. (2018). "SMAP soil moisture improves global evapotranspiration." Remote Sensing of Environment, 219, 1-14. https://doi.org/10.1016/j.rse.2018.09.023

Fisher, J. B., Tu, K. P., & Baldocchi, D. D. (2008). "Global estimates of the land–atmosphere water flux based on monthly AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites." Remote Sensing of Environment, 112(3), 901-919. https://doi.org/10.1016/j.rse.2007.06.025
