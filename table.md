# Processing Tables with PTJPLSM

This guide explains how to process tabular data (e.g., site-level or point measurements, or extracted pixel values) using the PT-JPL-SM model in Python. The workflow is suitable for batch processing, sensitivity analysis, and is compatible with ECOSTRESS Cal-Val or similar datasets.

## Overview

The function `process_PTJPLSM_table` takes a pandas DataFrame as input, prepares all required variables, handles missing or alternative column names, computes derived variables as needed, and runs the PT-JPL-SM model. The model outputs are appended as new columns to the DataFrame.

## Required and Optional Columns

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

## Typical Workflow

1. **Prepare your DataFrame**
   - Ensure your data includes the required columns listed above.
   - If you do not have `Rn`, you can compute it using the `verma_net_radiation_table` function.

2. **Process the Table**
   - Use `process_PTJPLSM_table` to process your DataFrame and run the PT-JPL-SM model.

3. **Analyze the Output**
   - The output DataFrame will include the original columns plus new columns for each model output.

## Example

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

## Output Columns

- `G`: Soil heat flux
- `Rn_soil`: Net radiation of the soil
- `LE_soil`: Soil evaporation
- `Rn_canopy`: Net radiation of the canopy
- `PET`: Potential evapotranspiration
- `LE_canopy`: Canopy transpiration
- `LE_interception`: Interception evaporation
- `LE`: Total instantaneous evapotranspiration (constrained between 0 and PET)

## Notes
- If any required columns are missing, a KeyError will be raised.
- If geometry is not provided, latitude and longitude columns are required to construct spatial context.
- All input columns should be numeric and of compatible shape.
- This function is suitable for batch-processing site-level or point data tables for ET partitioning and for use in sensitivity analysis workflows (see the PTJPLSM Sensitivity notebook for an example).
