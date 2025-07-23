# PTJPLSM Function Usage Guide

The `PTJPLSM` function implements the PT-JPL-SM (Priestley-Taylor Jet Propulsion Laboratory with Soil Moisture) model for partitioning evapotranspiration (ET) into its components using remote sensing and meteorological data. This function is designed for use with gridded (raster) or array-based environmental datasets, such as those from ECOSTRESS or similar sources.

## Function Signature

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

## Parameters

- **NDVI** (`Raster` or `np.ndarray`): Normalized Difference Vegetation Index.
- **Rn_Wm2** (`Raster` or `np.ndarray`): Net radiation (W/m²).
- **geometry** (`RasterGeometry`, optional): Geospatial metadata for rasters.
- **time_UTC** (`datetime`, optional): Observation time.
- **hour_of_day**, **day_of_year** (`np.ndarray`, optional): Time information.
- **GEOS5FP_connection** (`GEOS5FP`, optional): Meteorology data connection.
- **ST_C** (`Raster` or `np.ndarray`, optional): Surface temperature (°C).
- **emissivity** (`Raster` or `np.ndarray`, optional): Surface emissivity.
- **albedo** (`Raster` or `np.ndarray`, optional): Surface albedo.
- **G** (`Raster` or `np.ndarray`, optional): Soil heat flux.
- **Ta_C** (`Raster` or `np.ndarray`, optional): Air temperature (°C).
- **RH** (`Raster` or `np.ndarray`, optional): Relative humidity (0-1).
- **soil_moisture** (`Raster` or `np.ndarray`, optional): Soil moisture.
- **field_capacity**, **wilting_point** (`Raster` or `np.ndarray`, optional): Soil properties.
- **Topt** (`Raster` or `np.ndarray`, optional): Optimal plant temperature.
- **fAPARmax** (`Raster` or `np.ndarray`, optional): Maximum fAPAR.
- **canopy_height_meters** (`Raster` or `np.ndarray`, optional): Canopy height.
- **delta_Pa**, **gamma_Pa**, **epsilon**, **beta_Pa**, **PT_alpha**: Model parameters.
- **field_capacity_scale**, **minimum_Topt**: Soil and plant parameter constraints.
- **field_capacity_directory**, **wilting_point_directory**, **canopy_height_directory**: Data directories.
- **floor_Topt** (`bool`): Whether to floor Topt to Ta_C.
- **resampling** (`str`): Resampling method for raster data.

## Returns

A dictionary with the following keys (each value is a `Raster` or `np.ndarray`):

- `'G'`: Soil heat flux
- `'Rn_soil'`: Net radiation of the soil
- `'LE_soil'`: Soil evaporation
- `'Rn_canopy'`: Net radiation of the canopy
- `'PET'`: Potential evapotranspiration
- `'LE_canopy'`: Canopy transpiration
- `'LE_interception'`: Interception evaporation
- `'LE'`: Total instantaneous evapotranspiration (constrained between 0 and PET)

## Example Usage

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

## Notes

- If some inputs (e.g., `Ta_C`, `RH`, `soil_moisture`, `field_capacity`, `wilting_point`, `canopy_height_meters`, `Topt`, `fAPARmax`) are not provided, the function will attempt to load or compute them using the provided `geometry` and `time_UTC`.
- All input rasters/arrays should be spatially aligned and have the same shape.
- The function is robust to missing optional parameters, but key variables (e.g., `NDVI`, `Rn_Wm2`, `Ta_C`, `RH`, `soil_moisture`) must be provided or derivable.

## References

- Purdy et al. (2018), "PT-JPL-SM: A PT-JPL model variant incorporating soil moisture stress for improved global evapotranspiration partitioning."
