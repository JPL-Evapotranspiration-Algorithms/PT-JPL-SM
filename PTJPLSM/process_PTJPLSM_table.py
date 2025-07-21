"""
Module: process_PTJPLSM_table.py

This module provides a function to process input data for the PT-JPL-SM (Priestley-Taylor Jet Propulsion Laboratory Soil Moisture) model.
It prepares the required variables from a pandas DataFrame, handles missing or alternative column names, computes derived variables as needed, and runs the PTJPLSM model to generate output variables, which are appended to the input DataFrame.
"""
import logging

import numpy as np
import rasters as rt
from rasters import MultiPoint, WGS84

from dateutil import parser
from pandas import DataFrame

from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from PTJPL import load_Topt
from PTJPL import load_fAPARmax

from .model import PTJPLSM

logger = logging.getLogger(__name__)

# FIXME include additional inputs required by PT-JPL-SM that were not required by PT-JPL


def process_PTJPLSM_table(input_df: DataFrame) -> DataFrame:
    """
    Processes an input DataFrame to prepare all required variables for the PT-JPL-SM model,
    runs the model, and returns a DataFrame with the model outputs appended as new columns.

    Args:
        input_df (DataFrame): Input data containing all necessary columns for PT-JPL-SM.

    Returns:
        DataFrame: The input DataFrame with PT-JPL-SM model outputs added as columns.
    """
    # Extract and typecast surface temperature (ST_C) and NDVI
    ST_C = np.array(input_df.ST_C).astype(np.float64)
    NDVI = np.array(input_df.NDVI).astype(np.float64)

    # Mask NDVI values below threshold (0.06) as NaN
    NDVI = np.where(NDVI > 0.06, NDVI, np.nan).astype(np.float64)

    # Extract and typecast albedo
    albedo = np.array(input_df.albedo).astype(np.float64)
    
    # Handle air temperature column name differences (Ta_C or Ta)
    if "Ta_C" in input_df:
        Ta_C = np.array(input_df.Ta_C).astype(np.float64)
    elif "Ta" in input_df:
        Ta_C = np.array(input_df.Ta).astype(np.float64)
    else:
        raise KeyError("Input DataFrame must contain either 'Ta_C' or 'Ta' column.")

    # Extract and typecast relative humidity, soil moisture, net radiation, Topt, and fAPARmax
    RH = np.array(input_df.RH).astype(np.float64)
    soil_moisture = np.array(input_df.SM).astype(np.float64)
    Rn_Wm2 = np.array(input_df.Rn).astype(np.float64)

    if "Topt" in input_df:
        # If Topt is provided, use it directly
        Topt = np.array(input_df.Topt).astype(np.float64)
    else:
        Topt = None
    
    if "fAPARmax" in input_df:
        # If fAPARmax is provided, use it directly
        fAPARmax = np.array(input_df.fAPARmax).astype(np.float64)
    else:
        fAPARmax = None

    if "canopy_height_meters" in input_df:
        # If canopy height is provided, use it directly
        canopy_height_meters = np.array(input_df.canopy_height_meters).astype(np.float64)
    else:
        canopy_height_meters = None

    if "field_capacity" in input_df:
        # If field capacity is provided, use it directly
        field_capacity = np.array(input_df.field_capacity).astype(np.float64)
    else:
        field_capacity = None

    if "wilting_point" in input_df:
        # If wilting point is provided, use it directly
        wilting_point = np.array(input_df.wilting_point).astype(np.float64)
    else:
        wilting_point = None

    # Soil heat flux (G): use provided column if available, otherwise calculate using SEBAL method
    if "G" in input_df:
        G = np.array(input_df.G).astype(np.float64)
    else:
        G = calculate_SEBAL_soil_heat_flux(
            Rn=Rn_Wm2,
            ST_C=ST_C,
            NDVI=NDVI,
            albedo=albedo
        ).astype(np.float64)

    if "geometry" in input_df:
        # If geometry is provided, use it directly
        geometry = input_df.geometry
    elif "lat" in input_df and "lon" in input_df:
        # If latitude and longitude are provided, create a MultiPoint geometry
        # Extract latitude and longitude, and create a geometry object for spatial context
        lat = np.array(input_df.lat).astype(np.float64)
        lon = np.array(input_df.lon).astype(np.float64)
        geometry = MultiPoint(x=lon, y=lat, crs=WGS84)
    elif Topt is None or fAPARmax is None or canopy_height_meters is None or field_capacity is None or wilting_point is None:
        raise KeyError("Input DataFrame must contain either 'geometry' or both 'lat' and 'lon' columns.")
    
    if Topt is None:
        # Load Topt if not provided
        Topt = load_Topt(geometry=geometry)
    
    if fAPARmax is None:
        # Load fAPARmax if not provided
        fAPARmax = load_fAPARmax(geometry=geometry)

    # If canopy height is not provided, load it
    if canopy_height_meters is None:
        from gedi_canopy_height import load_canopy_height
        canopy_height_meters = load_canopy_height(geometry=geometry)
    
    # If field capacity is not provided, load it
    if field_capacity is None:
        from soil_capacity_wilting import load_field_capacity
        field_capacity = load_field_capacity(geometry=geometry) 

    # If wilting point is not provided, load it
    if wilting_point is None:
        from soil_capacity_wilting import load_wilting_point
        wilting_point = load_wilting_point(geometry=geometry)

    # Mask fAPARmax values of zero as NaN
    fAPARmax = np.where(fAPARmax == 0, np.nan, fAPARmax).astype(np.float64)

    # Run the PTJPLSM model with all required inputs
    results = PTJPLSM(
        geometry=geometry,
        NDVI=NDVI,
        Ta_C=Ta_C,
        RH=RH,
        soil_moisture=soil_moisture,
        Rn_Wm2=Rn_Wm2,
        Topt=Topt,
        fAPARmax=fAPARmax,
        canopy_height_meters=canopy_height_meters,
        field_capacity=field_capacity,
        wilting_point=wilting_point,
        albedo=albedo,
        G=G
    )

    # Copy the input DataFrame to avoid modifying the original
    output_df = input_df.copy()

    # Append each model output as a new column in the DataFrame
    for key, value in results.items():
        output_df[key] = value

    return output_df
