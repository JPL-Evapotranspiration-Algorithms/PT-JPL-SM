"""
Module: generate_PTJPLSM_inputs.py
-------------------------------
This module provides a function to generate the required input DataFrame for the PT-JPL-SM model, extending the original PT-JPL input requirements. It processes a calibration/validation DataFrame, computes additional variables such as hour of day, day of year (doy), Topt, and fAPARmax, and appends them to the DataFrame. The function is robust to missing or problematic data, logging warnings and filling with NaN as needed.
"""
import logging

import numpy as np
import rasters as rt
from dateutil import parser
from pandas import DataFrame
from sentinel_tiles import sentinel_tiles
from solar_apparent_time import UTC_to_solar

from PTJPL import load_Topt
from PTJPL import load_fAPARmax

from .model import PTJPLSM

logger = logging.getLogger(__name__)

# FIXME include additional inputs required by PT-JPL-SM that were not required by PT-JPL

def generate_PTJPLSM_inputs(PTJPL_inputs_from_calval_df: DataFrame) -> DataFrame:
    """
    Generate a DataFrame with all required inputs for the PT-JPL-SM model.

    Parameters
    ----------
    PTJPL_inputs_from_calval_df : pandas.DataFrame
        DataFrame containing the columns: tower, lat, lon, time_UTC, albedo, elevation_km

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the original columns plus:
        - hour_of_day: int, hour of solar time at the site
        - doy: int, day of year
        - Topt: float, optimal temperature for photosynthesis (from spatial data)
        - fAPARmax: float, maximum fraction of absorbed photosynthetically active radiation (from spatial data)
        Additional columns may be added as required by the PT-JPL-SM model.

    Notes
    -----
    - This function is robust to missing or problematic spatial data; missing values are filled with np.nan.
    - The function logs progress and warnings for traceability.
    - If columns already exist in the input DataFrame, they are not overwritten.
    """
    # Copy input DataFrame to avoid modifying the original
    PTJPL_inputs_df = PTJPL_inputs_from_calval_df.copy()

    # Prepare lists to collect computed values
    hour_of_day = []
    doy = []
    Topt = []
    fAPARmax = []

    # Iterate over each row to compute additional variables
    for i, input_row in PTJPL_inputs_from_calval_df.iterrows():
        tower = input_row.tower
        lat = input_row.lat
        lon = input_row.lon
        time_UTC = input_row.time_UTC
        albedo = input_row.albedo
        elevation_km = input_row.elevation_km
        logger.info(f"collecting PTJPL inputs for tower {tower} lat {lat} lon {lon} time {time_UTC} UTC")
        # Parse time and convert to solar time
        time_UTC = parser.parse(str(time_UTC))
        time_solar = UTC_to_solar(time_UTC, lon)
        hour_of_day.append(time_solar.hour)
        doy.append(time_UTC.timetuple().tm_yday)
        date_UTC = time_UTC.date()
        
        try:
            # Get MGRS tile and grid for spatial data extraction
            tile = sentinel_tiles.toMGRS(lat, lon)[:5]
            tile_grid = sentinel_tiles.grid(tile=tile, cell_size=70)
        except Exception as e:
            logger.warning(e)
            Topt.append(np.nan)
            fAPARmax.append(np.nan)
            continue

        rows, cols = tile_grid.shape
        # Find the grid cell containing the point
        row, col = tile_grid.index_point(rt.Point(lon, lat))
        # Extract a 3x3 neighborhood around the point for robust statistics
        geometry = tile_grid[max(0, row - 1):min(row + 2, rows - 1),
                             max(0, col - 1):min(col + 2, cols - 1)]

        # Only compute Topt if not already present
        if not "Topt" in PTJPL_inputs_df.columns:
            try:
                logger.info("generating Topt")
                Topt_value = np.nanmedian(load_Topt(geometry=geometry))
                print(f"Topt: {Topt_value}")
                Topt.append(Topt_value)
            except Exception as e:
                Topt.append(np.nan)
                logger.exception(e)
        # Only compute fAPARmax if not already present
        if not "fAPARmax" in PTJPL_inputs_df.columns:
            try:
                logger.info("generating fAPARmax")
                fAPARmax_value = np.nanmedian(load_fAPARmax(geometry=geometry))
                print(f"fAPARmax: {fAPARmax_value}")
                fAPARmax.append(fAPARmax_value)
            except Exception as e:
                fAPARmax.append(np.nan)
                logger.exception(e)
    
    # Add computed columns to DataFrame if not already present
    if not "hour_of_day" in PTJPL_inputs_df.columns:
        PTJPL_inputs_df["hour_of_day"] = hour_of_day

    if not "doy" in PTJPL_inputs_df.columns:
        PTJPL_inputs_df["doy"] = doy
    
    if not "Topt" in PTJPL_inputs_df.columns:
        PTJPL_inputs_df["Topt"] = Topt
    
    if not "fAPARmax" in PTJPL_inputs_df.columns:
        PTJPL_inputs_df["fAPARmax"] = fAPARmax
    
    # Rename temperature column if needed for model compatibility
    if "Ta" in PTJPL_inputs_df and "Ta_C" not in PTJPL_inputs_df:
        PTJPL_inputs_df.rename({"Ta": "Ta_C"}, inplace=True)
    
    return PTJPL_inputs_df
