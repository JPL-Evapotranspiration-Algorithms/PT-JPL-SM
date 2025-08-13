"""
PTJPLSM Model Implementation
---------------------------
This module provides the PTJPLSM function, which implements the PT-JPL-SM (Priestley-Taylor Jet Propulsion Laboratory with Soil Moisture) model for partitioning evapotranspiration into its components using remote sensing and meteorological data.

Main Function:
    PTJPLSM(...):
        Computes soil evaporation, canopy transpiration, interception evaporation, and total evapotranspiration (LE) using a variety of biophysical and meteorological inputs. Handles missing data by attempting to load or compute required variables from provided geometry and time information.

Returns:
    Dict[str, Union[Raster, np.ndarray]]: Dictionary containing the following keys:
        - 'G': Soil heat flux
        - 'Rn_soil': Net radiation of the soil
        - 'LE_soil': Soil evaporation
        - 'Rn_canopy': Net radiation of the canopy
        - 'PET': Potential evapotranspiration
        - 'LE_canopy': Canopy transpiration
        - 'LE_interception': Interception evaporation
        - 'LE': Total instantaneous evapotranspiration (constrained between 0 and PET)

References:
    - Purdy et al. (2018), "PT-JPL-SM: A PT-JPL model variant incorporating soil moisture stress for improved global evapotranspiration partitioning."
"""

from typing import Union, Dict
from pytictoc import TicToc
import warnings
import numpy as np
import pandas as pd
from datetime import datetime
import logging

import rasters as rt
from rasters import Raster, RasterGeometry
from GEOS5FP import GEOS5FP

from check_distribution import check_distribution
from soil_capacity_wilting import DEFAULT_DOWNLOAD_DIRECTORY as SOIL_CAPACITY_DIRECTORY
from soil_capacity_wilting import load_field_capacity, load_wilting_point
from gedi_canopy_height import GEDI_DOWNLOAD_DIRECTORY
from gedi_canopy_height import load_canopy_height

from carlson_leaf_area_index import carlson_leaf_area_index
from sun_angles import calculate_daylight
from verma_net_radiation import verma_net_radiation, daily_Rn_integration_verma
from daily_evapotranspiration_upscaling import daily_ET_from_instantaneous, daily_ET_from_daily_LE

from PTJPL import GAMMA_PA
from PTJPL import BETA_PA
from PTJPL import PT_ALPHA
from PTJPL import MINIMUM_TOPT
from PTJPL import FLOOR_TOPT

from PTJPL import SVP_Pa_from_Ta_C
from PTJPL import delta_Pa_from_Ta_C
from PTJPL import SAVI_from_NDVI
from PTJPL import fAPAR_from_SAVI
from PTJPL import fIPAR_from_NDVI
from PTJPL import calculate_relative_surface_wetness
from PTJPL import calculate_green_canopy_fraction
from PTJPL import calculate_plant_moisture_constraint
from PTJPL import calculate_plant_temperature_constraint
from PTJPL import calculate_soil_net_radiation
from PTJPL import calculate_interception
from PTJPL import load_Topt
from PTJPL import load_fAPARmax
from PTJPL import calculate_SEBAL_soil_heat_flux

from .constants import *
from .partitioning import (
    calculate_fREW, calculate_fTRM,
    calculate_soil_latent_heat_flux, calculate_canopy_latent_heat_flux
)

logger = logging.getLogger(__name__)

def PTJPLSM(
        NDVI: Union[Raster, np.ndarray, float],
        Rn_Wm2: Union[Raster, np.ndarray, float] = None,
        Rn_daily_Wm2: Union[Raster, np.ndarray, float] = None,
        geometry: RasterGeometry = None,
        time_UTC: datetime = None,
        hour_of_day: np.ndarray = None,
        day_of_year: np.ndarray = None,
        GEOS5FP_connection: GEOS5FP = None,
        ST_C: Union[Raster, np.ndarray] = None,
        emissivity: Union[Raster, np.ndarray, float] = None,
        albedo: Union[Raster, np.ndarray, float] = None,
        G_Wm2: Union[Raster, np.ndarray, float] = None,
        SWin_Wm2: Union[Raster, np.ndarray, float] = None,
        Ta_C: Union[Raster, np.ndarray, float] = None,
        RH: Union[Raster, np.ndarray, float] = None,
        soil_moisture: Union[Raster, np.ndarray, float] = None,
        field_capacity: Union[Raster, np.ndarray, float] = None,
        wilting_point: Union[Raster, np.ndarray, float] = None,
        Topt_C: Union[Raster, np.ndarray, float] = None,
        fAPARmax: Union[Raster, np.ndarray, float] = None,
        canopy_height_meters: Union[Raster, np.ndarray, float] = None,
        delta_Pa: Union[Raster, np.ndarray, float] = None,
        gamma_Pa: Union[Raster, np.ndarray, float] = GAMMA_PA,
        epsilon=None,
        beta_Pa: float = BETA_PA,
        PT_alpha: float = PT_ALPHA,
        field_capacity_scale: float = FIELD_CAPACITY_SCALE,
        minimum_Topt: float = MINIMUM_TOPT,
        field_capacity_directory: str = SOIL_CAPACITY_DIRECTORY,
        wilting_point_directory: str = SOIL_CAPACITY_DIRECTORY,
        canopy_height_directory: str = GEDI_DOWNLOAD_DIRECTORY,
        floor_Topt: bool = FLOOR_TOPT,
        upscale_to_daily: bool = UPSCALE_TO_DAILY,
        regenerate_net_radiation: bool = False,
        resampling: str = RESAMPLING) -> Dict[str, Union[Raster, np.ndarray]]:
    """
    PTJPLSM: Compute partitioned evapotranspiration using the PT-JPL-SM model.

    Parameters:
        NDVI: Normalized Difference Vegetation Index (Raster or np.ndarray)
        Rn_Wm2: Net radiation (W/m^2) (Raster or np.ndarray)
        geometry: RasterGeometry object (optional)
        time_UTC: Datetime object for the observation (optional)
        hour_of_day: Hour of day (np.ndarray, optional)
        day_of_year: Day of year (np.ndarray, optional)
        GEOS5FP_connection: GEOS5FP meteorology connection (optional)
        ST_C: Surface temperature in Celsius (Raster or np.ndarray, optional)
        emissivity: Surface emissivity (Raster or np.ndarray, optional)
        albedo: Surface albedo (Raster or np.ndarray, optional)
        G: Soil heat flux (Raster or np.ndarray, optional)
        Ta_C: Air temperature in Celsius (Raster or np.ndarray, optional)
        RH: Relative humidity (0-1) (Raster or np.ndarray, optional)
        soil_moisture: Soil moisture (Raster or np.ndarray, optional)
        field_capacity: Soil field capacity (Raster or np.ndarray, optional)
        wilting_point: Soil wilting point (Raster or np.ndarray, optional)
        Topt: Optimal plant temperature (Raster or np.ndarray, optional)
        fAPARmax: Maximum fAPAR (Raster or np.ndarray, optional)
        canopy_height_meters: Canopy height (Raster or np.ndarray, optional)
        delta_Pa: Slope of SVP curve (Pa/degC, optional)
        gamma_Pa: Psychrometric constant (Pa, optional)
        epsilon: Ratio delta/(delta+gamma) (optional)
        beta_Pa: Model parameter (float, optional)
        PT_alpha: Priestley-Taylor alpha (float, optional)
        field_capacity_scale: Field capacity scaling factor (float, optional)
        minimum_Topt: Minimum allowed Topt (float, optional)
        field_capacity_directory: Directory for field capacity data (str, optional)
        wilting_point_directory: Directory for wilting point data (str, optional)
        canopy_height_directory: Directory for canopy height data (str, optional)
        floor_Topt: Whether to floor Topt to Ta_C (bool, optional)
        resampling: Resampling method (str, optional)

    Returns:
        Dictionary with keys: 'G', 'Rn_soil', 'LE_soil', 'Rn_canopy', 'PET', 'LE_canopy', 'LE_interception', 'LE'

    Example:
        The following example demonstrates how to use PTJPLSM with ECOSTRESS data:

        ```python
        from PTJPLSM import PTJPLSM
        # Assume you have already loaded the following variables from ECOSTRESS granules:
        # geometry, time_UTC, NDVI, Ta_C, RH, Rn, ST_C, albedo

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

    References:
        - Purdy et al. (2018), "PT-JPL-SM: A PT-JPL model variant incorporating soil moisture stress for improved global evapotranspiration partitioning."
    """
    results = {}

    t = TicToc()
    t.tic()
    logger.info("starting PT-JPL-SM model run")

    # If geometry is not provided, try to extract from NDVI raster
    if geometry is None and isinstance(NDVI, Raster):
        geometry = NDVI.geometry

    # Load Topt and fAPARmax if not provided
    if Topt_C is None and geometry is not None:
        logger.info("loading optimum temperature (Topt_C)")
        Topt_C = load_Topt(geometry)
    elif Topt_C is not None:
        logger.info("using given optimum temperature (Topt_C)")

    check_distribution(Topt_C, "Topt_C")

    if fAPARmax is None and geometry is not None:
        logger.info("loading maximum fAPARmax")
        fAPARmax = load_fAPARmax(geometry)

    check_distribution(fAPARmax, "fAPARmax")

    # Create GEOS5FP connection if not provided
    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    # Load air temperature if not provided
    if Ta_C is None and geometry is not None and time_UTC is not None:
        logger.info("retrieving air temperature (Ta_C) from GEOS-5 FP")
        Ta_C = GEOS5FP_connection.Ta_C(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )
    elif Ta_C is not None:
        logger.info("using given air temperature (Ta_C)")
    elif Ta_C is None:
        raise ValueError("air temperature (Ta_C) not given")

    check_distribution(Ta_C, "Ta_C")

    # Load relative humidity if not provided
    if RH is None and geometry is not None and time_UTC is not None:
        logger.info("retrieving relative humidity (RH) from GEOS-5 FP")
        RH = GEOS5FP_connection.RH(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )
    elif RH is not None:
        logger.info("using given relative humidity (RH)")

    if RH is None:
        raise ValueError("relative humidity (RH) not given")

    check_distribution(RH, "RH")

    # Load soil moisture if not provided
    if soil_moisture is None and geometry is not None and time_UTC is not None:
        logger.info("retrieving soil moisture (SM) from GEOS-5 FP")
        soil_moisture = GEOS5FP_connection.SM(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )
    elif soil_moisture is not None:
        logger.info("using given soil moisture (SM)")

    if soil_moisture is None:
        raise ValueError("soil moisture not given")

    check_distribution(soil_moisture, "soil_moisture")

    # Load field capacity if not provided
    if field_capacity is None and geometry is not None:
        logger.info("loading field capacity")
        field_capacity = load_field_capacity(
            geometry=geometry,
            directory=field_capacity_directory,
            resampling=resampling
        )
    elif field_capacity is not None:
        logger.info("using given field capacity")

    check_distribution(field_capacity, "field_capacity")

    # Load wilting point if not provided
    if wilting_point is None and geometry is not None:
        logger.info("loading wilting point")
        wilting_point = load_wilting_point(
            geometry=geometry, 
            directory=wilting_point_directory,
            resampling=resampling
        )
    elif wilting_point is not None:
        logger.info("using given wilting point")

    check_distribution(wilting_point, "wilting_point")

    # Load canopy height if not provided
    if canopy_height_meters is None and geometry is not None:
        logger.info("loading canopy height")
        canopy_height_meters = load_canopy_height(
            geometry=geometry, 
            source_directory=canopy_height_directory,
            resampling=resampling
        )
    elif canopy_height_meters is not None:
        logger.info("using given canopy height")

    check_distribution(canopy_height_meters, "canopy_height_meters")

    # If net radiation is not provided, compute from components
    if regenerate_net_radiation or (Rn_Wm2 is None and albedo is not None and ST_C is not None and emissivity is not None):
        if SWin_Wm2 is None and geometry is not None and time_UTC is not None:
            logger.info("retrieving shortwave radiation (SWin_Wm2) from GEOS-5 FP")
            SWin_Wm2 = GEOS5FP_connection.SWin(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )
        elif SWin_Wm2 is not None:
            logger.info("using given shortwave radiation (SWin_Wm2)")
        
        if upscale_to_daily:
            logger.info("running Verma net radiation with daily upscaling")
        else:
            logger.info("running instantaneous Verma net radiation")

        Rn_results = verma_net_radiation(
            SWin_Wm2=SWin_Wm2,
            albedo=albedo,
            ST_C=ST_C,
            emissivity=emissivity,
            Ta_C=Ta_C,
            RH=RH,
            upscale_to_daily=upscale_to_daily,
        )

        Rn_Wm2 = Rn_results["Rn_Wm2"]
        
        if "Rn_daily_Wm2" in Rn_results:
            Rn_daily_Wm2 = Rn_results["Rn_daily_Wm2"]

    elif Rn_Wm2 is not None:
        logger.info("using given net radiation (Rn_Wm2) for PT-JPL-SM processing")

    if Rn_Wm2 is None:
            missing_vars = []
            if albedo is None:
                missing_vars.append('albedo')
            if ST_C is None:
                missing_vars.append('ST_C')
            if emissivity is None:
                missing_vars.append('emissivity')
            if missing_vars:
                raise ValueError(f"net radiation (Rn_Wm2) not given, and missing required variables to calculate: {', '.join(missing_vars)}")
            else:
                raise ValueError("net radiation (Rn_Wm2) not given and cannot be calculated")

    check_distribution(Rn_Wm2, "Rn_Wm2")
    results["Rn_Wm2"] = Rn_Wm2

    # Compute soil heat flux if not provided
    if G_Wm2 is None and Rn_Wm2 is not None and ST_C is not None and NDVI is not None and albedo is not None:
        logger.info("calculating soil heat flux (G_Wm2)")
        G_Wm2 = calculate_SEBAL_soil_heat_flux(
            Rn=Rn_Wm2,
            ST_C=ST_C,
            NDVI=NDVI,
            albedo=albedo
        )
    elif G_Wm2 is not None:
        logger.info("using given soil heat flux (G_Wm2)")

    if G_Wm2 is None:
        raise ValueError("soil heat flux (G_Wm2) not given, no Rn_Wm2, ST_C, NDVI, and albedo to calculate")
    
    check_distribution(G_Wm2, "G")
    results["G_Wm2"] = G_Wm2

    # --- Meteorological Calculations ---
    # Calculate saturation vapor pressure (SVP) from air temperature
    SVP_Pa = SVP_Pa_from_Ta_C(Ta_C)
    # Constrain RH between 0 and 1
    RH = rt.clip(RH, 0, 1)
    # Calculate actual vapor pressure
    Ea_Pa = RH * SVP_Pa
    # Calculate vapor pressure deficit (VPD)
    VPD_Pa = rt.clip(SVP_Pa - Ea_Pa, 0, None)
    # Calculate relative surface wetness
    fwet = calculate_relative_surface_wetness(RH)

    check_distribution(fwet, "fwet")

    # --- Vegetation Calculations ---
    # Convert NDVI to SAVI
    SAVI = SAVI_from_NDVI(NDVI)
    # Calculate fAPAR from SAVI
    fAPAR = fAPAR_from_SAVI(SAVI)
    # Calculate fIPAR from NDVI
    fIPAR = fIPAR_from_NDVI(NDVI)
    # Replace zero fIPAR with NaN
    fIPAR = np.where(fIPAR == 0, np.nan, fIPAR)
    # Calculate green canopy fraction (fg)
    fg = calculate_green_canopy_fraction(fAPAR, fIPAR)
    # Calculate plant moisture constraint (fM)
    fM = calculate_plant_moisture_constraint(fAPAR, fAPARmax)
    # Calculate soil moisture constraint (fREW)
    fREW = calculate_fREW(soil_moisture, field_capacity, wilting_point, field_capacity_scale)

    check_distribution(fREW, "fREW")

    check_distribution(Topt_C, "Topt_C")

    # Floor Topt to Ta_C if requested, then clip to minimum_Topt
    if floor_Topt:
        Topt_C = rt.where(Ta_C > Topt_C, Ta_C, Topt_C)

    Topt_C = rt.clip(Topt_C, minimum_Topt, None)

    check_distribution(Topt_C, "Topt_C")

    # Calculate plant temperature constraint (fT)
    fT = calculate_plant_temperature_constraint(Ta_C, Topt_C)
    
    # Calculate LAI from NDVI
    LAI = carlson_leaf_area_index(NDVI)

    # --- Partitioning Calculations ---
    # Calculate epsilon if not provided
    if epsilon is None:
        # If delta in Pascals is not provided, calculate from air temperature in Celcius
        if delta_Pa is None:
            delta_Pa = delta_Pa_from_Ta_C(Ta_C)

        epsilon = delta_Pa / (delta_Pa + gamma_Pa)

    check_distribution(epsilon, "epsilon")

    # --- Soil Evaporation ---
    # Net radiation of the soil
    Rn_soil_Wm2 = calculate_soil_net_radiation(Rn_Wm2, LAI)
    check_distribution(Rn_soil_Wm2, "Rn_soil_Wm2")
    results["Rn_soil_Wm2"] = Rn_soil_Wm2

    # Soil evaporation (LEs)
    LE_soil_Wm2 = calculate_soil_latent_heat_flux(Rn_soil_Wm2, G_Wm2, epsilon, fwet, fREW, PT_alpha)
    check_distribution(LE_soil_Wm2, "LE_soil_Wm2")
    results["LE_soil_Wm2"] = LE_soil_Wm2

    # --- Canopy Transpiration ---
    # Net radiation of the canopy
    Rn_canopy_Wm2 = Rn_Wm2 - Rn_soil_Wm2
    check_distribution(Rn_canopy_Wm2, "Rn_canopy_Wm2")
    results["Rn_canopy_Wm2"] = Rn_canopy_Wm2
    # Potential evapotranspiration (PET)
    PET_Wm2 = PT_alpha * epsilon * (Rn_Wm2 - G_Wm2)
    check_distribution(PET_Wm2, "PET")
    results["PET_Wm2"] = PET_Wm2
    # Canopy moisture constraint (fTRM)
    fTRM = calculate_fTRM(PET_Wm2, RH, canopy_height_meters, soil_moisture, field_capacity, wilting_point, fM)
    check_distribution(fTRM, "fTRM")
    # Canopy transpiration (LEc)
    LE_canopy_Wm2 = calculate_canopy_latent_heat_flux(Rn_canopy_Wm2, epsilon, fwet, fg, fT, fTRM, PT_alpha)
    check_distribution(LE_canopy_Wm2, "LE_canopy_Wm2")
    results["LE_canopy_Wm2"] = LE_canopy_Wm2

    # --- Interception Evaporation ---
    # Interception evaporation (LEi)
    LE_interception_Wm2 = calculate_interception(Rn_canopy_Wm2, epsilon, fwet, PT_alpha)
    check_distribution(LE_interception_Wm2, "LE_interception_Wm2")
    results["LE_interception_Wm2"] = LE_interception_Wm2

    # --- Combined Evapotranspiration ---
    # Total instantaneous evapotranspiration (LE)
    LE_Wm2 = LE_soil_Wm2 + LE_canopy_Wm2 + LE_interception_Wm2
    # Constrain LE between 0 and PET
    LE_Wm2 = np.clip(LE_Wm2, 0, PET_Wm2)
    check_distribution(LE_Wm2, "LE_Wm2")
    results["LE_Wm2"] = LE_Wm2

    if upscale_to_daily and time_UTC is not None:
        logger.info("started daily ET upscaling")
        t_et = TicToc()
        t_et.tic()

        if Rn_daily_Wm2 is None:
            logger.info("running daily net radiation integration")
            Rn_daily_Wm2 = daily_Rn_integration_verma(
                Rn_Wm2=Rn_Wm2,
                time_UTC=time_UTC,
                geometry=geometry
            )

        check_distribution(Rn_daily_Wm2, "Rn_daily_Rn_Wm2")
        results["Rn_daily_Wm2"] = Rn_daily_Wm2

        EF = rt.where((LE_Wm2 == 0) | ((Rn_Wm2 - G_Wm2) == 0), 0, LE_Wm2 / (Rn_Wm2 - G_Wm2))
        check_distribution(EF, "EF")
        results["EF"] = EF

        # Calculate latent heat flux during daylight
        LE_daylight_Wm2 = EF * Rn_daily_Wm2
        check_distribution(LE_daylight_Wm2, "LE_daylight_Wm2")
        results["LE_daylight_Wm2"] = LE_daylight_Wm2

        # Calculate daily ET
        # ET = daily_ET_from_daily_LE(LE_daylight_Wm2, datetime_UTC=time_UTC, geometry=geometry)

        daylight_hours = calculate_daylight(day_of_year=day_of_year, time_UTC=time_UTC, geometry=geometry)

        # convert length of day in hours to seconds
        daylight_seconds = daylight_hours * 3600.0

        LAMBDA_JKG_WATER_20C = 2450000.0

        # factor seconds out of watts to get joules and divide by latent heat of vaporization to get kilograms
        ET_daily_kg = rt.clip(LE_daylight_Wm2 * daylight_seconds / LAMBDA_JKG_WATER_20C, 0.0, None)

        check_distribution(ET_daily_kg, "ET_daily_kg")
        results["ET_daily_kg"] = ET_daily_kg

        elapsed_et = t_et.tocvalue()
        logger.info(f"completed daily ET upscaling (elapsed: {elapsed_et:.2f} seconds)")

    elapsed = t.tocvalue()
    logger.info(f"PT-JPL-SM model run complete (elapsed: {elapsed:.2f} seconds)")

    return results
