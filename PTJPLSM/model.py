from typing import Union, Dict
import warnings
import numpy as np
import pandas as pd
from datetime import datetime
import rasters as rt
from rasters import Raster, RasterGeometry
from GEOS5FP import GEOS5FP

from soil_capacity_wilting import DEFAULT_DOWNLOAD_DIRECTORY as SOIL_CAPACITY_DIRECTORY
from soil_capacity_wilting import load_field_capacity, load_wilting_point
from gedi_canopy_height import GEDI_DOWNLOAD_DIRECTORY
from gedi_canopy_height import load_canopy_height

from PTJPL import GAMMA_PA
from PTJPL import BETA_PA
from PTJPL import PT_ALPHA
from PTJPL import MINIMUM_TOPT
from PTJPL import FLOOR_TOPT

from PTJPL import SVP_Pa_from_Ta_C
from PTJPL import delta_Pa_from_Ta_C
from PTJPL import LAI_from_NDVI
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
from PTJPL import process_verma_net_radiation
from PTJPL import calculate_SEBAL_soil_heat_flux

from .constants import *
from .partitioning import (
    calculate_fREW, calculate_fTRM,
    calculate_soil_latent_heat_flux, calculate_canopy_latent_heat_flux
)

def PTJPLSM(
        NDVI: Union[Raster, np.ndarray],
        Rn_Wm2: Union[Raster, np.ndarray],
        geometry: RasterGeometry = None,
        time_UTC: datetime = None,
        hour_of_day: np.ndarray = None,
        day_of_year: np.ndarray = None,
        GEOS5FP_connection: GEOS5FP = None,
        ST_C: Union[Raster, np.ndarray] = None,
        emissivity: Union[Raster, np.ndarray] = None,
        albedo: Union[Raster, np.ndarray] = None,
        G: Union[Raster, np.ndarray] = None,
        Ta_C: Union[Raster, np.ndarray] = None,
        RH: Union[Raster, np.ndarray] = None,
        soil_moisture: Union[Raster, np.ndarray] = None,
        field_capacity: Union[Raster, np.ndarray] = None,
        wilting_point: Union[Raster, np.ndarray] = None,
        Topt: Union[Raster, np.ndarray] = None,
        fAPARmax: Union[Raster, np.ndarray] = None,
        canopy_height_meters: Union[Raster, np.ndarray] = None,
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
        resampling: str = RESAMPLING) -> Dict[str, Union[Raster, np.ndarray]]:
    results = {}

    if geometry is None and isinstance(NDVI, Raster):
        geometry = NDVI.geometry

    if Topt is None and geometry is not None:
        Topt = load_Topt(geometry)

    if fAPARmax is None and geometry is not None:
        fAPARmax = load_fAPARmax(geometry)

    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    if Ta_C is None and geometry is not None and time_UTC is not None:
        Ta_C = GEOS5FP_connection.Ta_C(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    if Ta_C is None:
        raise ValueError("air temperature (Ta_C) not given")
    
    if RH is None and geometry is not None and time_UTC is not None:
        RH = GEOS5FP_connection.RH(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    if RH is None:
        raise ValueError("relative humidity (RH) not given")
    
    if soil_moisture is None and geometry is not None and time_UTC is not None:
        soil_moisture = GEOS5FP_connection.SM(
            time_UTC=time_UTC,
            geometry=geometry,
            resampling=resampling
        )

    if soil_moisture is None:
        raise ValueError("soil moisture not given")

    if field_capacity is None and geometry is not None:
        field_capacity = load_field_capacity(
            geometry=geometry,
            directory=filed_capacity_directory,
            resampling=resampling
        )
    
    if wilting_point is None and geometry is not None:
        wilting_point = load_wilting_point(
            geometry=geometry, 
            directory=wilting_point_directory,
            resampling=resampling
        )

    if canopy_height_meters is None and geometry is not None:
        canopy_height_meters = load_canopy_height(
            geometry=geometry, 
            directory=canopy_height_directory,
            resampling=resampling
        )

    if Rn_Wm2 is None and albedo is not None and ST_C is not None and emissivity is not None:
        if SWin is None and geometry is not None and time_UTC is not None:
            SWin = GEOS5FP_connection.SWin(
                time_UTC=time_UTC,
                geometry=geometry,
                resampling=resampling
            )

        Rn_results = process_verma_net_radiation(
            SWin=SWin,
            albedo=albedo,
            ST_C=ST_C,
            emissivity=emissivity,
            Ta_C=Ta_C,
            RH=RH
        )

        Rn_Wm2 = Rn_results["Rn"]

    if Rn_Wm2 is None:
        raise ValueError("net radiation (Rn) not given")

    if G is None and Rn_Wm2 is not None and ST_C is not None and NDVI is not None and albedo is not None:
        G = calculate_SEBAL_soil_heat_flux(
            Rn=Rn_Wm2,
            ST_C=ST_C,
            NDVI=NDVI,
            albedo=albedo
        )

    if G is None:
        raise ValueError("soil heat flux (G) not given")
    
    results["G"] = G

    # calculate meteorology

    # calculate saturation vapor pressure in kPa from air temperature in celsius
    # floor saturation vapor pressure at 1
    SVP_Pa = SVP_Pa_from_Ta_C(Ta_C)

    # constrain relative humidity between 0 and 1
    RH = rt.clip(RH, 0, 1)

    # calculate water vapor pressure in Pascals from relative humidity and saturation vapor pressure
    Ea_Pa = RH * SVP_Pa

    # calculate vapor pressure deficit from water vapor pressure
    VPD_Pa = rt.clip(SVP_Pa - Ea_Pa, 0, None)

    # calculate relative surface wetness from relative humidity
    fwet = calculate_relative_surface_wetness(RH)

    # calculate vegetation values

    # convert normalized difference vegetation index to soil-adjusted vegetation index
    SAVI = SAVI_from_NDVI(NDVI)

    # calculate fraction of absorbed photosynthetically active radiation from soil-adjusted vegetation index
    fAPAR = fAPAR_from_SAVI(SAVI)

    # calculate fIPAR from NDVI
    fIPAR = fIPAR_from_NDVI(NDVI)

    # replace zero fIPAR with NaN
    fIPAR = np.where(fIPAR == 0, np.nan, fIPAR)

    # calculate green canopy fraction (fg) from fAPAR and fIPAR, constrained between zero and one
    fg = calculate_green_canopy_fraction(fAPAR, fIPAR)

    # calculate plant moisture constraint (fM) from fraction of photosynthetically active radiation, constrained between zero and one
    fM = calculate_plant_moisture_constraint(fAPAR, fAPARmax)

    fREW = calculate_fREW(soil_moisture, field_capacity, wilting_point, field_capacity_scale)

    if floor_Topt:
        Topt = rt.where(Ta_C > Topt, Ta_C, Topt)

    Topt = rt.clip(Topt, minimum_Topt, None)

    # calculate plant temperature constraint (fT) from optimal phenology
    fT = calculate_plant_temperature_constraint(Ta_C, Topt)

    LAI = LAI_from_NDVI(NDVI)

    # calculate delta / (delta + gamma) term if it's not given
    if epsilon is None:
        # calculate delta if it's not given
        if delta_Pa is None:
            # calculate slope of saturation to vapor pressure curve in kiloPascal per degree Celsius
            delta_Pa = delta_Pa_from_Ta_C(Ta_C)

        # calculate delta / (delta + gamma)
        epsilon = delta_Pa / (delta_Pa + gamma_Pa)

    # soil evaporation

    # calculate net radiation of the soil from leaf area index
    Rn_soil = calculate_soil_net_radiation(Rn_Wm2, LAI)
    results["Rn_soil"] = Rn_soil

    # calculate soil evaporation (LEs) from relative surface wetness, soil moisture constraint,
    # priestley taylor coefficient, epsilon = delta / (delta + gamma), net radiation of the soil,
    # and soil heat flux
    LE_soil = calculate_soil_latent_heat_flux(Rn_soil, G, epsilon, fwet, fREW, PT_alpha)
    results["LE_soil"] = LE_soil

    # canopy transpiration

    # calculate net radiation of the canopy from net radiation of the soil
    Rn_canopy = Rn_Wm2 - Rn_soil
    results["Rn_canopy"] = Rn_canopy
    
    # calculate potential evapotranspiration (pET) from net radiation, and soil heat flux
    PET = PT_alpha * epsilon * (Rn_Wm2 - G)
    results["PET"] = PET

    fTRM = calculate_fTRM(PET, RH, canopy_height_meters, soil_moisture, field_capacity, wilting_point, fM)
    
    # calculate canopy transpiration (LEc) from priestley taylor, relative surface wetness,
    # green canopy fraction, plant temperature constraint, plant moisture constraint,
    # epsilon = delta / (delta + gamma), and net radiation of the canopy
    LE_canopy = calculate_canopy_latent_heat_flux(Rn_canopy, epsilon, fwet, fg, fT, fTRM, PT_alpha)
    results["LE_canopy"] = LE_canopy

    # interception evaporation

    # calculate interception evaporation (LEi) from relative surface wetness and net radiation of the canopy
    LE_interception = calculate_interception(Rn_canopy, epsilon, fwet, PT_alpha)
    results["LE_interception"] = LE_interception

    # combined evapotranspiration

    # combine soil evaporation (LEs), canopy transpiration (LEc), and interception evaporation (LEi)
    # into instantaneous evapotranspiration (LE)
    LE = LE_soil + LE_canopy + LE_interception
    
    # constrain instantaneous evapotranspiration between zero and potential evapotranspiration
    LE = np.clip(LE, 0, PET)
    results["LE"] = LE

    return results
