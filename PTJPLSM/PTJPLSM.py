"""
This module contains the core physics of the PT-JPL algorithm.
"""
import logging
import warnings
from datetime import datetime
from os.path import join, abspath, expanduser
from typing import Callable, Dict, List

import numpy as np

import rasters as rt

__author__ = "Gregory Halverson"

from BESS import DEFAULT_DOWNSCALE_AIR, DEFAULT_DOWNSCALE_HUMIDITY, DEFAULT_DOWNSCALE_MOISTURE ## separate pip package? 

from gedi_canopy_height import GEDICanopyHeight ## separate pip package?

from geos5fp import GEOS5FP
from modisci import MODISCI
from SRTM import SRTM

from PTJPL import PTJPL ## assuming this will now be pulling in the PTJPL pip package?
from soil_capacity_wilting import SoilGrids

from rasters import Raster, RasterGeometry

DEFAULT_WORKING_DIRECTORY = "."
DEFAULT_GEOS5FP_DOWNLOAD = "GEOS5FP_download_directory"
DEFAULT_GEOS5FP_PRODUCTS = "GEOS5FP_products"
DEFAULT_PTJPL_INTERMEDIATE = "PTJPL_intermediate"

DEFAULT_RESAMPLING = "cubic"
DEFAULT_PREVIEW_QUALITY = 20

DEFAULT_OUTPUT_VARIABLES = [
    "Rn",
    "LE",
    "ETc",
    "ETi",
    "ETs",
    "ET",
    "ESI",
    "WUE",
    "SM"
]

FLOOR_TOPT = True

# Priestley-Taylor coefficient alpha
PT_ALPHA = 1.26
BETA = 1.0

# maximum proportion of soil heat flux out of net radiation to the soil
G_MAX_PROPORTION = 0.35
W_MAX_PROPORTION = 0.5

# psychrometric constant gamma in pascals over kelvin
# same as value for ventilated (Asmann type) psychrometers, with an air movement of some 5 m/s
# http://www.fao.org/docrep/x0490e/x0490e07.htm
PSYCHROMETRIC_GAMMA = 0.0662  # Pa/K

KRN = 0.6
KPAR = 0.5


class GEOS5FPNotAvailableError(IOError):
    pass


class PTJPLSM(PTJPL):
    logger = logging.getLogger(__name__)

    def __init__(
            self,
            working_directory: str = None,
            static_directory: str = None,
            SRTM_connection: SRTM = None,
            SRTM_download: str = None,
            GEOS5FP_connection: GEOS5FP = None,
            GEOS5FP_download: str = None,
            GEOS5FP_products: str = None,
            GEDI_connection: GEDICanopyHeight = None,
            GEDI_download: str = None,
            ORNL_connection: MODISCI = None,
            CI_directory: str = None,
            soil_grids_connection: SoilGrids = None,
            soil_grids_download: str = None,
            intermediate_directory=None,
            preview_quality: int = DEFAULT_PREVIEW_QUALITY,
            ANN_model: Callable = None,
            ANN_model_filename: str = None,
            resampling: str = DEFAULT_RESAMPLING,
            downscale_air: bool = DEFAULT_DOWNSCALE_AIR,
            downscale_humidity: bool = DEFAULT_DOWNSCALE_HUMIDITY,
            downscale_moisture: bool = DEFAULT_DOWNSCALE_MOISTURE,
            floor_Topt: bool = FLOOR_TOPT,
            save_intermediate: bool = False,
            include_preview: bool = True,
            show_distribution: bool = True):
        if working_directory is None:
            working_directory = DEFAULT_WORKING_DIRECTORY

        working_directory = abspath(expanduser(working_directory))

        if static_directory is None:
            static_directory = working_directory

        static_directory = abspath(expanduser(static_directory))

        if intermediate_directory is None:
            intermediate_directory = join(working_directory, DEFAULT_PTJPL_INTERMEDIATE)

        if soil_grids_connection is None:
            soil_grids_connection = SoilGrids(
                working_directory=static_directory,
                source_directory=soil_grids_download
            )

        super(PTJPLSM, self).__init__(
            working_directory=working_directory,
            static_directory=static_directory,
            SRTM_connection=SRTM_connection,
            SRTM_download=SRTM_download,
            GEOS5FP_connection=GEOS5FP_connection,
            GEOS5FP_download=GEOS5FP_download,
            GEOS5FP_products=GEOS5FP_products,
            GEDI_connection=GEDI_connection,
            GEDI_download=GEDI_download,
            ORNL_connection=ORNL_connection,
            CI_directory=CI_directory,
            intermediate_directory=intermediate_directory,
            preview_quality=preview_quality,
            ANN_model=ANN_model,
            ANN_model_filename=ANN_model_filename,
            resampling=resampling,
            save_intermediate=save_intermediate,
            show_distribution=show_distribution,
            include_preview=include_preview,
            downscale_air=downscale_air,
            downscale_humidity=downscale_humidity,
            downscale_moisture=downscale_moisture,
            floor_Topt=floor_Topt
        )

        self.soil_grids = soil_grids_connection

    def FC(self, geometry: RasterGeometry, resampling: str = None):
        return self.soil_grids.FC(geometry=geometry, resampling=resampling)

    def WP(self, geometry: RasterGeometry, resampling: str = None):
        return self.soil_grids.WP(geometry=geometry, resampling=resampling)

    def fREW(self, SM: Raster, FC: Raster, WP: Raster, FC_scale: float = 0.7) -> Raster:
        # SMWP = SM - WP
        # FCWP = FC * FC_scale - WP
        # fREW = rt.clip(rt.where(FCWP == 0, 1, SMWP / FCWP), 0, 1)
        SMWP = rt.clip(SM - WP, 0, 1)
        FCWP = rt.clip(FC * FC_scale - WP, 0, 1)
        fREW = rt.clip(rt.where(FCWP == 0, 0, SMWP / FCWP), 0, 1)

        return fREW

    def PTJPL(
            self,
            geometry: RasterGeometry,
            target: str,
            time_UTC: datetime or str,
            ST_C: Raster,
            emissivity: Raster,
            NDVI: Raster,
            albedo: Raster,
            water: Raster = None,
            cloud_mask: Raster = None,
            SWin: Raster = None,
            Ta_C: Raster = None,
            RH: Raster = None,
            Ea_kPa: Raster = None,
            Topt: Raster = None,
            fAPARmax: Raster = None,
            Rn: Raster = None,
            Rn_daily: Raster = None,
            G: Raster = None,
            SM: Raster = None,
            wind_speed: Raster = None,
            output_variables: List[str] = DEFAULT_OUTPUT_VARIABLES) -> Dict[str, Raster]:

        STEFAN_BOLTZMAN_CONSTANT = 5.67036713e-8  # SI units watts per square meter per kelvin to the fourth
        warnings.filterwarnings('ignore')

        results = {}

        # calculate time
        hour_of_day = self.hour_of_day(time_UTC=time_UTC, geometry=geometry)
        day_of_year = self.day_of_year(time_UTC=time_UTC, geometry=geometry)
        SHA = self.SHA_deg_from_doy_lat(day_of_year, geometry.lat)
        sunrise_hour = self.sunrise_from_sha(SHA)
        daylight_hours = self.daylight_from_sha(SHA)
        date_UTC = time_UTC.date()

        self.diagnostic(ST_C, "ST_C", date_UTC, target)

        # calculate meteorology
        ST_K = ST_C + 273.15

        if Ta_C is None:
            Ta_C = self.Ta_C(time_UTC=time_UTC, geometry=geometry, ST_K=ST_K)

        self.diagnostic(Ta_C, "Ta_C", date_UTC, target)

        # calculate saturation vapor pressure in kPa from air temperature in celsius
        # floor saturation vapor pressure at 1
        SVP_kPa = rt.clip(self.SVP_kPa_from_Ta_C(Ta_C), 1, None)
        self.diagnostic(SVP_kPa, "SVP_kPa", date_UTC, target)

        if Ea_kPa is None and RH is None:
            Ea_kPa = self.Ea_kPa(time_UTC=time_UTC, geometry=geometry, ST_K=ST_K)
        elif Ea_kPa is None and RH is not None:
            RH = rt.clip(RH, 0, 1)
            Ea_kPa = RH * SVP_kPa

        self.diagnostic(Ea_kPa, "Ea_kPa", date_UTC, target)

        # calculate vapor pressure deficit from water vapor pressure
        VPD_kPa = SVP_kPa - Ea_kPa

        # lower bound of vapor pressure deficit is zero, negative values replaced with nodata
        VPD_kPa = rt.where(VPD_kPa < 0, np.nan, VPD_kPa)

        # calculate relative humidity from water vapor pressure and saturation vapor pressure
        # upper bound of relative humidity is one, results higher than one are capped at one
        if RH is None:
            RH = rt.clip(Ea_kPa / SVP_kPa, 0, 1)

        self.diagnostic(RH, "RH", date_UTC, target)

        if SWin is None:
            # SWin = self.SWin(time_UTC=time_UTC, geometry=geometry)
            Ra, Rg, UV, VIS, NIR, VISdiff, NIRdiff, VISdir, NIRdir = self.FLiES(
                geometry=geometry,
                target=target,
                time_UTC=time_UTC,
                albedo=albedo
            )

            SWin = Rg

        if "SWin" in output_variables:
            results["SWin"] = SWin

        self.diagnostic(SWin, "SWin", date_UTC, target)

        albedo = rt.clip(albedo, 0, 1)
        self.diagnostic(albedo, "albedo", date_UTC, target)

        # calculate outgoing shortwave from incoming shortwave and albedo
        SWout = rt.clip(SWin * albedo, 0, None)
        self.diagnostic(SWout, "SWout", date_UTC, target)

        # calculate instantaneous net radiation from components
        SWnet = rt.clip(SWin - SWout, 0, None)
        self.diagnostic(SWnet, "SWnet", date_UTC, target)

        if Rn is None:
            Rn = self.Rn(
                date_UTC=date_UTC,
                target=target,
                SWin=SWin,
                albedo=albedo,
                ST_C=ST_C,
                emissivity=emissivity,
                Ea_kPa=Ea_kPa,
                Ta_C=Ta_C,
                cloud_mask=cloud_mask
            )

        self.diagnostic(Rn, "Rn", date_UTC, target)

        if "Rn" in output_variables:
            results["Rn"] = Rn

        if Rn_daily is None:
            # integrate net radiation to daily value
            Rn_daily = self.Rn_daily(
                Rn,
                hour_of_day,
                sunrise_hour,
                daylight_hours
            )

            # constrain negative values of daily integrated net radiation
            Rn_daily = rt.clip(Rn_daily, 0, None)

        self.diagnostic(Rn_daily, "Rn_daily", date_UTC, target)

        if "Rn_daily" in output_variables:
            results["Rn_daily"] = Rn_daily

        # calculate relative surface wetness from relative humidity
        fwet = self.fwet_from_RH(RH)
        self.diagnostic(fwet, "fwet", date_UTC, target)

        if "fwet" in output_variables:
            results["fwet"] = fwet

        # calculate slope of saturation to vapor pressure curve Pa/K
        delta = self.delta_from_Ta(Ta_C)
        self.diagnostic(delta, "delta", date_UTC, target)

        # calculate vegetation values

        self.diagnostic(NDVI, "NDVI", date_UTC, target)

        if water is None:
            water = NDVI < 0

        # water-mask NDVI
        NDVI = rt.where(water, np.nan, NDVI)

        # calculate fAPAR from NDVI
        fAPAR = self.fAPAR_from_NDVI(NDVI)
        self.diagnostic(fAPAR, "fAPAR", date_UTC, target)

        # calculate fIPAR from NDVI
        fIPAR = self.fIPAR_from_NDVI(NDVI)
        self.diagnostic(fIPAR, "fIPAR", date_UTC, target)

        # calculate green canopy fraction (fg) from fAPAR and fIPAR, constrained between zero and one
        fg = rt.clip(fAPAR / fIPAR, 0, 1)
        self.diagnostic(fg, "fg", date_UTC, target)

        if "fg" in output_variables:
            results["fg"] = fg

        if fAPARmax is None:
            fAPARmax = self.load_fAPARmax(geometry=geometry)

        self.diagnostic(fAPARmax, "fAPARmax", date_UTC, target)

        # calculate plant moisture constraint (fM) from fraction of photosynthetically active radiation,
        # constrained between zero and one
        fM = rt.clip(fAPAR / fAPARmax, 0.0, 1.0)
        self.diagnostic(fM, "fM", date_UTC, target)

        if "fM" in output_variables:
            results["fM"] = fM

        if SM is None:
            SM = self.SM(
                time_UTC=time_UTC,
                geometry=geometry,
                ST_fine=ST_C.mask(~water),
                NDVI_fine=NDVI.mask(~water),
                water=water
            )

            # # calculate soil moisture constraint from mean relative humidity and vapor pressure deficit,
            # # constrained between zero and one
            # fSM = rt.clip(RH ** (VPD_kPa / BETA), 0.0, 1.0)

        if "SM" in output_variables:
            results["SM"] = SM

        self.diagnostic(SM, "SM", date_UTC, target)

        FC = self.FC(geometry=geometry)
        self.diagnostic(FC, "FC", date_UTC, target)
        WP = self.WP(geometry=geometry)
        self.diagnostic(WP, "WP", date_UTC, target)

        fREW = self.fREW(SM=SM, FC=FC, WP=WP)
        # fREW = (SM - WP) / (FC - WP)
        self.diagnostic(fREW, "fREW", date_UTC, target)

        if "fREW" in output_variables:
            results["fREW"] = fREW

        if Topt is None:
            Topt = self.load_Topt(geometry=geometry)

            if self.floor_Topt:
                Topt = rt.where(Ta_C > Topt, Ta_C, Topt)

        self.diagnostic(Topt, "Topt", date_UTC, target)

        if "Topt" in output_variables:
            results["Topt"] = Topt

        # calculate plant temperature constraint (fT) from optimal phenology
        fT = np.exp(-(((Ta_C - Topt) / Topt) ** 2))

        self.diagnostic(fT, "fT", date_UTC, target)

        if "fT" in output_variables:
            results["fT"] = fT

        # calculate leaf area index
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            LAI = -np.log(1 - fIPAR) * (1 / KPAR)

        self.diagnostic(LAI, "LAI", date_UTC, target)

        # calculate delta / (delta + gamma)
        epsilon = delta / (delta + PSYCHROMETRIC_GAMMA)
        self.diagnostic(epsilon, "epsilon", date_UTC, target)

        # soil evaporation

        # caluclate net radiation of the soil from leaf area index
        Rn_soil = Rn * np.exp(-KRN * LAI)
        self.diagnostic(Rn_soil, "Rn_soil", date_UTC, target)

        if G is None:
            # calculate instantaneous soil heat flux from net radiation and fractional vegetation cover
            G = Rn * (0.05 + (1 - rt.clip(fIPAR, 0, 1)) * 0.265)

            # constrain soil heat flux
            G = rt.where(np.isnan(G), np.nan, rt.clip(G, 0, G_MAX_PROPORTION * Rn))

        G = G.mask(~water)
        self.diagnostic(G, "G", date_UTC, target)

        # dew-point temperature in Celsius
        Td_C = Ta_C - ((100 - RH * 100) / 5.0)
        self.diagnostic(Td_C, "Td_C", date_UTC, target)

        if wind_speed is None:
            wind_speed = self.wind_speed(time_UTC=time_UTC, geometry=geometry)

        # water heat flux
        W = self.water_heat_flux(
            ST_C,
            Td_C,
            wind_speed,
            SWnet
        )

        W = rt.clip(W, 0, W_MAX_PROPORTION * Rn)
        W = W.mask(water)
        self.diagnostic(W, "W", date_UTC, target, blank_OK=True)

        # calculate soil evaporation (LEs) from relative surface wetness, soil moisture constraint,
        # priestley taylor coefficient, epsilon = delta / (delta + gamma), net radiation of the soil,
        # and soil heat flux
        LE_soil = rt.clip((fwet + fREW * (1 - fwet)) * PT_ALPHA * epsilon * (Rn_soil - G), 0, None)
        self.diagnostic(LE_soil, "LE_soil", date_UTC, target)

        # canopy transpiration

        # calculate net radiation of the canopy from net radiation of the soil
        Rn_canopy = Rn - Rn_soil
        self.diagnostic(Rn_canopy, "Rn_canopy", date_UTC, target)

        # calculate potential evapotranspiration (pET) from net radiation, and soil heat flux

        PET_water = PT_ALPHA * epsilon * (Rn - W)
        self.diagnostic(PET_water, "PET_water", date_UTC, target, blank_OK=True)
        PET_land = PT_ALPHA * epsilon * (Rn - G)
        self.diagnostic(PET_land, "PET_land", date_UTC, target)

        PET = rt.where(
            water,
            PET_water,
            PET_land
        )

        if "PET" in output_variables:
            results["PET"] = PET

        self.diagnostic(PET, "PET", date_UTC, target)

        CH = self.canopy_height_meters(geometry=geometry)
        self.diagnostic(CH, "CH", date_UTC, target)
        a = 0.1
        p = (1 / (1 + PET)) - (a / (1 + CH))
        self.diagnostic(p, "p", date_UTC, target)
        CHscalar = np.sqrt(CH)
        self.diagnostic(CHscalar, "CHscalar", date_UTC, target)
        WPCH = rt.clip(rt.where(CHscalar == 0, 0, WP / CHscalar), 0, 1)
        self.diagnostic(WPCH, "WPCH", date_UTC, target)
        CR = (1 - p) * (FC - WPCH) + WPCH
        self.diagnostic(CR, "CR", date_UTC, target)

        fTREW = rt.clip(1 - ((CR - SM) / (CR - WPCH)) ** CHscalar, 0, 1)
        self.diagnostic(fTREW, "fTREW", date_UTC, target)

        if "fTREW" in output_variables:
            results["fTREW"] = fTREW

        RHSM = RH ** (4 * (1 - SM) * (1 - RH))

        fTRM = (1 - RHSM) * fM + RHSM * rt.where(np.isnan(fTREW), 0, fTREW)
        self.diagnostic(fTRM, "fTRM", date_UTC, target)

        if "fTRM" in output_variables:
            results["fTRM"] = fTRM

        # calculate canopy transpiration (LEc) from priestley taylor, relative surface wetness,
        # green canopy fraction, plant temperature constraint, plant moisture constraint,
        # epsilon = delta / (delta + gamma), and net radiation of the canopy
        LE_canopy = rt.clip(PT_ALPHA * (1 - fwet) * fg * fT * fTRM * epsilon * Rn_canopy, 0, None)
        self.diagnostic(LE_canopy, "LE_canopy", date_UTC, target)

        if "LE_canopy" in output_variables:
            results["LE_canopy"] = LE_canopy

        # interception evaporation

        # calculate interception evaporation (LEi) from relative surface wetness and net radiation of the canopy
        LE_interception = rt.clip(fwet * PT_ALPHA * epsilon * Rn_canopy, 0, None)
        self.diagnostic(LE_interception, "LE_interception", date_UTC, target)

        # combined evapotranspiration

        # combine soil evaporation (LEs), canopy transpiration (LEc), and interception evaporation (LEi)
        # into instantaneous evapotranspiration (LE)
        LE = LE_soil + LE_canopy + LE_interception
        soil_proportion = rt.where((LE_soil == 0) | (LE == 0), 0, LE_soil / LE)
        self.diagnostic(soil_proportion, "soil_proportion", date_UTC, target)

        if "soil_proportion" in output_variables:
            results["soil_proportion"] = soil_proportion

        canopy_proportion = rt.where((LE_canopy == 0) | (LE == 0), 0, LE_canopy / LE)
        self.diagnostic(canopy_proportion, "canopy_proportion", date_UTC, target)

        if "canopy_proportion" in output_variables:
            results["canopy_proportion"] = canopy_proportion

        interception_proportion = rt.where((LE_interception == 0) | (LE == 0), 0, LE_interception / LE)
        self.diagnostic(interception_proportion, "interception_proportion", date_UTC, target)

        if "interception_proportion" in output_variables:
            results["interception_proportion"] = interception_proportion

        LE = rt.where(water, PET, LE)
        LE = rt.where(np.isnan(LE), PET, LE)
        LE = rt.clip(LE, 0, PET)

        if "LE" in output_variables:
            results["LE"] = LE

        self.diagnostic(LE, "LE", date_UTC, target)

        # daily evapotranspiration

        # calculate evaporative fraction (EF) from evapotranspiration, net radiation, and soil heat flux

        EF = rt.clip(rt.where((LE == 0) | (Rn == 0), 0, LE / Rn), 0, 1)
        self.diagnostic(EF, "EF", date_UTC, target)

        # calculate daily latent heat flux from evaporative fraction and daily net radiation with minimum of zero
        LE_daily = rt.clip(EF * Rn_daily, 0, Rn_daily)
        self.diagnostic(LE_daily, "LE_daily", date_UTC, target)

        # calculate instantaneous evaporative stress index as the ratio of latent heat flux to potential evapotranspiration
        ESI = rt.clip(rt.where((LE == 0) | (PET == 0), 0, LE / PET), 0, 1)
        ESI = ESI.mask(~water)

        if "ESI" in output_variables:
            results["ESI"] = ESI

        self.diagnostic(ESI, "ESI", date_UTC, target)

        # calculate daytime daily evapotranspiration in kilograms equivalent to millimeters
        ET = self.calculate_vapor(LE_daily, daylight_hours)

        if "ET" in output_variables:
            results["ET"] = ET

        self.diagnostic(ET, "ET", date_UTC, target)

        warnings.resetwarnings()

        return results
