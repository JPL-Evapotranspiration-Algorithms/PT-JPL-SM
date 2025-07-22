from typing import Union
import numpy as np
import rasters as rt
from rasters import Raster

def calculate_fTRM(
        PET: Union[Raster, np.ndarray], 
        RH: Union[Raster, np.ndarray], 
        canopy_height_meters: Union[Raster, np.ndarray], 
        soil_moisture: Union[Raster, np.ndarray], 
        field_capacity: Union[Raster, np.ndarray], 
        wilting_point: Union[Raster, np.ndarray], 
        fM: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """
    calculates PT-JPL-SM fTRM term as an update to PT-JPL fM for canopy latent heat flux

    PET: potential evapotranspiration in watts per square meter
    RH: relative humidity between 0 and 1
    canopy_height_meters: height of the plant canopy in meters
    soil_moisture: soil moisture in cubic meters per cubic meters
    field_capacity: soil field capacity in cubic meters per cubic meters
    wilting_point: soil wilting point in cubic meters per cubic meters
    fM: PT-JPL plant moisture constraint
    """
    a = 0.1
    p = (1 / (1 + PET)) - (a / (1 + canopy_height_meters))
    CHscalar = np.sqrt(canopy_height_meters)
    with np.errstate(divide='ignore', invalid='ignore'):
        WPCH = rt.clip(rt.where(CHscalar == 0, 0, wilting_point / CHscalar), 0, 1)
        CR = (1 - p) * (field_capacity - WPCH) + WPCH
        fTREW = rt.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)
    RHSM = RH ** (4 * (1 - soil_moisture) * (1 - RH))
    fTRM = (1 - RHSM) * fM + RHSM * rt.where(np.isnan(fTREW), 0, fTREW)

    return fTRM
