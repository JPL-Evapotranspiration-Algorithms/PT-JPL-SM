from typing import Union
import numpy as np
import rasters as rt
from rasters import Raster

from ..constants import FIELD_CAPACITY_SCALE

def calculate_fREW(
        soil_moisture: Union[Raster, np.ndarray], 
        field_capacity: Union[Raster, np.ndarray], 
        wilting_point: Union[Raster, np.ndarray], 
        field_capacity_scale: float = FIELD_CAPACITY_SCALE) -> Union[Raster, np.ndarray]:
    """
    Calculate the fraction of relative extractable water (fREW) for the PT-JPL-SM evapotranspiration model.

    This function is based on the equation 6 from Purdy et al., 2018. It computes the fREW, which represents the 
    proportion of water in the soil that plants can extract relative to the total available water between the 
    wilting point and field capacity.

    This constraint is the PT-JPL-SM replacement for fSM in the original PT-JPL.

    Parameters:
    soil_moisture (Union[Raster, np.ndarray]): soil moisture in cubic meters per cubic meters
    field_capacity (Union[Raster, np.ndarray]): soil field capacity in cubic meters per cubic meters
    wilting_point (Union[Raster, np.ndarray]): soil-plant wilting-point in cubic meters per cubic meters
    field_capacity_scale (float, optional): A scaling factor applied to the field capacity. Defaults to 0.7.

    Returns:
    Union[Raster, np.ndarray]: The calculated fREW, clipped between 0 and 1.
    """
    # calculate difference of soil moisture and wilting point
    SMWP = rt.clip(soil_moisture - wilting_point, 0, 1)
    # calculate difference of field capacity and wilting point
    FCWP = rt.clip(field_capacity * field_capacity_scale - wilting_point, 0, 1)
    # calculate fraction of relative extractable water (fREW)
    with np.errstate(divide='ignore', invalid='ignore'):
        fREW = rt.clip(rt.where(FCWP == 0, 0, SMWP / FCWP), 0, 1)

    return fREW
