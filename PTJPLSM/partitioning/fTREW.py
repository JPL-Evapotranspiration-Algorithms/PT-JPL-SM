from typing import Union

import numpy as np
import rasters as rt
from rasters import Raster

CANOPY_BUFFER_SENSITIVITY = 0.1


def calculate_fTREW(
        PET: Union[Raster, np.ndarray],
        canopy_height_meters: Union[Raster, np.ndarray],
        soil_moisture: Union[Raster, np.ndarray],
        field_capacity: Union[Raster, np.ndarray],
        wilting_point: Union[Raster, np.ndarray],
        canopy_buffer_sensitivity: float = CANOPY_BUFFER_SENSITIVITY) -> Union[Raster, np.ndarray]:
    r"""
    Calculate the transpiration-side relative extractable water constraint (fTREW).

    Parameters
    ----------
    PET : Union[Raster, np.ndarray]
        Potential evapotranspiration in watts per square meter (W/m²).
    canopy_height_meters : Union[Raster, np.ndarray]
        Height of the plant canopy in meters.
    soil_moisture : Union[Raster, np.ndarray]
        Volumetric soil moisture in cubic meters per cubic meter (m³/m³).
    field_capacity : Union[Raster, np.ndarray]
        Soil field capacity in cubic meters per cubic meter (m³/m³).
    wilting_point : Union[Raster, np.ndarray]
        Soil wilting point in cubic meters per cubic meter (m³/m³).
    canopy_buffer_sensitivity : float, optional
        Empirical parameter adjusting stress onset based on canopy height and atmospheric demand.

    Returns
    -------
    Union[Raster, np.ndarray]
        fTREW scalar clipped to [0, 1]. Any undefined values are corrected to 0.
    """
    if not (0.0 <= canopy_buffer_sensitivity <= 1.0):
        raise ValueError(
            f"Invalid canopy_buffer_sensitivity ({canopy_buffer_sensitivity}). "
            f"Parameter must be bounded between 0.0 and 1.0 to preserve "
            f"eco-hydrological physical constraints."
        )

    stress_onset_weight = (1 / (1 + PET)) - (canopy_buffer_sensitivity / (1 + canopy_height_meters))
    CHscalar = np.sqrt(canopy_height_meters)

    with np.errstate(divide='ignore', invalid='ignore'):
        WPCH = rt.clip(rt.where(CHscalar == 0, 0, wilting_point / CHscalar), 0, 1)
        CR = (1 - stress_onset_weight) * (field_capacity - WPCH) + WPCH
        fTREW = rt.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)

    # Final correction: replace undefined fTREW values with a conservative stressed condition.
    return rt.where(np.isnan(fTREW), 0, fTREW)
