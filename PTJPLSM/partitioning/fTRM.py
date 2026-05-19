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
    Calculates the PT-JPL-SM Transpiration Reduction Modifier (fTRM) term.
    This serves as an update to the standard PT-JPL plant moisture constraint (fM)
    for canopy latent heat flux by integrating explicit soil moisture dynamics.

    Parameters:
    ----------
    PET : Union[Raster, np.ndarray]
        Potential evapotranspiration in watts per square meter (W/m²).
    RH : Union[Raster, np.ndarray]
        Relative humidity scaled between 0 and 1.
    canopy_height_meters : Union[Raster, np.ndarray]
        Height of the plant canopy in meters.
    soil_moisture : Union[Raster, np.ndarray]
        Volumetric soil moisture in cubic meters per cubic meter (m³/m³).
    field_capacity : Union[Raster, np.ndarray]
        Soil field capacity in cubic meters per cubic meter (m³/m³).
    wilting_point : Union[Raster, np.ndarray]
        Soil wilting point in cubic meters per cubic meter (m³/m³).
    fM : Union[Raster, np.ndarray]
        Original PT-JPL plant moisture constraint based on atmospheric demand.

    Returns:
    -------
    Union[Raster, np.ndarray]
        The updated canopy moisture constraint (fTRM), scaled between 0 and 1.
    """
    # -------------------------------------------------------------------------
    # 1. CANOPY HEIGHT SCALING & ATMOSPHERIC SENSITIVITY
    # -------------------------------------------------------------------------
    # 'p' is an empirical parameter adjusting the soil moisture stress threshold 
    # based on atmospheric demand (PET) and physical vegetation stature. Higher 
    # atmospheric demand shifts the soil moisture threshold for plant stress.
    a = 0.1
    p = (1 / (1 + PET)) - (a / (1 + canopy_height_meters))
    
    # 'CHscalar' acts as a proxy for aerodynamic resistance and hydraulic capacitance.
    # The square root of canopy height is a standard scaling factor in these models.
    CHscalar = np.sqrt(canopy_height_meters)
    
    # Suppress runtime warnings caused by division by zero or NaN values over 
    # invalid pixels (e.g., open water body pixels or non-vegetated regions).
    with np.errstate(divide='ignore', invalid='ignore'):
        
        # 'WPCH' scales the wilting point downward for taller canopies, accounting for 
        # deeper root networks and greater water extraction capabilities under tension.
        # rt.where safely sets bare soil scenarios (CHscalar == 0) to 0.
        WPCH = rt.clip(rt.where(CHscalar == 0, 0, wilting_point / CHscalar), 0, 1)
        
        # ---------------------------------------------------------------------
        # 2. CRITICAL MOISTURE POINT (CR)
        # ---------------------------------------------------------------------
        # Represents the specific soil moisture threshold below which the vegetation 
        # begins to actively experience transpiration reduction/moisture stress.
        CR = (1 - p) * (field_capacity - WPCH) + WPCH
        
        # ---------------------------------------------------------------------
        # 3. TRANSPIRATION REDUCTION EVAPORATIVE WATER STRESS (fTREW)
        # ---------------------------------------------------------------------
        # Calculates raw soil moisture stress. As soil moisture drops below CR, 
        # fTREW scales down toward 0. The rate of this drop-off is non-linearly 
        # driven by CHscalar as an exponent.
        fTREW = rt.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)
        
    # -------------------------------------------------------------------------
    # 4. RELATIVE HUMIDITY BUFFERING & FINAL INTEGRATION
    # -------------------------------------------------------------------------
    # 'RHSM' serves as a dynamic weighting factor. It determines the balance 
    # between atmospheric control (fM) and direct soil moisture control (fTREW).
    RHSM = RH ** (4 * (1 - soil_moisture) * (1 - RH))
    
    # 'fTRM' blends the original atmospheric constraint (fM) with the new soil moisture
    # constraint (fTREW) using RHSM as the slider. NaNs in fTREW (e.g., where CR == WPCH)
    # are caught and defaulted safely to a highly-stressed zero condition.
    fTRM = (1 - RHSM) * fM + RHSM * rt.where(np.isnan(fTREW), 0, fTREW)

    return fTRM