from typing import Union
import numpy as np
from rasters import Raster

from .fTREW import CANOPY_BUFFER_SENSITIVITY, calculate_fTREW

def calculate_fTRM(
    PET_Wm2: Union[Raster, np.ndarray], 
        RH: Union[Raster, np.ndarray], 
        canopy_height_meters: Union[Raster, np.ndarray], 
        soil_moisture: Union[Raster, np.ndarray], 
        field_capacity: Union[Raster, np.ndarray], 
        wilting_point: Union[Raster, np.ndarray], 
        fM: Union[Raster, np.ndarray],
        canopy_buffer_sensitivity: float = CANOPY_BUFFER_SENSITIVITY) -> Union[Raster, np.ndarray]:
    r"""
    Calculates the PT-JPL-SM Transpiration Reduction Modifier (fTRM) term.
    This serves as an update to the standard PT-JPL plant moisture constraint (fM)
    for canopy latent heat flux by integrating explicit soil moisture dynamics.

    Parameters:
    ----------
    PET_Wm2 : Union[Raster, np.ndarray]
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
    canopy_buffer_sensitivity : float, optional
        Empirical parameter adjusting the soil moisture stress threshold based on canopy height and atmospheric demand.

    Returns:
    -------
    Union[Raster, np.ndarray]
        The updated canopy moisture constraint (fTRM), scaled between 0 and 1.

    Mathematical Formulations & Sequential Logic:
    ----------------------------------------------
    The algorithm first computes the transpiration-side soil moisture stress scalar 
    ($f_{TREW}$, represented by `fTREW`) by delegating to `calculate_fTREW`, which applies 
    the canopy-height-adjusted soil moisture stress formulation documented in `fTREW.py`.

    To prevent a disconnect between topsoil drought and atmospheric humidity constraints, 
    a dynamic relative humidity soil moisture weighting factor ($RHSM$, represented by `RHSM`) 
    is introduced to balance atmospheric and direct soil moisture controls based on relative humidity 
    ($RH$) and volumetric soil moisture ($\theta_{obs}$):
    $$RHSM = RH^{4(1 - \theta_{obs})(1 - RH)}$$

    Finally, the original PT-JPL atmospheric plant moisture constraint ($f_M$, represented by `fM`) 
    and the new soil moisture stress constraint ($f_{TREW}$) are dynamically blended using $RHSM$ 
    as the slider to arrive at the comprehensive Transpiration Reduction Modifier ($f_{TRM}$, 
    represented by `fTRM`):
    $$f_{TRM} = (1 - RHSM) \times f_M + RHSM \times f_{TREW}$$

    Parameter Constraints & Eco-hydrological Bounds:
    -------------------------------------------------
    The canopy_buffer_sensitivity parameter ($a$) must be bounded within [0.0, 1.0]:
        * Setting to 0.0: Nullifies physical vegetation buffering in the delegated 
            $f_{TREW}$ calculation. The stress onset threshold becomes driven solely by 
            atmospheric demand (PET_Wm2), meaning a 30-meter forest and a 10-centimeter 
            grassland respond identically to topsoil drying.
    * Setting too high (> 1.0): Overpowers the atmospheric demand term and can drive 
            the delegated stress onset weight ($p$) negative. This unphysically pushes the 
            Critical Moisture Point ($\theta_{CR}$) above Field Capacity ($\theta_{FC}$), 
            causing the model to falsely simulate severe transpiration stress in fully 
            saturated soils.

    References:
    -----------
    1. Purdy, A. J., Fisher, J. B., Goulden, M. L., Colliander, A., Halverson, G. H., 
       Tu, K., & Famiglietti, J. S. (2018). SMAP soil moisture improves global 
       evapotranspiration. Remote Sensing of Environment, 219, 1-14. 
       https://doi.org/10.1016/j.rse.2018.09.023

    2. Fisher, J. B., Tu, K., & Baldocchi, D. D. (2008). Global estimates of the 
       land-atmosphere water flux based on monthly AVHRR and ISLSCP-II data, 
       validated at 16 FLUXNET sites. Remote Sensing of Environment, 112(3), 901-919.
       https://doi.org/10.1016/j.rse.2007.06.025

    3. Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat 
       flux and evaporation using large-scale parameters. Monthly Weather Review, 
       100(2), 81-92. https://doi.org/10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2
    """
    # Parameter Validity Check
    if not (0.0 <= canopy_buffer_sensitivity <= 1.0):
        raise ValueError(
            f"Invalid canopy_buffer_sensitivity ({canopy_buffer_sensitivity}). "
            f"Parameter must be bounded between 0.0 and 1.0 to preserve "
            f"eco-hydrological physical constraints."
        )

    # Compute the transpiration-side soil moisture stress scalar.
    fTREW = calculate_fTREW(
        PET_Wm2=PET_Wm2,
        canopy_height_meters=canopy_height_meters,
        soil_moisture=soil_moisture,
        field_capacity=field_capacity,
        wilting_point=wilting_point,
        canopy_buffer_sensitivity=canopy_buffer_sensitivity,
    )
        
    # Relative Humidity Buffering & Final Integration
    # 'RHSM' serves as a dynamic weighting factor based on RH and volumetric water content (\theta_obs).
    # It determines the balance between atmospheric control (fM) and direct soil moisture control (fTREW).
    RHSM = RH ** (4 * (1 - soil_moisture) * (1 - RH))
    
    # 'fTRM' (f_TRM) blends the original atmospheric constraint (fM) with the new soil moisture
    # constraint (fTREW) using RHSM as the slider.
    fTRM = (1 - RHSM) * fM + RHSM * fTREW

    return fTRM