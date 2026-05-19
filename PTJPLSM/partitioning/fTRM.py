from typing import Union
import numpy as np
import rasters as rt
from rasters import Raster

CANOPY_BUFFER_SENSITIVITY = 0.1

def calculate_fTRM(
        PET: Union[Raster, np.ndarray], 
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
    canopy_buffer_sensitivity : float, optional
        Empirical parameter adjusting the soil moisture stress threshold based on canopy height and atmospheric demand.

    Returns:
    -------
    Union[Raster, np.ndarray]
        The updated canopy moisture constraint (fTRM), scaled between 0 and 1.

    Mathematical Formulations:
    --------------------------
    1. Stress Onset Weight (p):
       $$stress\_onset\_weight = \frac{1}{1 + PET} - \frac{canopy\_buffer\_sensitivity}{1 + canopy\_height\_meters}$$

    2. Canopy Height Scalar:
       $$CHscalar = \sqrt{canopy\_height\_meters}$$

    3. Canopy-Scaled Wilting Point (WPCH):
       $$WPCH = \begin{cases} 0, & \text{if } CHscalar = 0 \\ \text{clip}\left(\frac{wilting\_point}{CHscalar}, 0, 1\right), & \text{otherwise} \end{cases}$$

    4. Critical Moisture Point (CR):
       $$CR = (1 - stress\_onset\_weight) \times (field\_capacity - WPCH) + WPCH$$

    5. Transpiration Reduction Evaporative Water Stress (fTREW):
       $$fTREW = \text{clip}\left(1 - \left(\frac{CR - soil\_moisture}{CR - WPCH}\right)^{CHscalar}, 0, 1\right)$$

    6. Relative Humidity Soil Moisture Weighting Factor (RHSM):
       $$RHSM = RH^{4 \times (1 - soil\_moisture) \times (1 - RH)}$$

    7. Final Transpiration Reduction Modifier (fTRM):
       $$fTRM = (1 - RHSM) \times fM + RHSM \times fTREW$$

    Parameter Constraints & Eco-hydrological Bounds:
    -------------------------------------------------
    The canopy_buffer_sensitivity parameter must be bounded within [0.0, 1.0]:
    * Setting to 0.0: Nullifies physical vegetation buffering. The stress onset 
      threshold becomes driven solely by atmospheric demand (PET), meaning a 
      30-meter forest and a 10-centimeter grassland respond identically to topsoil drying.
    * Setting too high (> 1.0): Overpowers the atmospheric demand term and can drive 
      stress_onset_weight negative. This unphysically pushes the Critical Moisture 
      Point (CR) above Field Capacity, causing the model to falsely simulate severe 
      transpiration stress in fully saturated soils.

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

    # Canopy Height Scaling & Atmospheric Sensitivity
    # 'stress_onset_weight' is an empirical parameter adjusting the soil moisture stress threshold 
    # based on atmospheric demand (PET) and physical vegetation stature. Higher 
    # atmospheric demand shifts the soil moisture threshold for plant stress.
    stress_onset_weight = (1 / (1 + PET)) - (canopy_buffer_sensitivity / (1 + canopy_height_meters))
    
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
        
        # Critical Moisture Point (CR)
        # Represents the specific soil moisture threshold below which the vegetation 
        # begins to actively experience transpiration reduction/moisture stress.
        CR = (1 - stress_onset_weight) * (field_capacity - WPCH) + WPCH
        
        # Transpiration Reduction Evaporative Water Stress (fTREW)
        # Calculates raw soil moisture stress. As soil moisture drops below CR, 
        # fTREW scales down toward 0. The rate of this drop-off is non-linearly 
        # driven by CHscalar as an exponent.
        fTREW = rt.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)
        
    # Relative Humidity Buffering & Final Integration
    # 'RHSM' serves as a dynamic weighting factor. It determines the balance 
    # between atmospheric control (fM) and direct soil moisture control (fTREW).
    RHSM = RH ** (4 * (1 - soil_moisture) * (1 - RH))
    
    # 'fTRM' blends the original atmospheric constraint (fM) with the new soil moisture
    # constraint (fTREW) using RHSM as the slider. NaNs in fTREW (e.g., where CR == WPCH)
    # are caught and defaulted safely to a highly-stressed zero condition.
    fTRM = (1 - RHSM) * fM + RHSM * rt.where(np.isnan(fTREW), 0, fTREW)

    return fTRM