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

    Mathematical Formulations & Sequential Logic:
    ----------------------------------------------
    The algorithm initiates by establishing how atmospheric demand and physical 
    canopy structure interact to modulate plant stress. First, a dynamic stress onset 
    weight parameter ($p$, represented by the Python variable `stress_onset_weight`) is derived 
    by balancing potential evapotranspiration ($PET$) against canopy height ($CH$, 
    represented by `canopy_height_meters`), scaled by the empirical canopy weight sensitivity 
    coefficient ($a$, represented by `canopy_buffer_sensitivity`):
    $$p = \frac{1}{1 + PET} - \frac{a}{1 + CH}$$

    Simultaneously, a structural canopy height scaling factor ($CH_{scalar}$, represented by 
    `CHscalar`) is calculated as $CH_{scalar} = \sqrt{CH}$. This factor dynamically lowers 
    the effective surface soil wilting point threshold ($\theta_{WP_{CH}}$, represented by `WPCH`) 
    relative to the baseline soil-plant wilting point ($\theta_{WP}$, represented by `wilting_point`) 
    for taller canopies, mathematically reflecting their deeper root networks and superior water 
    extraction capabilities under high suction tension:
    $$\theta_{WP_{CH}} = \begin{cases} 0, & \text{if } CH_{scalar} = 0 \\ \text{clip}\left(\frac{\theta_{WP}}{CH_{scalar}}, 0, 1\right), & \text{otherwise} \end{cases}$$

    Using these boundaries, the critical soil moisture point ($\theta_{CR}$, represented by `CR`)—
    the definitive volumetric soil moisture threshold below which vegetation begins to restrict 
    transpiration—is mapped via linear interpolation between the soil field capacity 
    ($\theta_{FC}$, represented by `field_capacity`) and the canopy-scaled wilting point ($\theta_{WP_{CH}}$):
    $$\theta_{CR} = (1 - p)(\theta_{FC} - \theta_{WP_{CH}}) + \theta_{WP_{CH}}$$

    The model then calculates the raw transpiration soil moisture constraint ($f_{TREW}$ or $STREW$, 
    represented by `fTREW`) relative to this critical threshold. The rate at which water stress 
    intensifies as observed volumetric soil moisture ($\theta_{obs}$ or $VWC$, represented by `soil_moisture`) 
    drops is driven non-linearly by $CH_{scalar}$ acting as an exponent:
    $$f_{TREW} = \text{clip}\left(1 - \left(\frac{\theta_{CR} - \theta_{obs}}{\theta_{CR} - \theta_{WP_{CH}}}\right)^{CH_{scalar}}, 0, 1\right)$$

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
    * Setting to 0.0: Nullifies physical vegetation buffering. The stress onset 
      threshold becomes driven solely by atmospheric demand (PET), meaning a 
      30-meter forest and a 10-centimeter grassland respond identically to topsoil drying.
    * Setting too high (> 1.0): Overpowers the atmospheric demand term and can drive 
      stress_onset_weight ($p$) negative. This unphysically pushes the Critical Moisture 
      Point ($\theta_{CR}$) above Field Capacity ($\theta_{FC}$), causing the model to falsely simulate 
      severe transpiration stress in fully saturated soils.

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
    # 'stress_onset_weight' (p) is an empirical parameter adjusting the soil moisture stress threshold 
    # based on atmospheric demand (PET) and physical vegetation stature (CH). Higher 
    # atmospheric demand shifts the soil moisture threshold for plant stress.
    stress_onset_weight = (1 / (1 + PET)) - (canopy_buffer_sensitivity / (1 + canopy_height_meters))
    
    # 'CHscalar' (CH_scalar) acts as a proxy for aerodynamic resistance and hydraulic capacitance.
    # The square root of canopy height is a standard scaling factor in these models.
    CHscalar = np.sqrt(canopy_height_meters)
    
    # Suppress runtime warnings caused by division by zero or NaN values over 
    # invalid pixels (e.g., open water body pixels or non-vegetated regions).
    with np.errstate(divide='ignore', invalid='ignore'):
        
        # 'WPCH' (\theta_WP_CH) scales the wilting point downward for taller canopies, accounting for 
        # deeper root networks and greater water extraction capabilities under tension.
        # rt.where safely sets bare soil scenarios (CHscalar == 0) to 0.
        WPCH = rt.clip(rt.where(CHscalar == 0, 0, wilting_point / CHscalar), 0, 1)
        
        # Critical Moisture Point (CR, \theta_CR)
        # Represents the specific soil moisture threshold below which the vegetation 
        # begins to actively experience transpiration reduction/moisture stress.
        CR = (1 - stress_onset_weight) * (field_capacity - WPCH) + WPCH
        
        # Transpiration Reduction Evaporative Water Stress (fTREW, STREW)
        # Calculates raw soil moisture stress. As observed soil moisture (\theta_obs) drops below CR, 
        # fTREW scales down toward 0. The rate of this drop-off is non-linearly 
        # driven by CHscalar as an exponent.
        fTREW = rt.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)
        
    # Relative Humidity Buffering & Final Integration
    # 'RHSM' serves as a dynamic weighting factor based on RH and volumetric water content (\theta_obs).
    # It determines the balance between atmospheric control (fM) and direct soil moisture control (fTREW).
    RHSM = RH ** (4 * (1 - soil_moisture) * (1 - RH))
    
    # 'fTRM' (f_TRM) blends the original atmospheric constraint (fM) with the new soil moisture
    # constraint (fTREW) using RHSM as the slider. NaNs in fTREW (e.g., where CR == WPCH)
    # are caught and defaulted safely to a highly-stressed zero condition.
    fTRM = (1 - RHSM) * fM + RHSM * rt.where(np.isnan(fTREW), 0, fTREW)

    return fTRM