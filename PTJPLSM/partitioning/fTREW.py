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
    Calculate the transpiration-side relative extractable water constraint (fTREW or STREW)[cite: 110, 184].

    This function isolates the direct surface soil moisture constraint on vegetation 
    transpiration as implemented in the PT-JPL_SM model framework[cite: 25, 150, 174].

    Parameters
    ----------
    PET : Union[Raster, np.ndarray]
        Potential evapotranspiration in watts per square meter (W/m²)[cite: 79, 197].
    canopy_height_meters : Union[Raster, np.ndarray]
        Height of the plant canopy (CH) in meters[cite: 103, 197].
    soil_moisture : Union[Raster, np.ndarray]
        Observed volumetric surface soil moisture (theta_obs or VWC) in m³/m³[cite: 167, 193].
    field_capacity : Union[Raster, np.ndarray]
        Soil field capacity (theta_FC) in m³/m³[cite: 167, 189].
    wilting_point : Union[Raster, np.ndarray]
        Baseline soil-plant wilting point (theta_WP) in m³/m³[cite: 167, 191].
    canopy_buffer_sensitivity : float, optional
        Empirical parameter (a) adjusting the weight of influence canopy height imposes 
        on the critical moisture point[cite: 197]. Must be bounded between 0.0 and 1.0.

    Returns
    -------
    Union[Raster, np.ndarray]
        fTREW stress scalar clipped to [0, 1], where undefined pixels default to 0[cite: 189].

    Mathematical Formulations & Sequential Logic:
    ----------------------------------------------
    The algorithm evaluates how atmospheric demand and physical canopy stature interact to 
    regulate plant water stress[cite: 174, 175]. Above-ground satellite-observable canopy 
    height characteristics are applied to implicitly represent plant resilience to surface soil 
    water deficits, operating under the premise that canopy height is directly related to 
    rooting depth and the potential to access water from deeper soil sources[cite: 178, 181, 182, 183].

    The calculation begins by determining a dynamic stress onset weight parameter ($p$, represented 
    by the Python variable `stress_onset_weight`) that quantifies the point at which soil water 
    availability begins to limit transpiration below the potential rate[cite: 197]. This parameter 
    balances atmospheric demand ($PET$) against canopy height ($CH$, represented by `canopy_height_meters`), 
    modulated by the empirical canopy weight sensitivity coefficient ($a$, represented by 
    `canopy_buffer_sensitivity`)[cite: 191, 197]:
    $$p = \frac{1}{1 + PET} - a\frac{1}{1 + CH}$$

    Concurrently, a structural canopy height scaling proxy ($CH_{scalar}$, represented by `CHscalar`) 
    is defined as the square root of the canopy height, which impacts the overall sensitivity to soil 
    water availability[cite: 189]:
    $$CH_{scalar} = \sqrt{CH}$$

    This canopy height scalar is used to compute the canopy height adjusted surface soil moisture 
    wilting point ($\theta_{WP_{CH}}$, represented by `WPCH`) by scaling the baseline soil-plant 
    wilting point ($\theta_{WP}$, represented by `wilting_point`) downward[cite: 191, 197]. This 
    mathematically accounts for deeper root networks and greater water extraction capabilities under 
    high tension, while safely handling non-vegetated regions or bare soil conditions via a conditional 
    assignment[cite: 189]:
    $$\theta_{WP_{CH}} = \begin{cases} 0, & \text{if } CH_{scalar} = 0 \\ \text{clip}\left(\frac{\theta_{WP}}{CH_{scalar}}, 0, 1\right), & \text{otherwise} \end{cases}$$

    Using these dynamic boundaries, the critical soil moisture point ($\theta_{CR}$, represented by 
    `CR`)—the explicit volumetric threshold below which vegetation begins to actively reduce its 
    transpiration rate—is established via linear interpolation between the soil field capacity 
    $$\theta_{CR} = (1 - p)(\theta_{FC} - \theta_{WP_{CH}}) + \theta_{WP_{CH}}$$

    Finally, the module evaluates the raw transpiration soil moisture constraint ($f_{TREW}$ or $STREW$, 
    represented by `fTREW`) by comparing the observed volumetric soil moisture ($\theta_{obs}$ or $VWC$, 
    represented by `soil_moisture`) against the calculated critical moisture threshold[cite: 186, 187, 189]. 
    The rate at which transpiration stress intensifies as soil water content drops below $\theta_{CR}$ 
    is non-linearly governed by $CH_{scalar}$ acting as an exponent[cite: 186, 187, 189]:
    $$f_{TREW} = \text{clip}\left(1 - \left(\frac{\theta_{CR} - \theta_{obs}}{\theta_{CR} - \theta_{WP_{CH}}}\right)^{CH_{scalar}}, 0, 1\right)$$

    Parameter Constraints & Eco-hydrological Bounds:
    -------------------------------------------------
    The canopy_buffer_sensitivity parameter ($a$) represents the weight of influence that canopy height 
    imposes on the critical moisture point boundary ($\theta_{CR}$)[cite: 197]. It must be physically 
    bounded within [0.0, 1.0] to preserve realistic eco-hydrological constraints:
    * Setting to 0.0: Nullifies physical vegetation buffering, driving the stress onset parameter 
      purely by atmospheric demand ($p = \frac{1}{1 + PET}$), meaning a 30-meter forest and a 
      10-centimeter grassland respond identically to topsoil drying[cite: 191].
    * Setting too high (> 1.0): Overpowers the atmospheric demand term and can drive the stress weight 
      ($p$) negative, unphysically pushing the Critical Moisture Point ($\theta_{CR}$) above Field 
      Capacity ($\theta_{FC}$) and falsely simulating severe stress in saturated soils[cite: 189, 191].

    References:
    -----------
    1. Purdy, A. J., Fisher, J. B., Goulden, M. L., Colliander, A., Halverson, G. H., 
       Tu, K., & Famiglietti, J. S. (2018). SMAP soil moisture improves global 
       evapotranspiration. Remote Sensing of Environment, 219, 1-14. 
       https://doi.org/10.1016/j.rse.2018.09.023 [cite: 2, 7, 41]

    2. Martens, B., Miralles, D. G., Lievens, H., Van Der Schalie, R., De Jeu, R. A. M., 
       Fernández-Prieto, D., Beck, H. E., Dorigo, W. A., Verhoest, N. E. C., 2017. 
       GLEAM v3: satellite-based land evaporation and root-zone soil moisture. 
       Geosci. Model Dev. 10, 1903-1925. https://doi.org/10.5194/gmd-10-1903-2017 [cite: 198, 867, 868]

    3. van Diepen, C.A., Wolf, J., van Keulen, H., Rappoldt, C., 1989. WOFOST: a simulation 
       model of crop production. Soil Use Manag. 5, 16-24. 
       https://doi.org/10.1114/j.1475-2743.1989.tb00755.x [cite: 200, 201, 943, 944]
    """
    # Parameter Validity Check
    if not (0.0 <= canopy_buffer_sensitivity <= 1.0):
        raise ValueError(
            f"Invalid canopy_buffer_sensitivity ({canopy_buffer_sensitivity}). "
            f"Parameter must be bounded between 0.0 and 1.0 to preserve "
            f"eco-hydrological physical constraints."
        )

    # 1. Stress Onset Weight Calculation (p)
    stress_onset_weight = (1 / (1 + PET)) - (canopy_buffer_sensitivity / (1 + canopy_height_meters))
    
    # 2. Structural Canopy Height Scaling (CH_scalar)
    CHscalar = np.sqrt(canopy_height_meters)

    # Suppress runtime warnings caused by division by zero over invalid pixels
    with np.errstate(divide='ignore', invalid='ignore'):
        
        # 3. Canopy-Scaled Wilting Point (\theta_WP_CH)
        WPCH = rt.clip(rt.where(CHscalar == 0, 0, wilting_point / CHscalar), 0, 1)
        
        # 4. Critical Moisture Point Determination (\theta_CR)
        CR = (1 - stress_onset_weight) * (field_capacity - WPCH) + WPCH
        
        # 5. Raw Transpiration Soil Moisture Stress Exponentiation (f_TREW / STREW)
        fTREW = rt.clip(1 - ((CR - soil_moisture) / (CR - WPCH)) ** CHscalar, 0, 1)

    # Final correction: replace undefined fTREW values with a conservative stressed condition.
    return rt.where(np.isnan(fTREW), 0, fTREW)