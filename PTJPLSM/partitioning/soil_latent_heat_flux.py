from typing import Union
import numpy as np
import rasters as rt
from rasters import Raster

import PTJPL
from ..constants import PT_ALPHA

def calculate_soil_latent_heat_flux(
        Rn_soil: Union[Raster, np.ndarray], 
        G: Union[Raster, np.ndarray], 
        epsilon: Union[Raster, np.ndarray], 
        fwet: Union[Raster, np.ndarray], 
        fREW: Union[Raster, np.ndarray], 
        PT_alpha: float = PT_ALPHA) -> Union[Raster, np.ndarray]:
    """
    Calculate the Priestley-Taylor Jet Propulsion Laboratory Soil Moisture (PT-JPL-SM) soil latent heat flux.

    Parameters:
    Rn_soil (Union[Raster, np.ndarray]): The soil net radiation in watts per square meter (W/m^2). It represents the balance of incoming and outgoing radiation.
    G (Union[Raster, np.ndarray]): The soil heat flux in watts per square meter (W/m^2). It represents the amount of heat moving into or out of the soil.
    epsilon (Union[Raster, np.ndarray]): The ratio of the slope of the saturation vapor pressure curve to the sum of the slope and psychrometric constant (delta / (delta + gamma)).
    fwet (Union[Raster, np.ndarray]): The relative surface wetness. It ranges from 0 (completely dry) to 1 (completely wet).
    fREW (Union[Raster, np.ndarray]): Fraction of relative extractable water, replacing fSM.
    PT_alpha (float, optional): The Priestley-Taylor alpha constant. It is a dimensionless empirical parameter typically equal to 1.26. Defaults to PT_ALPHA.

    Returns:
    Union[Raster, np.ndarray]: The calculated soil latent heat flux in watts per square meter (W/m^2). It represents the energy associated with the phase change of water in the soil.
    """
    return PTJPL.calculate_soil_latent_heat_flux(
        Rn_soil=Rn_soil, 
        G=G, 
        epsilon=epsilon, 
        fwet=fwet, 
        fSM=fREW, 
        PT_alpha=PT_alpha
    )
