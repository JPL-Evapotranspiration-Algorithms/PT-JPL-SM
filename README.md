# Priestley-Taylor Jet Propulsion Laboratory Soil Mositure (PT-JPL-SM) Evapotranspiration Model Python Implementation

[![CI](https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL-SM/actions/workflows/ci.yml/badge.svg)](https://github.com/JPL-Evapotranspiration-Algorithms/PT-JPL-SM/actions/workflows/ci.yml)

This software package is a Python implementation of the Priestley-Taylor Jet Propulsion Laboratory Soil Moisture (PT-JP-SM) model of evapotranspiration. It was re-implemented in Python by Gregory Halverson at Jet Propulsion Laboratory based on Python code developed by AJ Purdy. The original PT-JPL model was re-implemented from MATLAB code by Joshua Fisher. The PT-JPL model was designed for processing remote sensing data. It has the ability to partition latent heat flux into canopy transpiration, interception, and soil evaporation. Purdy et al., 2018 incorporated additional constraints from soil water availability on soil evaporation. Additional controls on transpiration are driven by soil water availability and canopy height include a weighting scheme to shift control on transpiration rates from soil water availability to atmospheric demand based on aridity.
 
The software was developed as part of a research grant by the NASA Research Opportunities in Space and Earth Sciences (ROSES) program. It was designed for use by the Ecosystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) mission as a precursor for the Surface Biology and Geology (SBG) mission. However, it may also be useful for general remote sensing and GIS projects in Python. This package can be utilized for remote sensing research in Jupyter notebooks and deployed for operations in data processing pipelines.
 
The software is being released according to the SPD-41 open-science requirements of NASA-funded ROSES projects.

Gregory H. Halverson (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
Lead developer<br>
NASA Jet Propulsion Laboratory 329G

Adam J. Purdy (he/him)<br>
[adpurdy@csumb.edu](mailto:adpurdy@csumb.edu)<br>
Algorithm inventor<br>
California State University Monterey Bay

Joshua B. Fisher (he/him)<br>
[jbfisher@chapman.edu](mailto:jbfisher@chapman.edu)<br>
Algorithm inventor<br>
Chapman University
 
Claire Villanueva-Weeks (she/her)<br>
[claire.s.villanueva-weeks@jpl.nasa.gov](mailto:claire.s.villanueva-weeks@jpl.nasa.gov)<br>
Code maintenance<br>
NASA Jet Propulsion Laboratory 329G
 
## Installation

```
pip install PTJPL
```

## Usage

```
import PTJPL
```
