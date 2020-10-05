# spatialpower
final code package for spatial power analysis and tissue generation paper. name subject to change. 


## Overview
This is a package for _in silico_ tissue generation and visualization. PackageName is useful for generating datasets for _in silico_ ground truthing, generation of spatial null models, testing of new algorithms for spatially-resolved single cell biology, and power analysis.

## Dependencies
PackageName is written for Python 3.7.3. 

It requires the following dependencies: 

- `numpy == 1.17.3`
- `matplotlib == 3.1.2`
- `pandas == 0.25.3`
- `networkx == 2.4`
- `joblib == 0.15.1`
- `scipy == 1.4.1`
- `jax == 0.1.62`

Example noteboooks make invoke additional dependencies, but are not required for the base functionality of PackageName. 

## Installation
1. Install appropriate dependencies. 
2. `pip install packageName` or `conda install -c [dev] packageName` 

We also make a Docker image available here. 

## Usage

We provide several notebooks as examples on potential uses of PackageName:
1. Tissue generation
2. Tissue labeling from existing scaffold
3. Using _in silico_ tissue for spatial power analysis [Citation]

PackageName can also be called from the command line for generation of tissues:
`example command goes here`. 

## Citing this code
If you use PackageName in your work, please cite: 

## Licence 

Copyright (C) 2020  Ethan Alexander García Baker

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

