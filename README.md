# GHF-Antarctica
This is an inversion code for geothermal heat flux.

The Bayesian inversion can invert for a set of parameters in order to fit two geotherms at certain depths. 
For the Antarctic continent, the depths of the Curie isotherm, the Moho and the lithosphere-asthenosphere boundary (LAB) (from An et al., 2015, https://doi.org/10.1002/2014JB011332 and https://doi.org/10.1002/2015JB011917) are collected in the "Depths.csv" file.
The code reads in these depth models and calculates the temperatures at the Curie and LAB depths with the one-dimensional heat conduction equation. Output parameters are crustal and mantle thermal conductivity, heat production, and mantle heat flux.
Heat production values for the Antarctic Peninsula (Peninsula_HP.csv) are from Burton-Johnson et al. 2017 (https://doi.org/10.1002/2017GL073596).

Get started by running do_grid.py. The toy.py module performs the Markov Chain Monte Carlo algorithm and Tpy.py calculates the isotherm temperatures. Initial code and toy.py were created by Wolfgang Szwillus (https://github.com/wokos).

More information will be given by Loesing et al. (2020) "Geothermal Heat Flux in Antarctica: Assessing Models and Observations by Bayesian Inversion". https://doi.org/10.3389/feart.2020.00105
