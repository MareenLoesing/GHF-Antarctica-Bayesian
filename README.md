# GHF-Antarctica
This is an inversion code for geothermal heat flux.

The Bayesian inversion can ivert for a set of parameters in order to fit two geotherms at certain depths. 
For the Antarctic continent, the depths of the Curie isotherm, the Moho and the lithosphere-asthenosphere boundary (LAB) (from An et al., 2015) are collected in the "Depths.csv" file.
The code reads in these depth models and calculates the temperatures at the Curie and LAB depths with the one-dimensional heat conduction equation. Output parameters are crustal and mantle thermal conductivity, heat production, and mantle heat flux.


More information will be given by Loesing et al. (2020) "Geothermal Heat Flux in Antarctica: Assessing Models and Observations by Bayesian Inversion".
