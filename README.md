# GHF-Antarctica
This is an inversion code for geothermal heat flux.
This Bayesian inversion can ivert for a set of parameters in order to fit certain geotherms at certain depths. 
For the antarctic continent, the depths of the Curie isotherm, the Moho and the lithosphere-asthenosphere boundary (LAB) (from An et al., 2015) are collected in the "Depths.txt" file.
The code reads in these depths and calculates the temperatures with the one-dimensional heat conduction equation at the Curie and LAB depths.
