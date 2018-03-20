# SpaceTimeClusteringDiseaseTrends

This repository contains R code for a Bayesian space-time model for clustering areal units based on their disease risk trends. This code is specifically for the Poisson variant of the model. An additional R file containing code on how to implement the model is also provided, with accompanying data relating to respiratory hospital admissions in the Greater Glasgow & Clyde region of Scotland between 2002 and 2011.

Descriptions of the available files are as follows:

CaseStudy2RCodeExample.R
Here, the example code for the respiratory hospital admissions in Greater Glasgow & Clyde between 2002 and 2011 can be found, which shows how to the data should be structured in order to use the model, as well as comments on how to choose which trends to include in the model. The code also shows how to obtain the necessary results from the model and plots of the data.

Poisson.LTM.General.R
This file contains the main function to implement the Poisson variant of the Bayesian space-time model for clustering areal units based on their trends in disease risk.

Poisson.LTM.General.cpp
This file contains accompanying functions used by the main function implemented using Rcpp for Metropolis updates used by the model.

respiratory admissions.csv and respiratory admissions expected counts.csv
These files contain the data for the respiratory hospital admissions data in Greater Glasgow & Clyde between 2002 and 2011.

GlasgowIZ.csv, Scotland study area.shp and Scotland study area.dbf
These files contain the intermediate zone codes for areas comprising Greater Glasgow & Clyde, as well as the shape files for mapping the data.

