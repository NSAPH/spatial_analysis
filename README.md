This repository contains the code used to generate the analysis in the paper "Bayesian Spatial Modeling of Mortality and Fine Particulate Air Matter Pollution in Massachusetts at the Zip Code Level". 

The files Data_processing.R produce shapefiles of MA for the analyses of the overall and stratified by age and sex populations merged with death outcomes from the Medicare Master Beneficiary Summary File (MBSF) merged with covariates from the US Census and the Behavioral Risk Factor Surveillance System and PM2.5 estimates from Di et al. (2019). 
  

Model_data_overall.R and Model_data_strat.R process the shapefiles to create a weight matrix and precompute the scaling factor for the BYM2 model and are run before the models.

The Nimble codes are in the file Nimble_code.R and the Stan codes are written in the .stan files.

The models are run in the files MA_deaths_[...].R and the plots and tables used in the paper are generated in the Plots.R file.
