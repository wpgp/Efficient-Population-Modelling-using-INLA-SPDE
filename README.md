# Efficient-Population-Modelling-using-INLA-SPDE
R codes for INLA-SPDE bottom-up population modelling (Cameroon Application)
Here, we provide the descriptions of the resources within this repository. It contains both the R scripts and all the relevant datasets, except the very large prediction grids covariates (prediction_data.RDS) which is available on Google Drive here: https://drive.google.com/file/d/1290hqUnBHhQS0I3iijddj34solTG-S_S/view?usp=sharing  

# R Scripts
There are three key scripts contained within the repository to help users reproduce the datasets and to facilitate the implementation of the methodology in other contexts:

1) efficient_pop_model_application_main.R - this is the main R script developed for the implementation of the methodology using real datasets based on five (5) nationally representative household listing datasets from Cameroon. It contains all the subscripts for initial data exploration, covariates selection through stepwise regression, Bayesian hierachical geostatistical model implementation, model selection, posterior simulations and grid cell predictions, uncertainty quantification, zonal statistics calculation, model cross-validation, and several plots.
   
2)  efficient_pop_sim_study.R - this contains the R codes used for the extensive simulation study conducted in the study. The simulation study explored the sensitivity of the methodology to changes in magnitudes of spatial autocorrelation and different proportions of missing data.

3) efficient_pop_sim_study_plots.R - this outlines the different codes used to produce the simulation plots within the main manuscript and the supplementary information.

# Datasets 
The input datasets used here are saved as .RData and it contains multiple datasets which are described in the respective scripts. 

Please kindly email the lead author on cc.nnanatu@soton.ac.uk or nnanatuchibuzor@gmail.com, should have any feedback or questions. This will be greatly appreciated. Thanks!
