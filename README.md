# Integrated-distribution-models-for-different-taxonomic-levels

This repository is for the paper titled "IDM for data from different taxonomic levels". The scripts for the analysis are stored in different folders to make the workflow easier.

## idmFunctions

This folder hosts scripts of functions needed to simulate dataset and manipulate results from fitted models to simulated and PoMS dataset. Specifically, the scripts and what they do are:

* parameterForSimulation.R : Contains a function that simulates parameters needed to simulate data for our simulation studies.
* functionForSimulationWithMissing.R: Contains function that simulates data with a given parameter values as inputs.
* functionForEstimation: Contains functions for extracting parameters from fitted MCMC summaries, estimating mean squared error (MSE) as well as other statistics needed.

## idmSimulations

This folder hosts scripts of functions and results for simulating dataset and fit MCMC to simulated dataset.

* simulations_interractions_with_missing.R: This script simulates data using the functions defined in idmFunctions folder. It returns simInterractionsNA50.RData as the simulated dataset.
* nimbleSimulations.R: This script fits MCMC to the simulated data using NIMBLE.
* IDM.R, ..., genusCov.R: This scripts runs the function in nimbleSimulations.R for a selected model. The model type used is the first word in the name of the script.
* plotSimulations.R: This script is used to plot the results for the simulation studies.

## CaseStudy

This folder contains scripts of functions and results for the data analysis on the PoMS dataset.

* dataFormat.R : contains function to format the PoMS dataset. It outputd data_idm.RData for use in the rest of the analysis
* nimbleCaseStudy.R: script to fit the MCMC to PoMS dataset
* bumblebees.R, solitarybees.R, hoverflies.R: run the MCMC for each insect group: bumblebees, solitarybees and hoverflies respectively using the function defined in nimbleCaseStudy.R
* plotCaseStudy.R: function to plot the results from the fitted model
* postPredCheck.R: script to run posterior predictive checks for the fitted models. Uses the results from the bumblebees.R, solitarybees.R, hoverflies.R script.
* bbPostPredCheck.R, sbPostPredCheck.R, hvPostPredCheck.R: run the posterior predictive for each insect group: bumblebees, solitarybees and hoverflies respectively using the function defined in postPredCheck.R

## crossValidation

This folder hosts scripts to run the 2-fold crossValidation in the study

* _cross_valid_funcion.R_: script that formats NIMBLE's runCrossValidation script to include modifications in the log-predictive density and also to include the samplers defined in the case study.
* _cross_validation.R_: script to run the two-fold cross vaildation
* _bumblebeesCorssValidation.R, hoverfliesCorssValidation.R, solitarybeesCorssValidation.R_: run the 2-fold crossvalidation for each insect group: bumblebees, solitarybees and hoverflies respectively using the function defined in cross_validation.R

## Figures

This filder contains the figures and results of manipulation of results.

## SupplementaryMaterials1
This folder contains the exploratory analysis done on the PoMS dataset.

## SupplementaryMaterials2
This folder contains the write-up of the supplementary figures referred to in the main paper.
