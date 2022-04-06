# Integrated-distribution-models-for-different-taxonomic-levels

This is the github repository that hosts the code for the simulation study and analysis of PoMs data for the paper.

## Simulation Study
This outlines how to reproduce the simulations and estimating the parameters.

### idm_functions Folder

This folder simulates the paramters used for the simulation and hosts the functions used to simulate the data. It contains the following Rscripts:

* Parameters for simulation.R - Simulates the parameters and outputs *parameters_10.RData* and *parameters_20.RData*.
* function_for_simulation.R - Fnction that simulates the data without missing species identification.
* function_for_simulation_with_missing.R - Fnction that simulates the data with missing species identification.

### idm_miss_na10 / idm_miss_na20
Each folder simulates the data and estimates the paramters from the IDM, species only and genus only model. It uses *parameters_10.RData* and *parameters_20.RData* respectively as inputs. It contains the following Rscripts:

* simulation_interractions.R - Simulates the data, and outputs *sim_interractions_na.RData* and *sim_input_na.RData*
* nimble_interrractions.R - Estimates the paramters and diversity estimates using the IDM. It outputs *estimate_inter.RData*
* nimble_species_specific.R - Estimates the paramters and diversity estimates using the Species only model. It outputs *estimate_species.RData*
* nimble_genus_specific.R - Estimates the paramters and diversity estimates using the Genus only model. It outputs *estimate_genus.RData*
