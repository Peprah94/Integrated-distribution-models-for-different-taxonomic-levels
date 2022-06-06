#Parameters needed for this function
load("/Volumes/kwakupa/idm_functions/parameters_10.RData")
source("/Volumes/kwakupa/idm_functions/function_for_simulation.R")
#Simulation of the data with replication at sites

#List of inpit parameters from idm_functions folder
input_list_na <- list(input10)

nreplicates <- 200 #Number of replicates
simulations_all <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z)
  }, cl=4)
}, cl=4)

#Put all replicates together
simulations_all_na <- flatten(simulations_all)

#save output
save(simulations_all_na, file="sim_interractions_na.RData")

#save input parameters to check if the right thing was done
#should be equal to input10.RData in the idm_functions folder
save(input_list_na, file="sim_input_na.RData")



