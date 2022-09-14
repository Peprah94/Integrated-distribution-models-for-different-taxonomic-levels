load("~/IDM_new/idm_simulations_new/sim_interractions_na.RData")

#load the Script for analysis
source("~/IDM_new/idm_simulations_new/nimble_simulations.R")

cl <- makeCluster(10)
setDefaultCluster(cl)
clusterExport(cl, c("run_nimble_model",
                    "incidence", "nimble_incidence",
                    "richness", "nimble_richness",
                    "hill_index", "nimble_hill_index"))

clusterEvalQ(cl, {
  library(ggplot2)
  library(stringr)
  library(nimble)
  library(parallel)
  library(foreach)
  library(coda)
  library(MCMCglmm)
  library(pbapply)
  library(ggmcmc)
  library(mcmcplots)
})

species_estimates_na <- pblapply(simulations_all_na[1:100], 
                                 function(x){run_nimble_model(x, method = "Spe", covariance_prior = "full", shared = "interractions")},
                                 cl=cl)

save(species_estimates_na, file="estimate_species_na.RData")
stopCluster(cl)
