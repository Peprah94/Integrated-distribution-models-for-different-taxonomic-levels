source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[1]], "Spe", crossvalidate = FALSE, shared = "all")


save(estimates, file="estimate_data_species_bum.RData")

