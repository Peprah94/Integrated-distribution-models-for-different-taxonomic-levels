source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[3]], "Spe", crossvalidate = FALSE, shared = "all")


save(estimates, file="estimate_data_species_sol1.RData")

