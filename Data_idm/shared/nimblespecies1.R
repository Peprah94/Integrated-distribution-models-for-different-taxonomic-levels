source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[2]], "Spe", crossvalidate = TRUE, shared = "all")


save(estimates, file="estimate_data_species_hov.RData")

