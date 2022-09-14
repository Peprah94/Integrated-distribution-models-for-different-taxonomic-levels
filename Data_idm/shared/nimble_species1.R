source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[2]], "Spe", crossvalidate = FALSE, shared = "all")


save(estimates, file="estimate_data_species_hov1.RData")

