source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[3]], "IG", crossvalidate = TRUE, shared = "covariate_inter")


save(estimates, file="estimate_data_genus_sol.RData")

