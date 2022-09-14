source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[2]], "IG", crossvalidate = FALSE, shared = "covariate_inter")


save(estimates, file="estimate_data_genus_hov.RData")

