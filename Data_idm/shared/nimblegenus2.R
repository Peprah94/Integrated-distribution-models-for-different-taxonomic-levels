source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[1]], "IG", crossvalidate = TRUE, shared = "all")

save(estimates, file="estimate_data_genus_bum.RData")

