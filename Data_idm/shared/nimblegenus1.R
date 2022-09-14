source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[2]], "IG", crossvalidate = TRUE, shared = "all")


save(estimates, file="estimate_data_genus_hov.RData")

