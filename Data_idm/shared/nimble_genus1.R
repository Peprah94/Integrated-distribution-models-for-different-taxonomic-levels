source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[2]], "IG", crossvalidate = FALSE, shared = "all")


save(estimates, file="estimate_data_genus_hov1.RData")

