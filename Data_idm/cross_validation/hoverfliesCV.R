load("data_idm.RData")
source("cross_validation.R")
library(doParallel)
library(parallel)
library(foreach)
library(iterators)

cl <- makeCluster(6)
#registerDoParallel(cl)
setDefaultCluster(cl)


#data <- 
sim <- function(i, data, method){
  #data <- simulations_all[[1]]
  #method = "IDM2"
  loss_method = as.list(rep(c("predictive", "RMSE"), each = 3))
  shared = as.list(rep(c("all", "covariate_inter", "interractions"), 2))
  ret <- run_nimble_model(simulations_all = data, method = method,
                          loss_method =  loss_method[[i]], 
                          shared = shared[[i]])
  return_value = c(loss_method[[i]], shared[[i]], method,ret)
  return(return_value)
}

clusterExport(cl, c("run_nimble_model","myrunCrossValidate",
                    "MSElossFunction", "mycalcCrossVal", "RMSElossFunction", "calcCrossValSD", 
                    "generateRandomFoldFunction", "mysum", "nimble_sum", "sim"))

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

hoverflies <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[2]], method = "IDM")
}, cl = cl)

save(hoverflies, file="hoverfliesIDMCV.RData")

hoverflies1 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[2]], method = "IDM1")
}, cl = cl)

save(hoverflies1, file="hoverfliesIDM1CV.RData")

hoverflies2 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[2]], method = "IDM2")
}, cl = cl)

save(hoverflies2, file="hoverfliesIDM2CV.RData")

hoverflies3 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[2]], method = "Spe")
}, cl = cl)

save(hoverflies3, file="hoverfliesSpeCV.RData")

hoverflies4 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[2]], method = "IG")
}, cl = cl)

save(hoverflies4, file="hoverfliesIGCV.RData")

stopCluster(cl)
