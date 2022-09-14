load("data_idm.RData")
source("cross_validation.R")
library(doParallel)
library(parallel)
library(foreach)
library(iterators)

cl <- makeCluster(6)
#registerDoParallel(cl)
setDefaultCluster(cl)


#data <- simulations_all[[1]]
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


bumblebees <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IDM")
}, cl = cl)

save(bumblebees, file="bumblebeesIDMCV.RData")

bumblebees1 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IDM1")
}, cl = cl)

save(bumblebees1, file="bumblebeesIDM1CV.RData")

bumblebees2 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IDM2")
}, cl = cl)

save(bumblebees2, file="bumblebeesIDM2CV.RData")

bumblebees3 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "Spe")
}, cl = cl)

save(bumblebees3, file="bumblebeesSpeCV.RData")

bumblebees4 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IG")
}, cl = cl)

save(bumblebees4, file="bumblebeesIGCV.RData")

stopCluster(cl)


