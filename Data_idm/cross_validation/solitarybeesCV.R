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

solitarybees <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IDM")
}, cl = cl)

save(solitarybees, file="solitarybeesIDMCV.RData")

solitarybees1 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IDM1")
}, cl = cl)

save(solitarybees1, file="solitarybeesIDM1CV.RData")

solitarybees2 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IDM2")
}, cl = cl)

save(solitarybees2, file="solitarybeesIDM2CV.RData")

solitarybees3 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "Spe")
}, cl = cl)

save(solitarybees3, file="solitarybeesSpeCV.RData")

solitarybees4 <- pblapply(as.list(1:6), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IG")
}, cl = cl)

save(solitarybees4, file="solitarybeesIGCV.RData")

stopCluster(cl)