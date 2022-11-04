#Load the required data
#load("data_idm.RData")
source("cross_valid_function.R")

require(coda)
require(nimble)
require(MCMCglmm)
library(parallel)
library(pbapply)

op <- pboptions(type="timer")

mysum <- function(x){
  # if(all(is.infinite(x)) == FALSE){
  # message("a realisation has -Inf")
  ## }else{
  ret <- sum(x[!is.infinite(x)], na.rm = TRUE)
  #}
  return(ret)
}

nimble_sum <- nimbleRcall(
  prototype = function(
    x=double(1)
    # beta is a vector
  ) {},
  returnType = double(0), # outcome is a vector
  Rfun = 'mysum'
)

#method = "IDM2"
#loss_method = "predictive"
#shared = "all"

run_nimble_model <- function(simulations_all, method, loss_method, shared, covariance_prior = "LV"){
  #### Load required packages ####
  require(coda)
  require(nimble)
  require(mcmcplots)
  require(MCMCglmm)
  require(ggmcmc)
  library(parallel)
  library(doParallel)
  
  start_time <- Sys.time()
  
  #write nimble code
  code <- nimbleCode({
    
    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################
    
    #specific specific intercept for species occupancy model
    #same for the group count model using the shared model
    for(spe.tag in 1:n.species){
      beta0[spe.tag]~ dnorm(0, sd = 100)   
    }
    
    # intercept for shared covariate and interractions
    intercept_lambda ~ dnorm(0, sd = 100) 
    
    #standard deviation of sites random effects
    sigma.obs ~ dunif(0.01,100)
    
    # covariance prior for interraction effect
    if(covariance_prior == "full"){
      omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
      Cov[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
    }
    
    if(covariance_prior == "LV"){
      delta[1] ~ dgamma(a1,1)
      for(factor in 2:NLatent) {
        delta[factor] ~ dgamma(a2,1)
      }
      
      for(l in 1:NLatent){
        tauDelta[l] <- prod(delta[1:l])
      }
      
      for(spe.tag in 1:n.species){
        sig[spe.tag] ~ dgamma(a.sigma, b.sigma)
      }
      Sigma[1:n.species, 1:n.species] <- diag(1/sig[1:n.species])
      
      for (spe.tag in 1:n.species) {
        for (l in 1:NLatent) {
          phi[spe.tag,l] ~ dgamma(nu/2,nu/2)
          tauFS[spe.tag, l] <- 1/(phi[spe.tag,l]*tauDelta[l])
          lamLatent[spe.tag, l] ~ dnorm(mean = 0, sd = sqrt(tauFS[spe.tag,l]))
        }
      }
      Cov[1:n.species, 1:n.species] <- lamLatent[1:n.species, 1:NLatent] %*% t(lamLatent[1:n.species, 1:NLatent]) + Sigma[1:n.species, 1:n.species]
      omega[1:n.species,1:n.species] <- inverse(Cov[1:n.species, 1:n.species])
    }
    
    
    for(site.tag in 1:n.sites){
      site.obs.vars[site.tag] ~ dnorm(0, sd = sigma.obs) # random site effect of observation process
      # species interraction effect
      eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],
                                              cov = Cov[1:n.species, 1:n.species])
    }  
    
    #correlation matrix
    for (spe.tag in 1:n.species) {
      for (ss in 1:n.species) {
        CorrIn[spe.tag, ss] <- Cov[spe.tag, ss]/sqrt(Cov[spe.tag, spe.tag] * Cov[ss, ss])
      }
    }
    
    #Link between the abundance and occupancy
    if(model == "shared"){
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){#loop over species
          mu[site.tag,spe.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag]
          log(lambda[site.tag, spe.tag]) <- mu[site.tag,spe.tag]
          cloglog(psi[site.tag, spe.tag]) <- mu[site.tag,spe.tag]
        } 
      }
    }
    
    if(model == "covariate"){
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){#loop over species
          mu[site.tag,spe.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] 
          cloglog(psi[site.tag, spe.tag]) <- mu[site.tag,spe.tag]
          log(lambda[site.tag, spe.tag]) <- intercept_lambda  + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] 
        }
      } 
    }
    
    if(model == "correlation"){
      #standard deviation of random effect of sites for group count data
      sigma.obs.lambda ~ dunif(0.01,100)
      #create extra site and  for lambda group count data
      for(site.tag in 1:n.sites){
        site.obs.vars.lambda[site.tag] ~ dnorm(0, sd = sigma.obs.lambda) # random site effect of observation process
      } 
      
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){#loop over species
          mu[site.tag,spe.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag]  
          cloglog(psi[site.tag, spe.tag]) <- mu[site.tag,spe.tag]
          log(lambda[site.tag, spe.tag]) <- intercept_lambda  + eta.lam[site.tag,spe.tag] + site.obs.vars.lambda[site.tag] 
        }
      } 
    }
    
    ###############################################
    #   LIKELIHOOD FOR INTEGRATED DISTRITION MODELS
    ################################################
    
    if(method == "IDM"|method == "IDM1"|method == "IDM2"){
      #Model for Group count data 
      for(site.tag in 1:n.sites){ # loop over sites
        lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
        for(year.tag in 1:n.years){ #loop over years
          for(visit.tag in 1:n.visit){ #loop over visits
            Y[site.tag, visit.tag, year.tag] ~ dpois(lambda.g[site.tag])
          }
        }
      }
      
      #   #Model for species occupancy data
      for(year.tag in 1:n.years){ # loop over years
        for(site.tag in 1:n.sites){ # loop over sites
          for(spe.tag in 1:n.species){ #loop over species
            for(visit.tag in 1:n.visit){ #loop over visits
              X[site.tag,spe.tag,visit.tag, year.tag] ~ dbin(psi[site.tag,spe.tag], n.replicates)
            }
          }
        }
      }
      
      # Relative proportion for diversity indices
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          pps[site.tag, spe.tag] <- exp(mu[site.tag, spe.tag])
          pis[site.tag, spe.tag] <- pps[site.tag, spe.tag]/sum(pps[site.tag, 1:n.species])
        }
      }
    }
    
    ###############################################
    #   LIKELIHOOD FOR GROUP COUNT  MODELS
    ################################################ 
    if(method == "IG"){
      for(site.tag in 1:n.sites){ # loop over sites
        lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
        for(year.tag in 1:n.years){ #loop over years
          for(visit.tag in 1:n.visit){ #loop over visits
            Y[site.tag, visit.tag, year.tag] ~ dpois(lambda.g[site.tag])
          }
        }
      }
      # Relative proportion for diversity indices
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          pps[site.tag, spe.tag] <- lambda[site.tag, spe.tag]
          pis[site.tag, spe.tag] <- pps[site.tag, spe.tag]/sum(pps[site.tag, 1:n.species])
        }
      }
    }
    
    ###############################################
    #   LIKELIHOOD FOR SPECIES OCCUPANCY MODELS
    ################################################ 
    if(method == "Spe"){
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){ #loop over species
          for(visit.tag in 1:n.visit){ #loop over visit
            for(year.tag in 1:n.years){ #loop over years
              X[site.tag,spe.tag,visit.tag, year.tag] ~ dbin(psi[site.tag,spe.tag], n.replicates)
            }
          }
        }
      }
      # Relative proportion for diversity indices
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          pps[site.tag, spe.tag] <- exp(mu[site.tag, spe.tag])
          pis[site.tag, spe.tag] <- pps[site.tag, spe.tag]/sum(pps[site.tag, 1:n.species])
        }
      }
    }
    
    #     #############################################################
    #     #               shannon Index                              #
    #############################################################
    
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        shans[site.tag, spe.tag] <- log(pis[site.tag, spe.tag])*(pis[site.tag, spe.tag])
      }
    }
    
    for(site.tag in 1:n.sites){
      shan[site.tag] <- - nimble_sum(shans[site.tag, 1:n.species])
    }
  })
  
  ###########################
  # DETAILS FOR THE NIMBLE MCMC
  ##########################
  
  # Extracting dimensions of data parameters
  data <- simulations_all
  dim_data <- dim(data[[1]]) 
  data_dim <- dim_data[2] 
  
  #constants for the model
  const <- list(n.sites = dim_data[1], 
                n.species= (dim_data[2]), 
                n.years = dim_data[4],
                n.replicates = 5,
                n.visit=dim_data[3],
                mu.eta = rep(0, dim_data[2]),
                nu = 3,
                NLatent=ceiling(data_dim/5),
                a1 = 30,
                a2= 30,
                a.sigma = 1,
                b.sigma = 1
  )
  
  #R and df for the inverse Wishart prior
  Rmat <- diag(const$n.species) # In the paper, this is initial covariance matrix
  df <- const$n.species + 1 
  
  #Data for the model
  idm_data <- list(Y = simulations_all[[2]], 
                   X = simulations_all[[1]], 
                   R = Rmat,
                   df = df)
  
  #For LV approach
  delta = rgamma(const$NLatent, 2,1)
  phi =matrix(rgamma(const$n.species*const$NLatent,1.5,1.5), ncol=const$NLatent, nrow=const$n.species)
  tau= cumprod(delta)
  lamLatent = matrix(NA, ncol=const$NLatent, nrow=const$n.species)
  for(spe.tag in 1:const$n.species){
    for(i in 1:const$NLatent){
      lamLatent[spe.tag,i] <- rnorm(1,0,sd=sqrt(1/(phi[spe.tag,i]*tau[i])))
    }
  }
  
  # Initial values for the model
  idm_inits <- function(){list(beta0= rnorm(const$n.species,0, 1),
                               sigma.obs = 1,
                               sigma.obs.lambda = 1,
                               intercept_lambda = rnorm(1,0,1),
                               site.obs.vars = rnorm(const$n.sites, 0, 1),
                               site.obs.vars.lambda = rnorm(const$n.sites, 0, 1),
                               eta.lam = array(0, dim = c(const$n.sites, const$n.species)),
                               omega = diag(1, const$n.species), 
                               delta=delta,
                               phi=phi,
                               lamLatent= lamLatent ,
                               sig = rgamma(const$n.species, 1,1)
  )
  }
  initsList <- idm_inits()  
  
  #############################################################
  #                 Compile Model                             #
  #############################################################
  
  mwtc <- nimbleModel(code,
                      name = "mwtc",
                      data = idm_data, 
                      constants = const, 
                      inits = initsList)
  
  Cmwtc <- compileNimble(mwtc)
  
  mcmcconf <- configureMCMC(Cmwtc, 
                            print = TRUE)
  
  #index for crossvalidation
  indx <- data$indx
  
  # Species only model
  if( method == "Spe"){
    fold_function <- function(i){
      if(i ==1 ){
        nodes <- mwtc$expandNodeNames(paste0('X[', indx[i: 37],', ,  ,]') )
      }
      if(i == 2){
        nodes <- mwtc$expandNodeNames(paste0('X[', indx[(i+36): 74],', , , ]') ) 
      }
      
      foldNodes_i <-mwtc$isData(nodes)
      ret <- nodes[foldNodes_i]
      return(ret)
    }
  }
  
  # Group count model
  if(method == "IG"){
    fold_function <- function(i){
      if(i ==1 ){
        nodes <- mwtc$expandNodeNames(paste0('Y[', indx[i: 37],',, ]') )
      }
      if(i == 2){
        nodes <- mwtc$expandNodeNames(paste0('Y[', indx[(i+36): 74],', ,]') ) 
      }
      foldNodes_i <- mwtc$isData(nodes)
      ret <- nodes[foldNodes_i]
      return(ret)
    }
  }
  
  # Unconditional crossvalidation
  if(method == "IDM"){
    fold_function <- function(i){
      if(i ==1 ){
        nodes <- mwtc$expandNodeNames(c(paste0('X[', indx[i: 37],', ,  ,]'), paste0('Y[', indx[i: 37],',, ]') ) )
      }
      if(i == 2){
        nodes <- mwtc$expandNodeNames(c(paste0('X[', indx[(i+36): 74],', , , ]'), paste0('Y[', indx[(i+36): 74],', , ]')) ) 
      }
      foldNodes_i <-mwtc$isData(nodes)
      ret <- nodes[foldNodes_i]
      return(ret)
    }
  }
  
  # crossvalidation given the group count data
  if(method == "IDM1"){
    fold_function <- function(i){
      if(i ==1 ){
        nodes <- mwtc$expandNodeNames(c(paste0('X[', indx[i: 37],', ,  ,]') ) )
      }
      if(i == 2){
        nodes <- mwtc$expandNodeNames(c(paste0('X[', indx[(i+36): 74],', , , ]')) ) 
      }
      foldNodes_i <-mwtc$isData(nodes)
      ret <- nodes[foldNodes_i]
      return(ret)
    }
  }
  
  # crossvalidation given the species occupancy data
  if(method == "IDM2"){
    fold_function <- function(i){
      if(i ==1 ){
        nodes <- mwtc$expandNodeNames(paste0('Y[', indx[i: 37],',, ]')  )
      }
      if(i == 2){
        nodes <- mwtc$expandNodeNames(paste0('Y[', indx[(i+36): 74],', , ]')) 
      }
      foldNodes_i <-mwtc$isData(nodes)
      ret <- nodes[foldNodes_i]
      return(ret)
    }
    # }
    
  }
  
  ret <- myrunCrossValidate(MCMCconfiguration = mcmcconf, 
                            k = 2,
                            nCores = 1,
                            lossFunction = loss_method,
                            MCMCcontrol = list(niter = 10000, 
                                               nburnin = 9500,
                                               inits = initsList
                            ), #5000
                            foldFunction = fold_function,
                            nBootReps = NA)
  
  ret1 <- ret$CVvalue

    return(ret)
  
  }

  