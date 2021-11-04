load("sim_interractions_na.RData")

source("fnx_for estimation.R")

require(coda)
require(nimble)
require(devtools)
require(mcmcplots)
require(MCMCglmm)
require(furrr)
require(purrr)
require(future.apply)
library(parallel)
library(foreach)
library(progressr)
library(pbapply)

handlers(global = TRUE)

run_nimble_model <- function(simulations_all){
  
  # Set up
  
  #### Load required packages ####
  require(coda)
  require(nimble)
  require(devtools)
  require(mcmcplots)
  require(MCMCglmm)
  require(furrr)
  require(purrr)
  require(future.apply)
  library(parallel)
  library(foreach)
  ############################################
  code <- nimbleCode({
    
    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################
    for(spe.tag in 1:n.species){
      beta0[spe.tag]~ dnorm(0, sd=1) 
      #mean.lambda[spe.tag] <- exp(beta0[spe.tag])
      #beta0[spe.tag] <- log(mean.lambda[spe.tag])
      #mean.lambda[spe.tag] ~ dunif(0,50)
      alpha0[spe.tag] ~  dnorm(0,sd=10)
      #mean.p[spe.tag] ~ dunif(0,1)
      beta1[spe.tag] ~ dnorm(0,sd=10)
      alpha1[spe.tag] ~ dnorm(0,sd=10)
    }
    
    for(spe.tag in 1:n.species){
      mu.eta[spe.tag] <- 0
    }
    
    for(site.tag in 1:n.sites){
      eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],cov= Cov[1:n.species, 1:n.species])
    }
    
    #Vague prior for variance covariance matrix
    #omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
    #sigma2[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
    
    #Priors for covariance matrix
    a1 ~ dgamma(2,1)
    a2 ~ dgamma(2,1)
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
        lamLatent[spe.tag, l] ~ dnorm(mean = 0, sd= sqrt(tauFS[spe.tag,l]))
      }
    }
    #Cov[1:n.species, 1:n.species] <- diag(n.species)
    
    
    Cov[1:n.species, 1:n.species] <- lamLatent[1:n.species, 1:NLatent] %*% t(lamLatent[1:n.species, 1:NLatent]) + Sigma[1:n.species, 1:n.species]
    for (spe.tag in 1:n.species) {
      for (ss in 1:n.species) {
        # Cov[spe.tag, ss] <- inprod(lamLatent[1:NLatent, spe.tag], lamLatent[1:NLatent, ss])+ Sigma[spe.tag,ss]
        CorrIn[spe.tag, ss] <- Cov[spe.tag, ss]/sqrt(Cov[spe.tag, spe.tag] * Cov[ss, ss])
      }
    }
    
    precmatrix[1:n.species, 1:n.species] <- inverse(Cov[1:n.species, 1:n.species])
    #Link between the abundance and occupancy
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        mu[site.tag,spe.tag] <- beta0[spe.tag] + beta1[spe.tag]*ecological_cov[site.tag]+eta.lam[site.tag,spe.tag]
        log(lambda[site.tag, spe.tag]) <-  mu[site.tag,spe.tag]
        #log(lambda[site.tag, spe.tag]) <-  cloglog(psi[site.tag, spe.tag])
        #cloglog(psi[site.tag, spe.tag]) <-log(lambda[site.tag, spe.tag])
        #cloglog(psi[site.tag, spe.tag]) <- mu[site.tag,spe.tag]
      }
    } 
    
    #Detection probability
    #for(site.tag in 1:n.sites){
    #  for(spe.tag in 1:n.species){
    #    logit(p.tag[site.tag, spe.tag]) <- alpha0[spe.tag] + alpha1[spe.tag]*detection_cov[site.tag]
    #  }
    #}
    

    # Likelihood: key definitions in the likelihood
    #############################################################
    #                 Abundance Model                           #
    #############################################################
    
    for(site.tag in 1:n.sites){
      lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
      for(k in 1:n.visit){
        Y[site.tag,k] ~ dpois(lambda.g[site.tag])
      }
    }
    
    #############################################################
    #                pis                           #
    #############################################################
    
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        pis[site.tag, spe.tag] <- (lambda[site.tag,spe.tag]/lambda.g[site.tag])
      }
    }
    #sumLogProb ~ dnorm(0,1)
  })
  
  data <- simulations_all
  #dimensions of the data
  dim_data <- dim(data[[1]]) 
  data_dim <- dim_data[2] 
  
  
  
  # Constants for the model 
  const <- list(n.sites = dim_data[1], 
                n.species= (dim_data[2]), 
                n.replicates = 5,
                n.visit=dim_data[3],
                nu = 3,
                NLatent=2,
                a.sigma = 1,
                b.sigma=1
  )
  
  
  #R and df
  Rmat <- diag(const$n.species)
  df <- const$n.species + 1 
  #Data for the model
  idm_data <- list(Y = simulations_all[[2]], 
                   X = simulations_all[[1]], 
                   detection_cov= data$detection_cov,
                   ecological_cov = data$ecological_cov,
                   R = Rmat,
                   df = df)
  
  #Retrieving the true Presence/Absence
  pa_data <- function(data){
    zst <- apply(data, c(1,2), max)
    zst[zst!= 0] <- 1
    zst[is.na(zst)] <- 0
    return(zst)
  }
  delta = rgamma(const$NLatent, 2,1)
  phi =matrix(rgamma(const$n.species*const$NLatent,1.5,1.5), ncol=const$NLatent, nrow=const$n.species)
  tau= cumprod(delta)
  lamLatent = matrix(NA, ncol=const$NLatent, nrow=const$n.species)
  
  for(spe.tag in 1:const$n.species){
    for(i in 1:const$NLatent){
      lamLatent[spe.tag,i] <- rnorm(1,0,sd=sqrt(1/(phi[spe.tag,i]*tau[i])))
    }
  }
  idm_inits <- function(){list(beta0= rnorm(const$n.species,0,1),
                               beta1= rnorm(const$n.species,0,1),#covariate effect
                               alpha0=rnorm(const$n.species,0,1), #detection probability
                               alpha1=rnorm(const$n.species,0,1),
                               #mean.p = rep(0.5,const$n.species),
                               omega= diag(const$n.species),
                               z=pa_data(idm_data$X),
                               a1=2,
                               a2=2,
                               delta=delta,
                               phi=phi,
                               lamLatent= lamLatent ,
                               sig = rgamma(const$n.species, 2,1)
  )#True Presence/Absence
  }
  initsList <- idm_inits() 
  
  modelInfo <- list(
    code = code,
    constants = const,
    data = idm_data,
    inits = idm_inits
  )
  
  #############################################################
  #                 Compile Model                             #
  #############################################################
  
  mwtc <- nimbleModel(code,
                      name = "mwtc",
                      data = idm_data, 
                      constants = const, 
                      inits = initsList)
  Cmwtc <- compileNimble(mwtc)
  
  mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0", "beta1", "alpha0","alpha1", "pis","CorrIn"))
  #mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0", "beta1", "alpha0","alpha1"))
  
  mcmcconf$removeSamplers(c("lamLatent", "phi","sig", "delta", "a1", "a2"))
  mcmcconf$addSampler(c("lamLatent", "phi","sig", "delta", "a1", "a2"), "RW_block")
  #mcmcconf$addSampler(c("lamLatent", "phi","sig", "delta", "a1", "a2"), "crossLevel")
  
  Rmcmc <- buildMCMC(mcmcconf)
  cmcmc <- compileNimble(Rmcmc, project = Cmwtc,resetFunctions = TRUE) 
  
  estimate <- runMCMC(cmcmc,
                      niter = 500000,
                      nburnin = 400000,
                      nchains=4,
                      thin = 10,
                      summary=TRUE,
                      samples=TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC=TRUE)
  
  estimation <- estimate$summary$all.chains
  est_fixed <- estimation[((const$n.species*const$n.species)+1):((const$n.species*const$n.species)+40),1]
  est_pis <- matrix(estimation[,1][((const$n.species*const$n.species)+((const$n.species)*4)+1) : ((const$n.species*const$n.species)+((const$n.species)*4)+(const$n.sites*const$n.species) )], nrow=const$n.sites, ncol=const$n.species, byrow=FALSE)
  return(list(est_pis, est_fixed, estimation))
}

## Run the simulations
#run_nimble_model(simulations_all)
cl <- makeCluster(5)
setDefaultCluster(cl)
genus_estimates_na <- pblapply(simulations_all_na, run_nimble_model, cl=cl)

save(genus_estimates_na, file="estimate_genus_na.RData")
