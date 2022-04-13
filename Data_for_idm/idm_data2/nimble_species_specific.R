load("data_idm.RData")

source("fnx_for estimation.R")




require(coda)
require(nimble)
require(devtools)
require(mcmcplots)
require(MCMCglmm)
#require(furrr)
#require(purrr)
#require(future.apply)
library(parallel)
#library(foreach)
library(progressr)
library(pbapply)

#handlers(global = TRUE)
op <- pboptions(type="timer")

run_nimble_model <- function(simulations_all){
  
  # Set up
  
  #### Load required packages ####
  require(coda)
  require(nimble)
  require(devtools)
  require(mcmcplots)
  require(MCMCglmm)
  #require(furrr)
  #require(purrr)
  #require(future.apply)
  library(parallel)
  library(foreach)
  library(ggmcmc)
  start_time <- Sys.time()
  ############################################
  code <- nimbleCode({
    
    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################
    for(spe.tag in 1:n.species){
      beta0[spe.tag]~ dnorm(0, sd=5) 
      #mean.lambda[spe.tag] <- exp(beta0[spe.tag])
      #beta0[spe.tag] <- log(mean.lambda[spe.tag])
      #mean.lambda[spe.tag] ~ dunif(0,50)
      alpha0[spe.tag] ~  dnorm(0,sd=100)
      #mean.p[spe.tag] ~ dunif(0,1)
      #beta1[spe.tag] ~ dnorm(0,sd=5)
      #alpha1[spe.tag] ~ dnorm(0,sd=100)
    }
    
    for(spe.tag in 1:n.species){
      mu.eta[spe.tag] <- 0
    }
    
    for(site.tag in 1:n.sites){
      eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],prec= omega[1:n.species, 1:n.species])
    }
    
    #Vague prior for variance covariance matrix
    omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
    Cov[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
    
    
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
        mu[site.tag,spe.tag] <- beta0[spe.tag] +eta.lam[site.tag,spe.tag]#+ beta1[spe.tag]*ecological_cov[site.tag]
        #log(lambda[site.tag, spe.tag]) <-  mu[site.tag,spe.tag]
        
        cloglog(psi[site.tag, spe.tag]) <-mu[site.tag,spe.tag]
        
        psi1[site.tag, spe.tag] <- psi[site.tag, spe.tag]-0.000001 
        #log(lambda[site.tag, spe.tag]) <-  cloglog(psi[site.tag, spe.tag])
        lambda[site.tag, spe.tag] <- -log( 1-psi1[site.tag, spe.tag])
        #cloglog(psi[site.tag, spe.tag]) <- log(lambda[site.tag, spe.tag])
        #psi[site.tag, spe.tag] <- 1-exp(- lambda[site.tag, spe.tag])
      }
    } 
    
    #Detection probability
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        logit(p.tag[site.tag, spe.tag]) <- alpha0[spe.tag] #+ alpha1[spe.tag]*detection_cov[site.tag]
      }
    }
    
    #############################################################
    #                 Occurence Model                           #
    #############################################################
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        z[site.tag,spe.tag] ~ dbern(psi[site.tag, spe.tag])
        for(k in 1:n.visit){
          X[site.tag,spe.tag,k] ~ dbin(z[site.tag,spe.tag]*p.tag[site.tag, spe.tag], n.replicates)
        }
      }
    }
    
    #############################################################
    #                pis                           #
    #############################################################


    for(site.tag in 1:n.sites){
      lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
      for(spe.tag in 1:n.species){
        pis[site.tag, spe.tag] <- (lambda[site.tag,spe.tag]/lambda.g[site.tag])
      }
    }
    #               shannon Index                              #
    #############################################################
    
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        shans[site.tag, spe.tag] <- log(pis[site.tag, spe.tag])*(pis[site.tag, spe.tag])
      }
    }
    
    for(site.tag in 1:n.sites){
      shan[site.tag] <- - sum(shans[site.tag, 1:n.species])
    }
    
    #sumLogProb ~ dnorm(0,1)
  })
  
  data <- simulations_all
  #dimensions of the data
  dim_data <- dim(data[[1]]) 
  data_dim <- dim_data[2] 
  
  n.species= (dim_data[2])
  
  # Constants for the model 
  const <- list(n.sites = dim_data[1], 
                n.species= (dim_data[2]), 
                n.replicates = 5,
                n.visit=dim_data[3],
                nu = 3,
                a.sigma = 1,
                b.sigma=1,
                NLatent=round((n.species/2),0)
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
  
  mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0", "alpha0", "shan"))
  #mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta0", "beta1", "alpha0","alpha1"))
  
  #mcmcconf$removeSamplers(c("lamLatent", "phi","sig", "delta", "a1", "a2"))
  #mcmcconf$addSampler(c("lamLatent", "phi","sig", "delta", "a1", "a2"), "RW_block")
  mcmcconf$removeSamplers(c("beta0", "alpha0"))
  mcmcconf$addSampler(c("beta0"), "RW_block")
  mcmcconf$addSampler(c("alpha0"), "RW_block")
  #mcmcconf$addSampler(c("lamLatent", "phi","sig", "delta", "a1", "a2"), "crossLevel")
  
  Rmcmc <- buildMCMC(mcmcconf)
  cmcmc <- compileNimble(Rmcmc, project = Cmwtc,resetFunctions = TRUE) 
  
  estimate <- runMCMC(cmcmc,
                      niter = 100000,
                      nburnin = 50000,
                      nchains=5,
                      thin = 10,
                      summary=TRUE,
                      samples=TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC=TRUE)
  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  
  
  mcmclist <- ggs(estimate$samples)
  print(ggs_density(mcmclist, family = "beta0"))
  print(ggs_density(mcmclist, family = "alpha0"))
  ggs_traceplot(mcmclist, family = "beta0")
  ggs_traceplot(mcmclist, family = "alpha0")
  
  
  subset_parameters <- c(paste0("beta0[",1:const$n.species,"]"),
                         paste0("beta1[",1:const$n.species,"]"),
                         paste0("alpha0[",1:const$n.species,"]"),
                         paste0("alpha1[",1:const$n.species,"]")
  )
  
  aa <-mcmclist%>%
    filter(Parameter %in% subset_parameters)%>%
    ggs_Rhat()
  Rhat_data <- aa$data[,5]
  all_rhat <- all(Rhat_data < 1.05) #Rhat is less than 1.05
  N_over_rhat <-length(which(Rhat_data > 1.05))/length(subset_parameters) #Rhats over 1.04 
  
  estimation <- estimate$summary$all.chains
  #est_fixed <- estimation[((const$n.species*const$n.species)+1):((const$n.species*const$n.species)+(const$n.species*2)),1]
  #est_pis <- matrix(estimation[,1][((const$n.species*const$n.species)+((const$n.species)*4)+1) : ((const$n.species*const$n.species)+((const$n.species)*4)+(const$n.sites*const$n.species) )], nrow=const$n.sites, ncol=const$n.species, byrow=FALSE)
  return(list(const, estimation,N_over_rhat))
}

## Run the simulations
#cl <- makeCluster(2)
#setDefaultCluster(cl)
#setDefaultCluster(cl)
#species_estimates_na <- pblapply(simulations_all, run_nimble_model, cl=cl)
#species_estimates_data <- foreach(i=1:6, .export=c('run_nimble_model'))%dopar% {
#  run_nimble_model(simulations_all[[i]])
#}
species_estimates_data <- run_nimble_model(simulations_all[[3]])
save(species_estimates_data, file="estimate_species_data.RData")
