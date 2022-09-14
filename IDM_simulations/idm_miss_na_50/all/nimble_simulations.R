#Runs the analyses of the simulated data
#takes input the simulated data, methods ("IDM", "Spe", "IG"), and 
#covariance_prior ("full" for sigma from invWishart), and "LV" for multiplicative
#shrinkage prior (a form of latent variable approach for sigma)
#and outputs the summary of the alpha's, beta's, z's, lambda's and correlation matrix
#as well as the values of alpha's and beta's that have converged

#
#load("/Volumes/kwakupa/IDM_new/idm_simulations_new/idm_miss_na_50/sim_interractions_na.RData")

#Packages required to run the package
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

incidence <- function(lambdas, prob, K=4){ #K is the number of visits
  incidence_mat <- lambdas*(1-(1-prob)^K)
  return(incidence_mat)
}

nimble_incidence <- nimbleRcall(
  prototype = function(
    lambda=double(2),
    prob = double(2),
    K = double(0)
    # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'incidence'
)  

richness <- function(z){
  rowSums(z)
}

nimble_richness <- nimbleRcall(
  prototype = function(
    z= double(2)
  ) {},
  returnType = double(1), # outcome is a vector
  Rfun = 'richness'
) 

hill_index <- function(q, incidence){
  pis <- proportions(incidence, margin = 1)
  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
    hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
  }
  return(hill)
}

nimble_hill_index <- nimbleRcall(
  prototype = function(
    q=double(0),
    incidence = double(2)
  ) {},
  returnType = double(1), # outcome is a vector
  Rfun = 'hill_index'
) 
run_nimble_model <- function(simulations_all, method, covariance_prior, shared){
  #### Load required packages ####
  require(coda)
  require(nimble)
  require(devtools)
  require(mcmcplots)
  require(MCMCglmm)
  require(furrr)
  require(ggmcmc)
  require(purrr)
  require(future.apply)
  library(parallel)
  library(foreach)
  
  #nimble functions for derived quantities
  incidence <- function(lambdas, prob, K=4){ #K is the number of visits
    incidence_mat <- lambdas*(1-(1-prob)^K)
    return(incidence_mat)
  }
  
  nimble_incidence <- nimbleRcall(
    prototype = function(
    lambda=double(2),
    prob = double(2),
    K = double(0)
    # beta is a vector
    ) {},
    returnType = double(2), # outcome is a vector
    Rfun = 'incidence'
  )  
  
  richness <- function(z){
    rowSums(z)
  }
  
  nimble_richness <- nimbleRcall(
    prototype = function(
    z= double(2)
    ) {},
    returnType = double(1), # outcome is a vector
    Rfun = 'richness'
  ) 
  
  hill_index <- function(q, incidence){
    pis <- proportions(incidence, margin = 1)
    if(q != 1){
      hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
    }else{
      hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
    }
    return(hill)
  }
  
  nimble_hill_index <- nimbleRcall(
    prototype = function(
    q=double(0),
    incidence = double(2)
    ) {},
    returnType = double(1), # outcome is a vector
    Rfun = 'hill_index'
  ) 
  
  ############################################
  start_time <- Sys.time()
  code <- nimbleCode({
    
    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################
    for(spe.tag in 1:n.species){
      beta0[spe.tag]~ dnorm(0, sd=3) 
      alpha0[spe.tag] ~  dnorm(0,sd=10)
      beta1[spe.tag] ~ dnorm(0,sd=3)
      alpha1[spe.tag] ~ dnorm(0,sd=10)
    }
    
    intercept_lambda ~ dnorm(0, sd = 10) # for shared covariate and interractions
    beta2 ~ dnorm(0, sd = 3)
    
    for(spe.tag in 1:n.species){
      mu.eta[spe.tag] <- 0
    }
    
    for(site.tag in 1:n.sites){
      eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],prec= omega[1:n.species, 1:n.species])
    }
    
    #Vague prior for variance covariance matrix
    if(covariance_prior == "full"){
    omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
    Cov[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
    }
    
    if(covariance_prior == "LV"){
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
      Cov[1:n.species, 1:n.species] <- lamLatent[1:n.species, 1:NLatent] %*% t(lamLatent[1:n.species, 1:NLatent]) + Sigma[1:n.species, 1:n.species]
      omega[1:n.species,1:n.species] <- inverse(Cov[1:n.species, 1:n.species])
      }
    
    #correlation matrix
    for (spe.tag in 1:n.species) {
      for (ss in 1:n.species) {
        CorrIn[spe.tag, ss] <- Cov[spe.tag, ss]/sqrt(Cov[spe.tag, spe.tag] * Cov[ss, ss])
      }
    }
    
    #Link between the abundance and occupancy
    if(shared == "all"){
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        mu[site.tag,spe.tag] <- beta0[spe.tag] + beta1[spe.tag]*ecological_cov[site.tag]+eta.lam[site.tag,spe.tag]
        cloglog(psi[site.tag, spe.tag]) <-mu[site.tag,spe.tag]
        log(lambda[site.tag, spe.tag]) <- mu[site.tag,spe.tag]  
        
      }
    } 
    }
    
    if(shared == "covariate_inter"){
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
      #mu[site.tag,spe.tag] <- beta0[spe.tag] + beta1[spe.tag]*ecological_cov[site.tag]+eta.lam[site.tag,spe.tag]
      cloglog(psi[site.tag, spe.tag]) <- beta0[spe.tag] + beta1[spe.tag]*ecological_cov[site.tag]+eta.lam[site.tag,spe.tag]
      log(lambda[site.tag, spe.tag]) <- intercept_lambda + beta1[spe.tag]*ecological_cov[site.tag]+eta.lam[site.tag,spe.tag]
        }
      }
    } 
    
    if(shared == "interractions"){
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
     # mu[site.tag,spe.tag] <- beta0[spe.tag] + beta1[spe.tag]*ecological_cov[site.tag]+eta.lam[site.tag,spe.tag]
      cloglog(psi[site.tag, spe.tag]) <- beta0[spe.tag] + beta1[spe.tag]*ecological_cov[site.tag] + eta.lam[site.tag,spe.tag]
      log(lambda[site.tag, spe.tag]) <- intercept_lambda + beta2*ecological_cov[site.tag] + eta.lam[site.tag,spe.tag]
    }
  }
} 
    
    #Detection probability
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        logit(p.tag[site.tag, spe.tag]) <- alpha0[spe.tag] + alpha1[spe.tag]*detection_cov[site.tag]
      }
    }
    
    
    # Likelihood: key definitions in the likelihood
    #############################################################
    #                 Abundance Model                           #
    #############################################################
    if(method == "IDM"){
      #Model for Group count data 
        for(site.tag in 1:n.sites){ # loop over sites
          lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
          for(visit.tag in 1:n.visit){ #loop over visits
            Y[site.tag, visit.tag] ~ dpois(lambda.g[site.tag])

        }
      }
      #   #Model for species occupancy data
        for(site.tag in 1:n.sites){ #loop over sites
          for(spe.tag in 1:n.species){ #loop over species
            z[site.tag,spe.tag] ~ dbern(psi[site.tag, spe.tag])
            for(visit.tag in 1:n.visit){ #loop over survey visits
              X[site.tag,spe.tag,visit.tag] ~ dbin(z[site.tag,spe.tag]*p.tag[site.tag, spe.tag], n.replicates)
            }
          }
        }
      
      #Derived quantities
      richness[1:n.sites] <- nimble_richness(z[1:n.sites, 1:n.species])
      hills_index0[1:n.sites] <- nimble_hill_index(0, lambda[1:n.sites, 1:n.species])
      hills_index1[1:n.sites] <- nimble_hill_index(1, lambda[1:n.sites, 1:n.species])
      hills_index2[1:n.sites] <- nimble_hill_index(2, lambda[1:n.sites, 1:n.species])
      
      }
    
    
    if(method == "IG"){
      #Model for Group counnt data 
        for(site.tag in 1:n.sites){ # loop over sites
          lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
          for(visit.tag in 1:n.visit){ #loop over visits
            Y[site.tag, visit.tag] ~ dpois(lambda.g[site.tag])
        }
        }
      
      #adding z for the parameters monitores
      #   #Model for species occupancy data. Is set to zero, but not used
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){ #loop over species
          z[site.tag,spe.tag] <- 0
        }
      }
      
      #Derived quantities
      richness[1:n.sites] <- nimble_richness(z[1:n.sites, 1:n.species])
      hills_index0[1:n.sites] <- nimble_hill_index(0, lambda[1:n.sites, 1:n.species])
      hills_index1[1:n.sites] <- nimble_hill_index(1, lambda[1:n.sites, 1:n.species])
      hills_index2[1:n.sites] <- nimble_hill_index(2, lambda[1:n.sites, 1:n.species])
    }
    
    if(method == "Spe"){
      #Model for Species occupancy data
        for(site.tag in 1:n.sites){ #loop over survey sites
          for(spe.tag in 1:n.species){ #loop over species
            z[site.tag,spe.tag] ~ dbern(psi[site.tag, spe.tag])
            for(visit.tag in 1:n.visit){ #loop over survey visits
              X[site.tag,spe.tag,visit.tag] ~ dbin(z[site.tag,spe.tag]*p.tag[site.tag, spe.tag], n.replicates)
            }
          }
        }
      
      #Derived quantities
      est_lambda[1:n.sites, 1:n.species] <- -log(1-psi[1:n.sites, 1:n.species])
      richness[1:n.sites] <- nimble_richness(z[1:n.sites, 1:n.species])
      hills_index0[1:n.sites] <- nimble_hill_index(0, est_lambda[1:n.sites, 1:n.species])
      hills_index1[1:n.sites] <- nimble_hill_index(1, est_lambda[1:n.sites, 1:n.species])
      hills_index2[1:n.sites] <- nimble_hill_index(2, est_lambda[1:n.sites, 1:n.species])
    } 
    

  })
  
  #set the data to be used
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
                NLatent=(data_dim/5),
                a.sigma = 1,
                b.sigma=1
  )
  
  
  #R and df
  Rmat <- diag(const$n.species)
  df <- const$n.species + 1 
  #omega <- LaplacesDemon:: rwishart(df, Rmat)
  
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
  
  # setting initial values
  idm_inits <- function(){list(beta0= rnorm(const$n.species,0,0.1),
                               beta1= rnorm(const$n.species,0,0.1),#covariate effect
                               beta2 = rnorm(1,0,0.1),
                               alpha0=rnorm(const$n.species,0,0.1), #detection probability
                               alpha1=rnorm(const$n.species,0,0.1),
                               intercept_lambda = rnorm(1,0,0.1),
                               z=pa_data(idm_data$X),
                               a1=2,
                               a2=3,
                               delta=delta,
                               phi=phi,
                               lamLatent= lamLatent ,
                               omega=diag(1, const$n.species),
                               sig = rgamma(const$n.species, 2,1),
                               eta.lam = array(0, dim = c(const$n.sites, const$n.species))
  )#True Presence/Absence
  }
  initsList <- idm_inits() 
  
  #Putting all together
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
  
  mcmcconf <- configureMCMC(Cmwtc,  print=TRUE,
                            monitors = c("beta0", "beta1", 
                                         "alpha0","alpha1", 
                                         "intercept_lambda",
                                         "richness", "hills_index1",
                                         "hills_index2", "Cov", 
                                         "hills_index0", "psi", "lambda"))

  Rmcmc <- buildMCMC(mcmcconf)
  cmcmc <- compileNimble(Rmcmc, project = Cmwtc,resetFunctions = TRUE) 
  
  estimate <- runMCMC(cmcmc,
                      niter = 10000,
                     nburnin = 5000,
                      inits = initsList,
                      nchains=3,
                      #thin = 50,
                      summary=TRUE,
                      samples=TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC=TRUE)
  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  
  
  mcmclist <- ggs(estimate$samples)
  aa <- mcmclist%>%
    filter(grepl('alpha|beta',Parameter))
  
  subset_parameters<- unique(aa$Parameter)
  
  # Rhat values
  subset_Rhat <- aa%>%
    ggs_Rhat()
  
  Rhat_data <- subset_Rhat$data[,5]
  all_rhat <- all(Rhat_data < 1.1) #Rhat is less than 1.05
  N_over_rhat <-length(which(Rhat_data > 1.1))/length(subset_parameters) #Rhats over 1.04 
  
  estimation <- estimate$summary$all.chains
  return(list(estimation,N_over_rhat))
}

#run_nimble_model(simulations_all_na[[1]], "IDM", "full","covariate_inter")


