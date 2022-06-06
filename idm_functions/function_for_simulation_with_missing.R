#Simulation of the data with replication at sites
#Function needed to run the simulation
#With missing species ID

##############################################################################################
#                         NOTES                                                              #
##############################################################################################
### We assume n.genus  g = 1,2,...,n.genus
# For each g genera, we assume n.species in that geunus
# We assume the data was collected at n.sites
#We also assume the total number of visits to each site is n.visit.
# We also assume that  some of the data are collected with species ID
# Most of them are collected with genus ID
#n.id = number of sites that have ID at the species level
#n.gen = number of sites that have ID at the genus level
#shan.index is the shannon's index
#gamma is the detection probability for the abundance model (genus data)
#p.tag is the detectiion probability for the occurence model (species data)
#beta is a n.species* nsites matrix
#includecov (T/F) indicates whether it is a model with a covariate (TRUE) or no covariate (FALSE)
##############################################################################################
set.seed(2020)
myclog <- function(psi){
  return(-log(1-psi))
}

myinvclog <- function(mu){
  return(1-exp(-exp(mu)))
}

#Packages needed to run this script
require(coda)
require(nimble)
require(devtools)
require(mcmcplots)
require(MCMCglmm)
require(furrr)
require(purrr)
require(dplyr)
require(stringr)
require(MASS)
require("LaplacesDemon")
#library(rWishart)
library(pbapply)

#simulating the covariance matrix 


sim <- function(input, seed){
  n.sites = input$constants$n.sites
  n.species = input$constants$n.species
  n.visit = input$constants$n.visit
  n.id = input$constants$n.id
  n.gen = input$constants$n.gen
  n.replicate= input$constants$n.replicate
  q = input$constants$q
  
  set.seed(seed)
  
  #### CHECKS
  
  # is anything missing
  
  if(is.null(n.sites)){stop("please specify number of sites (n.sites)")}
  if(is.null(n.species)){stop("please specify number of species (n.species)")}
  if(is.null(n.visit)){stop("please specify number of visits (n.visit)")}
  if(is.null(n.id)){stop("please specify number of sites that have ID at the species level (n.id)")}
  if(is.null(n.gen)){stop("please specify number of species that have ID at the species level(n.gen)")}
  #if(is.null(cov)){stop("please provide covariates (n.cov)")}
  if(is.null(n.replicate)){stop("please number of replicates for the observation process")}
  
  # are all the elements of the dimensions correct?
  
  if((length(n.sites) == length(n.species)) == FALSE)
  {stop("inputs(n.sites, n.species, n.visit) are not of the same length")}
  
  if((length(n.id) == length(n.gen)) == FALSE)
  {stop("inputs( n.id, n.gen) are not of the same length")}
  
  if((length(input$parameters$ecological$beta0) !=n.species) & (length(input$parameters$ecological$beta1)!= n.species))
  {stop("length of ecological parameters = number of species")}
  
  if((length(input$parameters$detection$alpha0) !=n.species) & (length(input$parameters$detection$alpha1)!= n.species))
  {stop("length of detection parameters = number of species")}
  
  if(dim(input$interaction)[1] != dim(input$interaction)[1])
  {stop("Dimension of interaction variance must be n.species * n.species")}
  
  # Dimension of matrix to store the various parameter values
  
  N <-log_lambda <- vector("numeric", n.sites)
  z <-mu <- lambda.s <- epsilon <- eta <- p.tag <- matrix(NA, nrow = n.sites, ncol = n.species)
  x <- array(NA, dim = c(n.sites, n.species, n.visit))
  y <- matrix(NA, nrow=n.sites, ncol=n.visit)
  
  # Simulation of interaction
  for(site.tag in 1:n.sites){
    epsilon[site.tag,] = mvrnorm(1,rep(0, n.species),input$interaction )
  }
  
  
  #### Simulation of data 
  #linear predictor for ecological process
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      mu[site.tag, spe.tag] = (input$parameters$ecological$beta0)[spe.tag] + (input$parameters$ecological$beta1)[spe.tag]* (input$covariate$ecological)[site.tag] + epsilon[site.tag,spe.tag]
    }
  }
  
  # estimation of occupancy probability
  psi.s <- 1-exp(-exp(mu)) - 0.000001
  # Subtracting 0.000001 ensuresthat psi.s gets very large but is never 1.
  #if psi.s = 1, then lambda.s is inf and does not generate a sample
  
  # mean abundance
  lambda.s <- -log(1-psi.s)
  
  #mean abundance for genus
  lambda.g <- rowSums(lambda.s) 
  
  #detection probability.
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      eta[site.tag, spe.tag]= (input$parameters$detection$alpha0)[spe.tag] + (input$parameters$detection$alpha1)[spe.tag] * (input$covariate$detection)[site.tag]
    }
  }
  p.tag <- invlogit(eta)
  
  #Simulation of observation process
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      #True presence or absence
      z[site.tag, spe.tag] <- rbinom(1,1, psi.s[site.tag, spe.tag]) 
      for(k in 1:n.visit){
        # Genus counts
        y[site.tag,k] <- rpois(1, lambda.g[site.tag])
        # Species Presence and absence observations
        x[site.tag, spe.tag,k] <- rbinom(1,n.replicate, z[site.tag, spe.tag]*p.tag[site.tag, spe.tag])
      }
    }
  }
  
  # Making provision for the missing data at sites and species
  
  for(k in 1:n.visit){
    index_missing_genus <- sample(1:n.sites, n.gen,replace=FALSE)
    index_missing_site_species <- sample(1:n.sites, n.id,replace=FALSE)
    y[-(index_missing_genus),k] <- 0
    for(site.tag in seq_along(index_missing_site_species)){
      index_missing_species <- sample(1:n.species, input$constants$nspecies_no_missing, replace = FALSE)
    x[(index_missing_site_species[site.tag]),-(index_missing_species) ,k] <- 0
    }
  }

  
  #proportions for shannon index
  pis <- lambda.s/lambda.g
  
  
  # Returning the results
  data <- list(mat.species=x,
               mat.genus = y,
               pis = pis,
               ecological_cov =input$covariate$ecological ,
               detection_cov = input$covariate$detection)
  return(data)
}
