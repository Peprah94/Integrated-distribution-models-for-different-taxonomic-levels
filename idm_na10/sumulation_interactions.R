#Simulation of the data with replication at sites

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


# Packages required to run this script
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
library(pbapply)

# Simulation of covariance matrix for interraction effect
sigma_sim <- function(n){
  # S <- toeplitz((n:1)/n)
  #S <- diag(n)
  # R <- rWishart(1,(n+1),S, covariance=TRUE)
  X <- matrix(rnorm(n^2,0,0.2), n, n) # Initialize empty matrix
  m <- (t(X) + X)/2   # Make symmetric "m"
  m_lower <- m[lower.tri(m, diag = FALSE)]
  #m_lower
  sig <- runif(n,0.2,1)
  m2 <- matrix(0, n, n)
  diag(m2) <- sig
  m2[upper.tri(m2, diag = FALSE)] <- m_lower
  cov <- t(m2) %*% m2
  return(cov)
}

sim <- function(input, seed){
  #extracting inputs for the functions
  n.sites = input$constants$n.sites #no. of sites
  n.species = input$constants$n.species #no. of species
  n.visit = input$constants$n.visit #no. of visits
  n.id = input$constants$n.id 
  n.gen = input$constants$n.gen 
  n.replicate= input$constants$n.replicate #no. of replicate eg. no of traps
  q = input$constants$q #hills index exponent
  
  #Set seed
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
  
  # Simulation of interaction effect
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
  psi.s <- 1-exp(-exp(mu)) 
  
  # mean abundance
  lambda.s <- exp(mu) 
  
  #mean abundance for genus
  lambda.g <- rowSums(lambda.s) 
  
  #detection probability.
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      eta[site.tag, spe.tag]= (input$parameters$detection$alpha0)[spe.tag] + (input$parameters$detection$alpha1)[spe.tag] * (input$covariate$detection)[site.tag]
    }
  }
  p.tag <- invlogit(eta)
  
  #Simulation of observation process for the species presence-absence data
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
    y[-(sample(1:n.sites, n.gen,replace=FALSE)),k] <- NA
    for(spe.tag in 1:n.species){
      x[-(sample(1:n.sites, n.id,replace=FALSE)), spe.tag,k] <- NA    
    }
  }
  
  pis <- lambda.s/lambda.g
  
  #I estimate the hills index from the pis when performing evaluations on the parameters
  
  # Hill's number D for a given q
  #hill_number <- (rowSums((lambda.s/lambda.g)^q))^(1/(1-q))
  
  #Shannon Index
  #shan.index <- -rowSums(log(lambda.s/lambda.g)*(lambda.s/lambda.g), na.rm = TRUE)
  
  #Species richness
  #richness <- log(rowSums((lambda.s/lambda.g)^0)^(1/(1-0)))
  
  #Eveness
  #eveness <- shan.index/richness
  
  # Returning the results
  data <- list(mat.species=x,
               mat.genus = y,
               pis = pis,
               ecological_cov =input$covariate$ecological ,
               detection_cov = input$covariate$detection)
  return(data)
}


#Simulation of replicates
input10 <- list(
  constants = list(
    n.sites = 75,
    n.species = 10,
    n.visit = 3,
    n.id= 45,
    n.gen= 55,
    n.replicate = 5,
    q = 2
  ) ,
  covariate =list(
    ecological = runif(75, -1,1),
    detection = runif(75,-1,1)
  ),
  parameters = list(
    ecological=list(
      beta0 = rnorm(10,0,1),
      beta1=rnorm(10,0,1)
    ),
    detection=list(
      alpha0 = rnorm(10,0,3),
      alpha1 = rnorm(10,0.3)
    )
  ),
  interaction = sigma_sim(10)
)


input_list_na <- list(input10)


niter <- 60 #No of simulations run
simulations_all <- pblapply(input_list_na, function(x){
  pblapply(1:niter, function(z){
  data <- sim(x, seed = z)
  }, cl=4)
}, cl=4)

# 
simulations_all_na <- flatten(simulations_all)

# save the results for later use
save(simulations_all_na, file="sim_interractions_na.RData")
save(input_list_na, file="sim_input_na.RData")



