#Script that hosts functions to format the MCMCoutput
#for summaries and plots
# x in all functions refers to the MCMCoutput

#Packages needed
library(unix)
require(dplyr)
require(ggcorrplot)
require(purrr)
require(ggpubr)
require(heatmaply)
require(reshape2)
library(ggplot2)
require(egg)
require(readr)
library(stringr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
# devtools::install_github("https://github.com/colinharrower/BRCmap")
library(BRCmap)
if(!require(devtools)) install.packages("devtools")
unix::rlimit_as(100*10^9)

# Extract Covariance matrix from the 
#MCMC output
extractCovariance <- function(x, #x is the MCMCoutput from the analysis
                               y #y is the species names
                               ){
  correlation_matrix <- matrix(x[[13]][grep("Cov", rownames(x[[13]])), 1], 
                               nrow=x$n.species,
                               ncol=x$n.species,
                               byrow = F)%>%
    cov2cor()
  
  rownames(correlation_matrix) <- y
  colnames(correlation_matrix) <- y
  return(correlation_matrix)
}

# Extract Shannon index from the MCMC output
extractShan <- function(x){
  shanIndex <- matrix(x[[13]][ grep("shan", rownames(x[[13]])), 1],
                       nrow= x$n.sites,
                       ncol=1,
                       byrow = F)
  return(shanIndex)
}

extractShanSD <- function(x){
  shanIndex <- matrix(1/(x[[13]][ grep("shan", rownames(x[[13]])), 3]),
                      nrow= x$n.sites,
                      ncol=1,
                      byrow = F)
  return(shanIndex)
}

#function to standard deviation of effects 
extractSigma <- function(x){
  sigma <- matrix(x[[13]][ grep("sigm", rownames(x[[13]])), 1],
                  nrow=1,
                  ncol=2,
                  byrow = T)
  return(sigma)
}

# Fundtion to extract and estimate hills indices for q = 0,1,2
hillsIndex <-function(x, q){
  pis = matrix((x[[13]][grep("lambda", rownames(x[[13]])), 1])[-1], 
               nrow=x$n.site,
               ncol=x$n.species,
               byrow = F)%>%
    proportions(., 1)
  
  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
    hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
  }
  
  return(hill)
}


# Extracting the beta0
extractBeta0 <- function(x, y){ #y is the species name
  ret<- x[[13]][ grep("beta0", rownames(x[[13]])), c(1,4,5)]%>%
    data.frame()%>%
    mutate(species_name = y)
  colnames(ret) <- c("mean", "low", "upper", "species")
  return(ret)
}

extractBeta0SD <- function(x, y){ #y is the species name
  ret<- (x[[13]][ grep("beta0", rownames(x[[13]])), 3])%>%
    data.frame()%>%
    mutate(species_name = y)
  colnames(ret) <- c("precision" , "species")
  return(ret)
}

#extracting the alpha0
extractAlpha0 <- function(x, y){
  ret<- x[[13]][ grep("alpha0", rownames(x[[13]])), c(1,4,5)]%>%
    data.frame()%>%
    mutate(species_name = y)
  colnames(ret) <- c("mean", "low", "upper", "species")
  return(ret)
}

extractAlpha0SD <- function(x, y){
  ret<- 1/(x[[13]][ grep("alpha0", rownames(x[[13]])), 3])%>%
    data.frame()%>%
    mutate(species_name = y)
  colnames(ret) <- c("precision", "species")
  return(ret)
}

#extracting the intercept of lambda
extractInterceptLambda <- function(x){
  ret<- x[[13]][ grep("intercept_lambda", rownames(x[[13]])), c(1,4,5)]
  return(ret)
}

extractInterceptLambdaSD <- function(x){
  ret<- 1/(x[[13]][ grep("intercept_lambda", rownames(x[[13]])), 3])
  return(ret)
}

#Extracting richness
extract_richness <- function(x){
  # constants <- x[[1]]
  richness <- matrix(x[[13]][ grep("richness", rownames(x[[13]])), 2],
                     nrow=x$n.sites,
                     ncol=2,
                     byrow = F)
  return(richness)
}
