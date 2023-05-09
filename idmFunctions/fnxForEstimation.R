#Hosts all the functions needed for post estimation tricks

library(unix)
unix::rlimit_as(100*10^9)

#define my cloglog function to be used for simulation
myclog <- function(psi){
  return(-log(1-psi))
}

myinvclog <- function(mu){
  return(1-exp(-exp(mu)))
}


#mean square error
mse <- function(true, estimate){
  ret <- mean((true-estimate)^2)
  return(ret)
}

###########################
# Diversity indices
#######################
incidence <- function(lambdas, prob, K=4){ #K is the number of visits
  incidence <- lambdas*(1-(1-prob)^K)
  return(incidence)
}


hill_index <- function(q, incidence){
  pis <- proportions(incidence, margin = 1)
  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
      hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
      }
  return(hill)
}

hill_index_pis <- function(q, pis){
  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
    hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
  }
  return(hill)
}

shan_index <- function(pis){
 shan <-  -rowSums(log(pis)*(pis), na.rm = TRUE)
  return(shan)
}

#Estimation of eveness
eveness <- function(pis){
  shan <- shan_index(pis)
  richness <- log(hill_index(0,pis))
  even <- shan/richness
return(even)
  }

hill_index_alt <- function(q, lambdas, zs){
  pis <- proportions(lambdas * zs, margin = 1)
  if(q != 0 & q != 1){
    hill <- rowSums((pis^q))^(1/(1-q)) 
  }else{
    if(q == 0){
      hill <- rowSums(zs)
    }else{if(q == 1){
      hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
    }
    }
    
  }
  return(hill)
}

subsetting_parameters <- function(x,name, method){
  if(method %in% c("mean", "sd")){
  if(method == "mean"){
  ret <- x[grep(name, rownames(x)),1]
  }
  if(method == "sd"){
    ret <- x[grep(name, rownames(x)),3]
  }
  }else{
      stop("method should have name mean or sd")
    }
  return(ret)
}


