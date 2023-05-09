#Script to plot the results
library(ggplot2)
library(MASS)
library(lme4)
library(rcompanion)

# Load data
load("data_idmNew.RData")
source("functionForEstimation.R")

# set theme for plots
theme_set(theme_bw())

exploratoryDataPlots <- function(i){
  #insect group
  group= list("bumblebees",  "hoverflies", "solitarybees")


groupDataPlot <- as.data.frame(cbind(ceiling(simulations_all[[i]]$mat.genus/simulations_all[[i]]$Nsurveys),
                                     simulations_all[[i]]$sites ))
colnames(groupDataPlot) <- c(paste0("visit", 1:8), "sites")

gg1 <- groupDataPlot%>%
  reshape2::melt(., id.vars = c("sites"))%>%
  #dplyr::filter(sites == "NC0429")%>%
  ggplot(mapping = aes(x = as.factor(as.numeric(variable)),
                       y = sapply(value, FUN=function(x) ifelse(x==0, -1,x)),
                       fill = variable))+
  geom_bar(stat = "identity")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~sites, ncol = 10)+
  xlab("visits")+
  ylab("counts")+
  ggtitle(paste("Distribution of ",group[i]," counts at each of the PoMS sites", sep = ""))+
  theme(text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom")


speciesDataPlot <- apply(simulations_all[[i]]$mat.species, c(1,3),
                         function(x){mean(x, na.rm = TRUE)})%>%
  as.data.frame()
speciesDataPlot <- cbind(speciesDataPlot, simulations_all[[1]]$sites )
colnames(speciesDataPlot) <- c(paste0("visit", 1:8), "sites")

gg2 <- speciesDataPlot%>%
  reshape2::melt(., id.vars = c("sites"))%>%
  #dplyr::filter(sites == "NC0429")%>%
  ggplot(mapping = aes(x = as.factor(as.numeric(variable)),
                       y = sapply(value, FUN=function(x) ifelse(x==0, -0.1,x)),
                                  fill = variable))+
  geom_bar(stat = "identity")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~sites, ncol = 10)+
  xlab("visits")+
  ylab("average occupancy probability")+
  ggtitle(paste("Distribution of ",group[i]," average occupancy prabability across all at each of the PoMS sites"))+
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "bottom")
return(list(gg1, gg2))
}

#rr <- exploratoryDataPlots(3)

#Fitting poisson regression with visits as randon variables
modelFitExploration <- function(index){
  groups = c("bumblebees",  "hoverflies", "solitarybees")

  #gr_ref <- read_csv("gr_ref.csv")
  #colnames(gr_ref)[1] <- "sites"
groupDataPlot <- as.data.frame(cbind(ceiling(simulations_all[[index]]$mat.genus/simulations_all[[index]]$Nsurveys),
                                     simulations_all[[index]]$sites,
                                     scale(simulations_all[[index]]$latitude)
                                     ))

colnames(groupDataPlot) <- c(paste0("visit", 1:8), "sites", "latitude")

#set the site with NA as 0
groupDataPlot$latitude[is.na(groupDataPlot$latitude)] <- 0

#format the data
groupDataFit <- groupDataPlot%>%
  reshape2::melt(., id.vars = c("sites", "latitude"))%>%
  dplyr::mutate(visit = as.numeric(variable))
colnames(groupDataFit) <- c("sites", "latitude","variable", "counts", "visit")

#fit intercept models to compare
#Why fitting the intercept Model and using it in comparison?
#let's write the formulas we are going to use for selecting the best models
# Poisson models for counts

pFit <- f <- aicVals <- bicVals <-  list()
f[[1]] <- as.formula(counts ~ 1)
f[[2]] <- counts ~ 1 + latitude
f[[3]] <- counts ~ latitude + (1|visit)

# fitting poisson models


for(i in 1:3){
  if(i ==3){
    pFit[[i]] <-glmer(f[[i]], family = poisson, data = groupDataFit)
    aicVals[[i]] <-  AIC(logLik(pFit[[i]]))
    bicVals[[i]] <-  BIC(logLik(pFit[[i]]))
  }else{
    pFit[[i]] <- glm(f[[i]], family = poisson, data = groupDataFit)
    aicVals[[i]] <-  AIC(logLik(pFit[[i]]))
    bicVals[[i]] <-  BIC(logLik(pFit[[i]]))
  }
}

#put all AIC and BIC values together
allAicVals <- do.call("c", aicVals)
allBicVals <- do.call("c", bicVals)

#Select the model with the smallest AIC
bestAicVals <- which.min(allAicVals)
bestBicVals <- which.min(allBicVals)

#Choose the best model
if(bestAicVals == bestAicVals){
  bestPoiModelCounts <- pFit[[bestAicVals]]
}else{
  bestPoiModelCounts <- pFit[[bestBicVals]]
}



# fitting NB models

nbFit <- aicVals <- bicVals <- thetaVals <- list()
for(i in 1:3){
  if(i ==3){
    nbFit[[i]] <- glmer.nb(f[[i]], data = groupDataFit)
    aicVals[[i]] <-  AIC(logLik(nbFit[[i]]))
    bicVals[[i]] <-  BIC(logLik(nbFit[[i]]))
    thetaVals[[i]] <- getME(nbFit[[i]], "glmer.nb.theta")
  }else{
    nbFit[[i]] <- glm.nb(f[[i]],  data = groupDataFit)
    aicVals[[i]] <-  AIC(logLik(nbFit[[i]]))
    bicVals[[i]] <-  BIC(logLik(nbFit[[i]]))
    thetaVals[[i]] <- nbFit[[i]]$theta
  }
}

#put all theta Values together
allThetaVals <- do.call("c", thetaVals)

#put all AIC and BIC values together
allAicVals <- do.call("c", aicVals)
allBicVals <- do.call("c", bicVals)

#Select the model with the smallest AIC
bestNbAicVals <- which.min(allAicVals)
bestNbBicVals <- which.min(allBicVals)

#Choose the best model for NB models
if(bestNbAicVals == bestNbAicVals){
  bestNBModelCounts <- nbFit[[bestNbAicVals]]
}else{
  bestNBModelCounts <- nbFit[[bestNbBicVals]]
}

#check if the same models were chosen for both Poisson and NB models
check1 <- bestNbAicVals == bestAicVals
check2 <- bestBicVals == bestNbBicVals

# Check for significance in overdispersion.
# Compare the best models for each one?

anova(pFit[[bestAicVals]], nbFit[[bestNbAicVals]], test = "Chisq")

pVal  <- pchisq(2 * (logLik(nbFit[[bestNbAicVals]]) - logLik(pFit[[bestNbAicVals]])),
                df = 1, lower.tail = FALSE)


if(pVal < 0.05){
  oversDisp <- TRUE
  finalCountModel <- nbFit[[bestNbAicVals]]
  indModel <- bestNbAicVals
  theta <- allThetaVals[bestNbAicVals]
  model <- "Negative Binomial"
}else{
  oversDisp <- FALSE
  finalCountModel <- pFit[[bestAicVals]]
  indModel <- bestAicVals
  theta <- NULL
  model <- "Poisson"
}

#extract Parameter values
coefBestModel <- summary(finalCountModel)$coefficients
if(indModel != 1){
intercept <- coefBestModel[1,1]
latitude <- coefBestModel[2,1]
interceptSE <- coefBestModel[1,2]
latitudeSE <- coefBestModel[2,2]
sigIntercept <- coefBestModel[1,4] < 0.05
siglatitude <- coefBestModel[1,4] < 0.05
}else{
  intercept <- coefBestModel[1,1]
  latitude <- NA
  interceptSE <- coefBestModel[1,2]
  latitudeSE <- NA
  sigIntercept <- coefBestModel[1,4] < 0.05
  siglatitude <- NA
}

if(indModel == 3){
  visitSd = as.numeric(summary(finalCountModel)$varcor)
}else{
  visitSd = NA
}
speciesSd = NA


resultsCounts <- data.frame(Model = model,
                            Data = "Counts",
                            BestModel = paste0("M",indModel-1),
                            Intercept = intercept,
                            Latitude = latitude,
                            VisitVar = visitSd,
                            SpeciesVar = speciesSd,
                            Psi = theta,
                            `P-value` = pVal)

extraResultsCounts <- data.frame(Model = model,
                                 intercept = sigIntercept,
                                 latitude = siglatitude,
                                 latitudeSE = latitudeSE,
                                 interceptSE = interceptSE,
                                 indModel = indModel,
                                 sigOverdisperse = oversDisp,
                                 check1 = check1,
                                 check2 = check2
                                 )

allResultsCounts <- list(resultsCounts = resultsCounts,
                         extraResultsCounts = extraResultsCounts,
                         finalCountModel = finalCountModel)

speciesDataPlot <- rbind(simulations_all[[index]]$mat.species[,,1],
                         simulations_all[[index]]$mat.species[,,2],
                         simulations_all[[index]]$mat.species[,,3],
                         simulations_all[[index]]$mat.species[,,4],
                         simulations_all[[index]]$mat.species[,,5],
                         simulations_all[[index]]$mat.species[,,6],
                         simulations_all[[index]]$mat.species[,,7],
                         simulations_all[[index]]$mat.species[,,8])%>%
  as.data.frame()

visit <- rep(1:8, each = nrow(simulations_all[[index]]$mat.species[,,1]))
nRep <- rep(5, nrow(simulations_all[[index]]$mat.species[,,1])*8)
speciesDataPlot <- cbind(speciesDataPlot,
                         simulations_all[[1]]$sites,
                         scale(simulations_all[[index]]$latitude),
                         visit,
                         nRep)
colnames(speciesDataPlot) <- c(paste0("species", 1:(ncol(speciesDataPlot)-4)), "sites", "latitude", "visit", "nRep")

#set the site with NA as 0
speciesDataPlot$latitude[is.na(speciesDataPlot$latitude)] <- 0

#format the data
speciesDataFit <- speciesDataPlot%>%
  reshape2::melt(., id.vars = c("sites", "latitude", "visit", "nRep"))%>%
  dplyr::mutate(species = as.numeric(variable))
colnames(speciesDataFit) <- c("sites", "latitude","visit", "nRep", "variable", "counts", "species")

binFit <- fBin <- aicBinVals <- bicBinVals <-  list()
fBin[[1]] <- as.formula(cbind(counts, nRep - counts) ~ 1)
fBin[[2]] <- as.formula(cbind(counts, nRep - counts) ~ 1 + latitude)
fBin[[3]] <- as.formula(cbind(counts, nRep - counts) ~ latitude + (1|species))
fBin[[4]] <- as.formula(cbind(counts, nRep - counts) ~ latitude + (1|visit) + (1|species))


# Fitting binomial models with cloglog link


for(i in 1:4){
  if(i > 3){
    binFit[[i]] <-glmer(fBin[[i]], family = binomial("cloglog"), data = speciesDataFit)
    aicBinVals[[i]] <-  AIC(logLik(binFit[[i]]))
    bicBinVals[[i]] <-  BIC(logLik(binFit[[i]]))
  }else{
    binFit[[i]] <- glm(fBin[[i]], family = binomial("cloglog"), data = speciesDataFit)
    aicBinVals[[i]] <-  AIC(logLik(binFit[[i]]))
    bicBinVals[[i]] <-  BIC(logLik(binFit[[i]]))
  }
}

#put all AIC and BIC values together
allBinAicVals <- do.call("c", aicBinVals)
allBinBicVals <- do.call("c", bicBinVals)

#Select the model with the smallest AIC
bestBinAicVals <- which.min(allBinAicVals)
bestBinBicVals <- which.min(allBinBicVals)

#Choose the best model
  bestBinModelOccupancy <- binFit[[bestBinBicVals]] #using BIC since it penalises for extra pars in the model
  indModel <- bestBinBicVals

  # Extract parameters
  coefBestModel <- summary( bestBinModelOccupancy)$coefficients
  if(indModel != 1){
    intercept <- coefBestModel[1,1]
    latitude <- coefBestModel[2,1]
    interceptSE <- coefBestModel[1,2]
    latitudeSE <- coefBestModel[2,2]
    sigIntercept <- coefBestModel[1,4] < 0.05
    siglatitude <- coefBestModel[1,4] < 0.05
  }else{
    intercept <- coefBestModel[1,1]
    latitude <- NA
    interceptSE <- coefBestModel[1,2]
    latitudeSE <- NA
    sigIntercept <- coefBestModel[1,4] < 0.05
    siglatitude <- NA
  }

  if(indModel == 3){
    speciesSd = as.numeric(summary(bestBinModelOccupancy)$varcor)
    visitSd = NA
  }else{
    if(indModel == 4){
      speciesSd = as.numeric(summary(bestBinModelOccupancy)$varcor)[1]
      visitSd = as.numeric(summary(bestBinModelOccupancy)$varcor)[2]
    }else{
      visitSd = NA
      speciesSd = NA
    }
  }

  model <- "Binomial"
  #groups = c("bumblebees",  "hoverflies", "solitarybees")
  resultsSpecies <- data.frame(Model = model,
                               Data = "Occupancy",
                               BestModel = paste0("M",indModel-1),
                              Intercept = intercept,
                              Latitude = latitude,
                              VisitVar = visitSd,
                              SpeciesVar = speciesSd,
                              Psi = NA,
                              `P-value` = NA)

  extraResultsSpecies <- data.frame(Model = model,
                                   intercept = sigIntercept,
                                   latitude = siglatitude,
                                   latitudeSE = latitudeSE,
                                   interceptSE = interceptSE,
                                   indModel = indModel
  )

  allResultsSpecies <- list(resultsSpecies = resultsSpecies,
                            extraResultsSpecies = extraResultsSpecies,
                           finalCountModel = bestBinModelOccupancy)


  #Put all the results for the main appendix together
allResults <- rbind(resultsSpecies,
                    resultsCounts)

return(list(species = allResultsSpecies,
            group = allResultsCounts,
            all = allResults))
}
