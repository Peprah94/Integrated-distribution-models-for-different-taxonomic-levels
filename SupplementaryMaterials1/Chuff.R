#Script to plot the results
library(ggplot2)
library(MASS)
library(lme4)
library(rcompanion)


SimCountProblem <- function(lambda, nsim=1e2) {
  x <- rpois(nsim, lambda)
  x0 <- x>0.5
  glm.c <- summary(glm(x ~ 1, family = poisson()))$coefficients
  glm.b <- summary(glm(x0 ~ 1, family = binomial("cloglog")))$coefficients

  r <- glm.c["(Intercept)", "Std. Error"]/glm.b["(Intercept)", "Std. Error"]

  c(lambda=lambda, p = mean(x0), r=r)
}

lambdas <- rep((1:20)/10, each=50)
Ratio <- sapply(lambdas, SimCountProblem)

plot(Ratio["lambda",], Ratio["r",])
plot(Ratio["p",], Ratio["r",])





# Load data
load("data_idm.RData")
source("functionForEstimation.R")

# set theme for plots
theme_set(theme_bw())

exploratoryDataPlots <- function(i){
  #insect group
  group= list("bumblebees",  "hoverflies", "solitarybees")


  groupDataPlot <- as.data.frame(cbind(simulations_all[[i]]$mat.genus,
                                       simulations_all[[i]]$sites ))
  colnames(groupDataPlot) <- c(paste0("visit", 1:8), "sites")

  gg1 <- groupDataPlot%>%
    reshape2::melt(., id.vars = c("sites"))%>%
    #dplyr::filter(sites == "NC0429")%>%
    ggplot(mapping = aes(x = as.factor(as.numeric(variable)), y = value, fill = variable))+
    geom_bar(stat = "identity")+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~sites, ncol = 10)+
    xlab("visits")+
    ylab("counts")+
    ggtitle(paste("Distribution of ",group[i]," counts at each of the PoMS sites", sep = ""))+
    theme(text = element_text(size = 20),
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
    ggplot(mapping = aes(x = as.factor(as.numeric(variable)), y = value, fill = variable))+
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
modelFitExploration <- function(i){
  group= list("bumblebees",  "hoverflies", "solitarybees")
  groupDataPlot <- as.data.frame(cbind(simulations_all[[i]]$mat.genus,
                                       simulations_all[[i]]$sites,
                                       scale(simulations_all[[i]]$latitude)
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

  pVal  <- pchisq(2 * (logLik(nbFit[[bestNbAicVals]]) - logLik(pFit[[bestAicVals]])),
                  df = 1, lower.tail = FALSE)


  if(pVal < 0.05){
    oversDisp <- TRUE
    finalCountModel <- nbFit[[bestNbAicVals]]
    indModel <- bestNbAicVals
    theta <- allThetaVals[bestNbAicVals]
  }else{
    oversDisp <- FALSE
    finalCountModel <- pFit[[bestAicVals]]
    indModel <- bestAicVals
    theta <- NULL
  }




  pFit <- glm(counts ~1, family = poisson, data = groupDataFit)
  pGLM <- summary(pFit)$coefficients
  nbfit <- glm.nb(counts ~1, data = groupDataFit)
  nbGLM <- summary(nbfit)$coefficients
  bGLM <- summary(glm((counts>0) ~1, family = binomial("cloglog"), data = groupDataFit))$coefficients

  r <- pGLM["(Intercept)", "Std. Error"]/bGLM["(Intercept)", "Std. Error"]
  rnb <- pGLM["(Intercept)", "Std. Error"]/nbGLM["(Intercept)", "Std. Error"]
  theta <-  nbfit$theta


  # Checking overparameterisation
  # https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
  intPval  <- pchisq(2 * (logLik(nbfit) - logLik(pFit)), df = 1, lower.tail = FALSE)
  #fit full models and do model comparisons
  m.poi <- glmer(counts ~ latitude + (1|visit), family = poisson, data = groupDataFit)
  pGLM.re <- summary(m.poi)$coefficients
  m.nb <- glmer.nb(counts ~ latitude + (1|visit), data = groupDataFit)
  nbGLM.re <- summary(m.nb)$coefficients
  b.nb <- glmer((counts>0) ~ latitude + (1|visit), family = binomial("cloglog"), data = groupDataFit)
  bGLM.re <- summary(b.nb)$coefficients

  theta.re <- getME(m.nb, "glmer.nb.theta")

  #ration of Poisson and Binomial using cloglog
  r.re <- pGLM.re["(Intercept)", "Std. Error"]/bGLM.re["(Intercept)", "Std. Error"]
  rlatitude.re <- pGLM.re["latitude", "Std. Error"]/bGLM.re["latitude", "Std. Error"]

  #ration of Poisson and NB
  rnb.re <- pGLM.re["(Intercept)", "Std. Error"]/nbGLM.re["(Intercept)", "Std. Error"]
  rnblatitude.re <- pGLM.re["latitude", "Std. Error"]/nbGLM.re["latitude", "Std. Error"]

  # model assumption for random effects
  intPval.re  <- pchisq(2 * (logLik(m.nb) - logLik(m.poi)), df = 1, lower.tail = FALSE)

  #variance of random effects
  nbVisitSd = as.numeric(summary(m.nb)$varcor)
  poiVisitSd = as.numeric(summary(m.poi)$varcor)
  binVisitSd = as.numeric(summary(b.nb)$varcor)

  model <- c("Poisson", "Negative Binomial","Poisson GLMM", "Negative Binomial GLMM")
  intercept <- c(pGLM[1] , nbGLM[1], pGLM.re[1], nbGLM.re[1])
  latitude <- c(NA, NA,pGLM.re[2], nbGLM.re[2])
  intSd <- c(NA, rnb, NA, rlatitude.re)
  latSd <- c(NA, NA, pGLM.re["latitude", "Std. Error"], nbGLM.re["latitude", "Std. Error"])
  thetaEst <- c(NA, theta, NA, rlatitude.re)
  visitSd <- c(NA, NA, poiVisitSd[1], nbVisitSd[1])
  pVal <- c(NA, intPval, NA, intPval.re )

  resultsCounts <- data.frame(Model = model,
                              Intercept = intercept,
                              Latitude = latitude,
                              `Standard Error of Intercept` = intSd,
                              `Standard Error of latitude` = latSd,
                              `Overdispersion parameter` = thetaEst,
                              `Standard Error of visit effect` = visitSd,
                              `P-value` = pVal)

  # resultsCounts <- data.frame(intRatio = r,
  #                       intNBRatio = rnb,
  #                       intRatioRE = r.re,
  #                       intNBRatioRE = rlatitude.re,
  #                       theta = theta,
  #                       thetaRE = rlatitude.re,
  #                       intPoi = pGLM[1],
  #                       intNB = nbGLM[1],
  #                       intBN = bGLM[1],
  #                       intPoiRe = pGLM.re[1],
  #                       intNbRe = nbGLM.re[1],
  #                       intBnRe = bGLM.re[1],
  #                       nbVisitSd = nbVisitSd[1],
  #                       poiVisitSd = poiVisitSd[1],
  #                       binVisitSd = binVisitSd[1],
  #                       group = group[i])

  speciesDataPlot <- rbind(simulations_all[[i]]$mat.species[,,1],
                           simulations_all[[i]]$mat.species[,,2],
                           simulations_all[[i]]$mat.species[,,3],
                           simulations_all[[i]]$mat.species[,,4],
                           simulations_all[[i]]$mat.species[,,5],
                           simulations_all[[i]]$mat.species[,,6],
                           simulations_all[[i]]$mat.species[,,7],
                           simulations_all[[i]]$mat.species[,,8])%>%
    as.data.frame()

  visit <- rep(1:8, each = nrow(simulations_all[[i]]$mat.species[,,1]))
  nRep <- rep(5, nrow(simulations_all[[i]]$mat.species[,,1])*8)
  speciesDataPlot <- cbind(speciesDataPlot,
                           simulations_all[[1]]$sites,
                           scale(simulations_all[[i]]$latitude),
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


  #Fit the models and compare
  speciesGLM <- summary(glm(cbind(counts, nRep - counts) ~1, family = binomial("cloglog"), data = speciesDataFit))$coefficients
  cloglogGLMM <- glmer(cbind(counts, nRep - counts) ~ latitude + (1|species) + (1|visit), family = binomial("cloglog"), data = speciesDataFit)
  speciesGLMM <-  summary(cloglogGLMM)$coefficients
  speciesVarGLMM <- as.numeric(summary(cloglogGLMM)$varcor)

  # resultsSpecies <- data.frame(intOnly = speciesGLM[,1],
  #                              intOnlysd = speciesGLM[,2],
  #                              glmmInt = speciesGLMM[1,1],
  #                              glmmLat = speciesGLMM[2,1],
  #                              glmmIntsd = speciesGLMM[1,2],
  #                              glmmLatsd = speciesGLMM[2,2],
  #                              speciesSD = speciesVarGLMM[1],
  #                              visitSD = speciesVarGLMM[2],
  #                              group = group[i])


  model <- c("Binomial", "Binomial GLMM")
  intercept <- c(speciesGLM[,1] , speciesGLMM[1,1])
  latitude <- c(NA,speciesGLMM[2,1])
  intSd <- c(speciesGLM[,2], speciesGLMM[1,2])
  latSd <- c(NA, speciesGLMM[2,2])
  visitSd <- c(NA, speciesVarGLMM[2])
  speciesSd <- c(NA, speciesVarGLMM[1])

  resultsSpecies <- data.frame(Model = model,
                               Intercept = intercept,
                               Latitude = latitude,
                               `SE of Intercept` = intSd,
                               `SE of latitude` = latSd,
                               `SE of visit effect` = visitSd,
                               `SE of species effect` = speciesSd)

  return(list(species = resultsSpecies,
              group = resultsCounts))
}


### presentation chuff
```{r, echo = FALSE, warning = FALSE, message = FALSE}
#| label: tbl-bbSpe
#| tbl-cap: Coefficient and standaed error estimates from the binomial GLMM fitted to the species occupancy data.
modelFitExploration(1)$species%>%
  kbl(format="latex", booktabs=TRUE) %>%
  kable_styling(full_width = F, font_size = 8)
```


```{r, echo = FALSE, warning = FALSE, message = FALSE}
#| label: tbl-bbGro
#| tbl-cap: Coefficient and standaed error estimates from the poisson and negative binomial GLMM fitted to the species occupancy data.
modelFitExploration(1)$group%>%
  kbl(format="latex", booktabs=TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))
```

```{r, echo = FALSE, warning = FALSE, message = FALSE}
#| label: tbl-hvSpe
#| tbl-cap: Coefficient and standaed error estimates from the binomial GLMM fitted to the species occupancy data.
modelFitExploration(2)$species%>%
  kbl(format="latex", booktabs=TRUE)%>%
  kable_styling(full_width = F, font_size = 8)
```


```{r, echo = FALSE, warning = FALSE, message = FALSE}
#| label: tbl-hvGro
#| tbl-cap: Coefficient and standaed error estimates from the poisson and negative binomial GLMM fitted to the species occupancy data.
modelFitExploration(2)$group%>%
  kbl(format="latex", booktabs=TRUE)%>%
  kable_styling(latex_options = c("striped", "scale_down"))
```

```{r, echo = FALSE, warning = FALSE, message = FALSE}
#| label: tbl-sbSpe
#| tbl-cap: Coefficient and standaed error estimates from the binomial GLMM fitted to the species occupancy data.
modelFitExploration(3)$species%>%
  kbl(format="latex", booktabs=TRUE)%>%
  kable_styling(full_width = F, font_size = 8)
```


```{r, echo = FALSE, warning = FALSE, message = FALSE}
#| label: tbl-sbGro
#| tbl-cap: Coefficient and standaed error estimates from the poisson and negative binomial GLMM fitted to the species occupancy data.
modelFitExploration(3)$group%>%
  kbl(format="latex", booktabs=TRUE)%>%
  kable_styling(latex_options = c("striped", "scale_down"))
```
