# Plot results from PoMS analysis

source("idmFunctions/functionForEstimation.R")
load("CaseStudy/data_idm.RData")

#Extracting values from formatted data for Table 1 in main paper
lapply(simulations_all, function(x){
  sumIsNA <- function(z){
    sum(!is.na(c(z)))/nrow(z)
  }
  missSitesSpecies <- sum(x$mat.species > 0 , na.rm = TRUE)/length(!is.na(x$mat.species))
  missSitesSCounts <- sum(x$mat.genus, na.rm = TRUE)/length(!is.na(x$mat.genus))
  missSitesSpeciesSD <- sd(x$mat.species > 0, na.rm = TRUE)
  missSitesCountsSD <- sd(x$mat.genus, na.rm = TRUE)
  ret <- c(missSitesSpecies,  missSitesSCounts,
           missSitesSpeciesSD, missSitesCountsSD)
  return(ret)
})


# Define parameters
models <- c("IDMSH" , "SOM", "IDMCO", "GCMSH", "GCMCO")
insectGrps <- c("Bumblebees", "Hoverflies", "Solitary bees")
nGroups = length(insectGrps)
nModels = length(models)

#########################
# Plot estimates from MCMC
##########################

#load data from Figures folder
parsEstimates <- read_csv("Figures/sigmaEstimates.csv")
shannonEstimates <- read_csv("Figures/shannonEstimates.csv")
rHatEstimates <- read_csv("Figures/rHatEstimates.csv")


#shannon Index

  averageShannon <- shannonEstimates%>%
    group_by(insects, model) %>%
    summarise(average = mean(mean),
              sd = sd(mean))%>%
    ungroup()%>%
    ggplot( aes(x=insects, y=average, col=as.factor(model))) +
    geom_point(position=position_dodge(0.5))+
    geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                  position=position_dodge(0.5))+
    #facet_grid(~model)+
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    xlab("Insect group")+
    ylab(expression(paste("Posterior mean" %+-% "SE")))+
    theme_classic()+
    coord_trans()+
    scale_color_discrete(labels=c('IDMCO', 'IDMSH', 'GCMCO', 'GCMSH','SOM'))+
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8),
          strip.placement = 'outside',
          strip.background = element_rect(colour = "black", fill = "white"))



  #save plot
  ggsave("Figures/averageShannonPlot.svg",
         plot = averageShannon,
         dpi = 320,
         height = 8,
         width = 10,
         units = "cm")

  
  # load data for MSE of Shannon index from simulation study
  # We intend to put the MSE together with the average Shannon Index plot
  shanEstimates0 <- read_csv("simEstimates/shanSimulationsEstimates0.csv")
  shanEstimates1 <- read_csv("simEstimates/shanSimulationsEstimates1.csv")
  shanEstimates2 <- read_csv("simEstimates/shanSimulationsEstimates2.csv")

  shanEstimates <- rbind(shanEstimates0,
                         shanEstimates1,
                         shanEstimates2)

  shanEst <- shanEstimates%>%
    ggplot()+
    geom_boxplot(aes(x = model, y = V1, fill = truth))+
    # facet_grid(~truth)+
    ylab("MSE(H')")+
    theme_classic()+
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 8))


  gg1 <- ggpubr::ggarrange(averageShannon, shanEst,
                           labels = c("A)", "B)"),
                           nrow = 2, ncol = 1,legend = "top")

  ggsave("Figures/simHills.svg",
         plot = gg1,
         dpi = 720,
         height = 14,
         width = 10,
         units = "cm")

  
  # Plot the Shannon Indices on the map
  # map of sites (represented at 10Km because too small to see at 1Km)
  library(BRCmap)
  coords <- gr_let2num(shannonEstimates$coords)
  shannonEstimates <- cbind(coords, shannonEstimates)
  UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
  UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

  #plot Shannon map for each group
  shannonMap <- lapply(as.list(c("Bumblebees", "Hoverflies" , "Solitary bees")), function(x){
    shannonEstimates%>%
      filter(insects == x)%>%
       mutate(model = ifelse(model == "IGCO", "GCMCO",
                             ifelse(model == "IGSH", "GCMSH", model)))%>%
   ggplot() +
      geom_polygon(data = UK_countries_lowres_df,
                   aes(x = long, y = lat, group = group),
                   fill = "gray93", colour = "black") +
      geom_point(aes(x = EASTING, y = NORTHING,
                     fill = mean), shape = 22, size = 1,
                 color = "gray") +
      facet_grid(~model) +
      cowplot::theme_map(font_size = 10) +
      theme(legend.position = "right",
            #legend.title = element_blank(),
            strip.placement = 'outside',
            strip.background = element_rect(colour = "black", fill = "white"),
            text = element_text(size = 10))+
      #coord_equal() +
      scale_fill_viridis_c(name = "H'")
  })
