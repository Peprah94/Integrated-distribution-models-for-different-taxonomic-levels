setwd("/Data_idm")

require(dplyr)
require(ggcorrplot)
require(purrr)
require(ggpubr)
require(heatmaply)
require(reshape2)
library(ggplot2)
require(egg)
require(readr)
if(!require(devtools)) install.packages("devtools")
if(!require(ggcorrplot)) devtools::install_github("kassambara/ggcorrplot")
if(!require(ggcorrplot)) install.packages("heatmaply")


#load data
load("Results/data_idm.RData")

#load estimates of parameters
load("Results/estimate_data.RData")

###############
# Exploration of data analysed
##############
all_data <- lapply(simulations_all, function(x){
  species_data <- x$mat.species%>%
    apply(.,c(1,2,3,4), function(x){x != 0})%>%
    apply(., c(1,4), function(x) mean(x, na.rm = TRUE))
  
  
  genus_data <-  x$mat.genus%>%
    apply(.,c(1,2,3), function(x){x != 0})%>%
    apply(., c(1,3), function(x) mean(x, na.rm = TRUE))
  
  return(list(species_data, genus_data))
})

all_data1 <- purrr::flatten(all_data)

all_data2 <- do.call("cbind", all_data1)%>%
  as.data.frame()%>%
  cbind(., sites = simulations_all[[1]]$sites$X1km_square)

colnames(all_data2) <- c("bumblebees_PA_17", "bumblebees_PA_18", 
                         "bumblebees_count_17", "bumblebees_count_18",
                         "hoverflies_PA_17", "hoverflies_PA_18", 
                         "hoverflies_count_17", "hoverflies_count_18",
                         "solitarybees_PA_17", "solitarybees_PA_18", 
                         "solitarybees_count_17", "solitarybees_count_18", 
                         "Sites")

rownames(all_data2) <- simulations_all[[1]]$sites$X1km_square
is.nan.data.frame <- function(x){
  return(do.call(cbind, lapply(x, is.nan)))
}
all_data2[is.nan(all_data2)] <- NA

melt_data <- reshape2::melt(all_data2, id.vars = c("Sites"))

ggplot(data = melt_data, mapping = aes(x= as.factor(variable), y= as.factor(Sites), fill = value))+
  geom_tile()

heatmaply(x = all_data2[,-13], 
          xlab = "Data type and year",
          ylab = "Site index",
          k_col = NA, 
          #subplot_margin = 1,
          k_row = NA,
          row_text_angle = 0,
          column_text_angle = 270,
          fontsize_row = 5,
          fontsize_col = 5,
          #hclust_method = "complete",
          plot_method = "plotly",
          dendrogram = "none",
          showticklabels = TRUE)%>% 
  layout(
    font = list(
      family = "Times New Roman",
      font ='italic'))


############
# Plot of estimates
########
#Extract the covariance estimates
extract_covariance <- function(x,y){
  constants <- x[[1]]
  correlation_matrix <- matrix(x[[2]][grep("Cov", rownames(x[[2]])), 1], 
                              nrow=constants$n.species,
                              ncol=constants$n.species,
                              byrow = F)%>%
    cov2cor()
  
  rownames(correlation_matrix) <- y
  colnames(correlation_matrix) <- y
  return(correlation_matrix)

}

# function to extract shannon index for all estimates
extract_shan <- function(x){
  constants <- x[[1]]
  shan_index <- matrix(x[[2]][ grep("shan", rownames(x[[2]])), 1],
                       nrow=constants$n.sites,
                       ncol=2,
                       byrow = T)
  return(shan_index)
}

#function to extract the sd of the site effects
extract_sigma <- function(x){
  sigma <- matrix(x[[2]][ grep("sigm", rownames(x[[2]])), 1],
                       nrow=1,
                       ncol=2,
                       byrow = T)
  return(sigma)
}

# Getting Shannon Index from data
shans <- lapply(estimates, extract_shan)%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(site = simulations_all[[1]]$sites$X1km_square)

 colnames(shans) <- c("bbIDM17", "bbIDM18", 
           "bbSPE17", "bbSPE18",
           "bbGRO17", "bbGRO18",
           "sbIDM17", "sbIDM18", 
           "sbSPE17", "sbSPE18",
           "sbGRO17", "sbGRO18",
           "hvIDM17", "hvIDM18", 
           "hvSPE17", "hvSPE18",
           "hvGRO17", "hvGRO18",
           "sites")

 names_list <- c(rep(list(simulations_all[[1]]$species_names),3),
                   rep(list(simulations_all[[2]]$species_names),3),
                 rep(list(simulations_all[[3]]$species_names),3)
                 )

 #Abbreviate the names of the species
abb_names_list <- lapply(names_list, function(x){
  abbreviate(x, minlength = 8)
})
                   
 correlation_matrices <- pmap(list(estimates, abb_names_list), extract_covariance)
   
 # plotting corr heatmap
   cor_plots1 <-  lapply(correlation_matrices, function(y){
                                      heatmaply_cor(x = y, 
                                     xlab = "Species list",
                                     ylab = "Species list",
                                     k_col = NA, 
                                     k_row = NA,
                                     row_text_angle = 0,
                                     column_text_angle = 90,
                                     fontsize_row = 3.5,
                                     fontsize_col = 3.5,
                                     #hclust_method = "complete",
                                     plot_method = "ggplot",
                                     dendrogram = "none",
                                     showticklabels = TRUE,
                                     margins = c(10, 10))%>% 
     layout( xaxis = list(titlefont = list(size = 10)),
             yaxis = list(titlefont = list(size = 10)) )
 })
 
 subplot(cor_plots1[[1]], cor_plots1[[2]],cor_plots1[[3]],
         cor_plots1[[4]], cor_plots1[[5]],cor_plots1[[6]],
         cor_plots1[[7]], cor_plots1[[8]],cor_plots1[[9]],
         margin = .06, 
         nrows = 3,
         #widths = c(0.2, 0.4, 0.4),
         heights = c(0.2,0.35,0.35),
         shareY = FALSE, 
         shareX = FALSE,
         titleX = TRUE, 
         titleY = TRUE)
 
 #save the file
 png(filename = "species_associations1.png",
     units = "in",
     width = 9.58,
     height = 9.18,
     res = 72)
 
 #extracting the sd of the site effects
 sigmas <- lapply(estimates, extract_sigma)%>%
   do.call("rbind", .)%>%
   as.data.frame()%>%
   mutate(method = rep(c("IDM" , "SPE", "GO"), 3),
          insects = rep(c("Bumblebees", "Solitary bees", "Hoverflies"),each= 3))

colnames(sigmas) <- c("Detection", "Observed", "Methods", "Insects")
write.csv(sigmas,"Results/sigmas.csv", row.names = FALSE)

# Plot variances
plot_sigmas <- melt(sigmas, id.vars = c("Methods", "Insects"))%>%
  ggplot(aes(x= Methods, y = value, group = Insects, col= Insects))+
  geom_point()+
  geom_line(linetype="dashed")+
  theme_article()+
  facet_wrap(~variable, nrow=1)
ggsave(filename = "Figure/site_variance.png",plot_sigmas)


#saving the shannon diversity with the region and site ID
gr_ref <- read_csv("Results/gr_ref.csv")
names(gr_ref)[1] <- "sites"
sites <- data.frame(sites = simulations_all[[1]]$sites$X1km_square)
merged_data_sites <- merge(x=sites, 
                           y=gr_ref, 
                           by="sites", 
                           all.x = TRUE, 
                           all.y = FALSE)%>%
  select(sites, region)%>%
  arrange(sites)

#Fill in the missing region for a location
merged_data_sites[62,2] <- "ENG"
all_data_estimates <- merge(shans, 
                            merged_data_sites, 
                            by = "sites", 
                            all.x = TRUE, 
                            all.y = FALSE)
  
write.csv(all_data_estimates,"Results/shannon_index.csv", row.names = FALSE)