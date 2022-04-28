# Packages needed to run the code
library(ggplot2)
library(pbapply)
library(reshape2)
library(purrr)
library(dplyr)
library(gganimate)
library(ggpubr)
library(pbapply)
library(ggplot2)
library(cowplot)

load("/idm_data/data_idm.RData")
source("/idm_data/fnx_for estimation.R")

#sites
site_id <- simulations_all[[1]]$sites$X1km_square

#######################
#BUMBLEBEES
############################
load("/idm_data/estimate_inter_data.RData")
load("/idm_data/estimate_species_data.RData")
load("/idm_data/estimate_genus_data.RData")
bb_inter_17 <- inter_estimates_data
bb_species_17 <- species_estimates_data
bb_genus_17 <- genus_estimates_data

load("/idm_data1/estimate_inter_data.RData")
load("/idm_data1/estimate_species_data.RData")
load("/idm_data1/estimate_genus_data.RData")
bb_inter_18 <- inter_estimates_data
bb_species_18 <- species_estimates_data
bb_genus_18 <- genus_estimates_data

# HOVERFLIES
load("/idm_data2/estimate_inter_data.RData")
load("/idm_data2/estimate_species_data.RData")
load("/idm_data2/estimate_genus_data.RData")
hv_inter_17 <- inter_estimates_data
hv_species_17 <- species_estimates_data
hv_genus_17 <- genus_estimates_data

load("/idm_data3/estimate_inter_data.RData")
load("/idm_data3/estimate_species_data.RData")
load("/idm_data3/estimate_genus_data.RData")
hv_inter_18 <- inter_estimates_data
hv_species_18 <- species_estimates_data
hv_genus_18 <- genus_estimates_data

# SOLITARY BEES
load("/idm_data4/estimate_inter_data.RData")
load("/idm_data4/estimate_species_data.RData")
load("/idm_data4/estimate_genus_data.RData")
sb_inter_17 <- inter_estimates_data
sb_species_17 <- species_estimates_data
sb_genus_17 <- genus_estimates_data

load("/idm_data5/estimate_inter_data.RData")
load("/idm_data5/estimate_species_data.RData")
load("/idm_data5/estimate_genus_data.RData")
sb_inter_18 <- inter_estimates_data
sb_species_18 <- species_estimates_data
sb_genus_18 <- genus_estimates_data


# estimating_shan <- function(x){
#   p <- matrix(x[[2]][(2*x[[1]]$n.species +1): (2*x[[1]]$n.species + (x[[1]]$n.sites*x[[1]]$n.species)),1],
#               nrow = x[[1]]$n.sites,
#               ncol=x[[1]]$n.species,
#               byrow=FALSE)
#   shan <- shan_index(p)
#   return(shan)
# }

estimating_shan <- function(x){
  p <- tail(x[[2]], x[[1]]$n.sites)
  #shan <- shan_index(p)
  return(p)
}


bb_shan_inter_17 <- estimating_shan(bb_inter_17)
bb_shan_species_17 <- estimating_shan(bb_species_17)
bb_shan_genus_17 <- estimating_shan(bb_genus_17)
bb_shan_inter_18 <- estimating_shan(bb_inter_18)
bb_shan_species_18 <- estimating_shan(bb_species_18)
bb_shan_genus_18 <- estimating_shan(bb_genus_18)

hv_shan_inter_17 <- estimating_shan(hv_inter_17)
hv_shan_species_17 <- estimating_shan(hv_species_17)
hv_shan_genus_17 <- estimating_shan(hv_genus_17)
hv_shan_inter_18 <- estimating_shan(hv_inter_18)
hv_shan_species_18 <- estimating_shan(hv_species_18)
hv_shan_genus_18 <- estimating_shan(hv_genus_18)

sb_shan_inter_17 <- estimating_shan(sb_inter_17)
sb_shan_species_17 <- estimating_shan(sb_species_17)
sb_shan_genus_17 <- estimating_shan(sb_genus_17)
sb_shan_inter_18 <- estimating_shan(sb_inter_18)
sb_shan_species_18 <- estimating_shan(sb_species_18)
sb_shan_genus_18 <- estimating_shan(sb_genus_18)

#Putting all data together

# all_data_17 <- data.frame(bb_shan_inter_17=bb_shan_inter_17,
#                        bb_shan_species_17=bb_shan_species_17,
#                        bb_shan_genus_17=bb_shan_genus_17,
#                        sites = site_id)

bb_all_data_17 <- data.frame(rbind(bb_shan_inter_17,
                                bb_shan_species_17,
                                bb_shan_genus_17),
                          sites = rep(site_id, 3),
                          model= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_17)),
                          year = "2017")

bb_all_data_18 <- data.frame(rbind(bb_shan_inter_18,
                                bb_shan_species_18,
                                bb_shan_genus_18),
                          sites = rep(site_id, 3),
                          model= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_18)),
                          year = "2018")

bb_all_data <- data.frame(rbind(bb_all_data_17,
                                bb_all_data_18),
                          group = "Bumblebees")




#####################
#   Hoverflies
##################
hv_all_data_17 <- data.frame(rbind(hv_shan_inter_17,
                                hv_shan_species_17,
                                hv_shan_genus_17),
                          sites = rep(site_id, 3),
                          model= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(hv_shan_inter_17)),
                          year = "2017")
hv_all_data_18 <- data.frame(rbind(hv_shan_inter_18,
                                hv_shan_species_18,
                                hv_shan_genus_18),
                          sites = rep(site_id, 3),
                          model= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(hv_shan_inter_18)),
                          year = "2018")

hv_all_data <- data.frame(rbind(hv_all_data_17,
                                hv_all_data_18),
                          group = "Hoverflies")





######################
#     Solitary bees
######################
sb_all_data_17 <- data.frame(rbind(sb_shan_inter_17,
                                sb_shan_species_17,
                                sb_shan_genus_17),
                          sites = rep(site_id, 3),
                          model= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(sb_shan_inter_17)),
                          year = "2017")
sb_all_data_18 <- data.frame(rbind(sb_shan_inter_18,
                                sb_shan_species_18,
                                sb_shan_genus_18),
                          sites = rep(site_id, 3),
                          model= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(sb_shan_inter_18)),
                          year = "2018")

sb_all_data <- data.frame(rbind(sb_all_data_17,
                                sb_all_data_18),
                          group = "Solitary bees")



all_data <- rbind(bb_all_data, hv_all_data, sb_all_data)


# old figure

gg1 <- ggplot(all_data, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(model),color=as.factor(model)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(model), fill=as.factor(model)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylab("Shannon Index")+
  xlab("Site ID")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,4)) +
  facet_wrap(~group+year, nrow = 3, ncol = 2)

# new figure

ggplot(data = all_data) +
  geom_boxplot(aes(x = group, y = Mean,
                   fill = model,
                   color = model),
               alpha = 0.2) +
  facet_wrap(~year) +
  ylab("Mean Shannon index across sites") +
  xlab("Genus") +
  theme_bw()+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                    name="Method") +
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                     name="Method")
  

# map of sites (represented at 10Km because too small to see at 1Km)

library(BRCmap)

# par(mar = c(0.1,0.1,1,0.1))
plot_GIS(UK[UK$REGION == "Great Britain",], xlab="", ylab="", show.axis = FALSE, show.grid = FALSE, 
         fill.col = "lightgrey", no.margin = TRUE)

plotUK_gr(all_data$sites, gr_prec = 10000, border = "black")


# map of shannon

coords <- gr_let2num(all_data$sites)

all_data <- cbind(all_data, coords)

UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

ggplot(data = all_data %>%
         filter(group == "Bumblebees")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = Mean), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")




