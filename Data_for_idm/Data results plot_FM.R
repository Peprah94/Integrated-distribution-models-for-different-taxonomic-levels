# Packages needed to run the code
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

shannon <- read.csv("./Data_idm/Results/shannon_index.csv")

names(shannon) <- c("bb_IDM_17", "bb_IDM_18", "bb_SPE_17", "bb_SPE_18",
                    "bb_GRO_17", "bb_GRO_18", "sb_IDM_17", "sb_IDM_18",
                    "sb_SPE_17", "sb_SPE_18", "sb_GRO_17", "sb_GRO_18",
                    "hv_IDM_17", "hv_IDM_18", "hv_SPE_17", "hv_SPE_18",
                    "hv_GRO_17", "hv_GRO_18", "sites")
#sites
site_id <- shannon$sites

shannon_long <- pivot_longer(shannon, cols = (-sites),
                             names_to = c("group", "model", "year"),
                             names_sep = "_")


# new figure

Shannon_boxplot <- ggplot(data = shannon_long) +
  geom_boxplot(aes(x = group, y = value,
                   fill = model,
                   color = model),
               alpha = 0.2) +
  facet_wrap(~year) +
  ylab("Mean Shannon index across sites") +
  xlab("Genus") +
  theme_bw()+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                    name="Method", labels = c("Genus only", "IDM", "Species only")) +
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                     name="Method", labels = c("Genus only", "IDM", "Species only")) +
  scale_x_discrete(labels= c("Bumblebees", "Hoverflies", "Solitary bees"))
  


png("Data_idm/Results/shannon_boxplot.png", 
    height = 8, width = 16, units = "cm", res = 320)

Shannon_boxplot

dev.off()


# map of sites (represented at 10Km because too small to see at 1Km)
# devtools::install_github("https://github.com/colinharrower/BRCmap")
library(BRCmap)

# par(mar = c(0.1,0.1,1,0.1))
plot_GIS(UK[UK$REGION == "Great Britain",], xlab="", ylab="", show.axis = FALSE, show.grid = FALSE, 
         fill.col = "lightgrey", no.margin = TRUE)

plotUK_gr(shannon_long$sites, gr_prec = 10000, border = "black")


# map of shannon

coords <- gr_let2num(shannon_long$sites)

all_data <- cbind(shannon_long, coords)

UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

all_data <- all_data %>%
  mutate(model = case_when(model == "GRO" ~ "Genus only",
                           model == "SPE" ~ "Species only",
                           model == "IDM" ~ "IDM"),
         year = case_when(year == 17 ~ 2017,
                          year == 18 ~ 2018))
bb_shannon_map <- ggplot(data = all_data %>%
         filter(group == "bb")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")

png("Data_idm/Results/bb_shannon_map.png", 
    height = 16, width = 14, units = "cm", res = 320)

bb_shannon_map

dev.off()



sb_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "sb")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")

png("Data_idm/Results/sb_shannon_map.png", 
    height = 16, width = 14, units = "cm", res = 320)

sb_shannon_map

dev.off()


hv_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "hv")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")

png("Data_idm/Results/hv_shannon_map.png", 
    height = 16, width = 14, units = "cm", res = 320)

hv_shannon_map

dev.off()
