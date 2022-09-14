#load data
 
library(tidyverse)
library(lubridate)
library(BRCmap)


#Reading the data

genus_counts <- read.csv("Data_for_idm/idm_data/FIT_counts.csv")
pan_traps <- read.csv("Data_for_idm/idm_data/pan_traps.csv")
species_lookup <- read.csv("Data_for_idm/idm_data/species_lookup.csv")

#Formating the date into year, month and day
pan_traps$date <- dmy_hms(pan_traps$date)
genus_counts$date <- dmy_hms(genus_counts$date)

###joining the insect group data to the species data
#Rename the column names as taxon aggregated to species
names(pan_traps)[names(pan_traps) == "taxon_aggregated"] <- "species" 

#Merge the species pantrap data to the species_lookup and match them to their taxon_group
update_pantraps <- merge(pan_traps, species_lookup, by= "species")

# Format the group taxon data (FSV counts)
update_genus <-genus_counts %>%
  filter(insect_group != "all_insects_total")

pantraps_surveys <- unique(
  update_pantraps[,2:3]) %>%
  mutate(pan_trap = "survey")

FIT_surveys <- unique(
  update_genus[,1:2]) %>%
  mutate(FIT_count = "survey")

all_surveys <- pantraps_surveys %>%
  full_join(FIT_surveys)

write.csv(all_surveys, "Data_for_idm/idm_data/all_surveys.csv",
          row.names = FALSE)


#######################
# attempt at formatting
#######################


update_pantraps_wide <- update_pantraps %>%
  select(-c("taxon_group", "n_traps", "insect_group")) %>%
  pivot_wider(id_cols = c("X1km_square", "date"), 
              names_from = species, values_from = n_traps_present) %>%
  mutate_at(vars(3:169), ~replace_na(., 0))


update_genus_wide <- update_genus %>%
  pivot_wider(id_cols = c("X1km_square", "date"), 
              names_from = insect_group, values_from = count) %>%
  mutate_at(vars(3:12), ~replace_na(., 0))


all_data_wide <- full_join(update_pantraps_wide, update_genus_wide,
                           by = c("X1km_square", "date"))


bumblebees <- all_data_wide %>%
  select(c(
    "X1km_square",
    "date",
    "bumblebees", 
    which(names(all_data_wide)%in%species_lookup[species_lookup$insect_group=="bumblebees", "species"] == TRUE)))



# get coordinates


coords <- gr_let2num(all_data_wide$X1km_square)

all_sites_coords <- unique(cbind(all_data_wide$X1km_square, coords))

write.csv(all_sites_coords, "Data_for_idm/idm_data/sites_coords.csv",
          row.names = FALSE)
