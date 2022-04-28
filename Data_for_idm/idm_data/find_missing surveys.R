#load data
 
library(tidyverse)
library(lubridate)


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
