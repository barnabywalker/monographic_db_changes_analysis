#############################################################################
#### script to generate data needed for plotting and analysis.       ########
#### i.e. assessments and changes in eoo                            #########
#############################################################################

library(here)

source("R/analysis_functions.R")
source("R/record_change_functions.R")

library(tidyverse)

# load in specimen and name change data ---------------------------------------

# if running as a stand-alone script, you will need to specify filepaths here

specimen_data <- read_csv(specimen_data_file)
name_changes <- read_csv(name_change_data_file)

# Make preliminary conservation assessmnets -----------------------------------
categories <- c("DD", "LC", "NT", "VU", "EN", "CR")

old_assessments <- 
  specimen_data %>%
  filter(dataset == 2007, 
         determined_to_species,
         !is.na(longitude),
         !is.na(latitude)) %>%
  select(-id) %>%
  assess_frame(cellsize = 2000) %>%
  as.tibble() %>%
  select(-aoo, -aoo_area) %>%
  mutate(eoo = ifelse(occurrences < 3, "DD", eoo),
         dataset="2007",
         eoo = factor(eoo, levels=categories, ordered=TRUE))
  
new_assessments <- 
  specimen_data %>%
  filter(dataset == 2017, 
         determined_to_species,
         !is.na(longitude),
         !is.na(latitude)) %>%
  select(-id) %>%
  assess_frame(cellsize = 2000) %>%
  as.tibble() %>%
  select(-aoo, -aoo_area) %>%
  mutate(eoo = ifelse(occurrences < 3, "DD", eoo),
         dataset="2017",
         eoo = factor(eoo, levels=categories, ordered=TRUE))

conservation_assessments <- rbind(old_assessments, new_assessments)

# Classify specimen changes ---------------------------------------------------
specimen_changes <- classify_changes(specimen_data, name_changes)


# need to force one type of specimen change on each specimen
changes <-
  specimen_changes %>%
  mutate(value = TRUE) %>%
  spread(change, value, fill = FALSE) %>%
  mutate(nomenclatural = ifelse((new_coordinates | corrected_coordinates), FALSE, nomenclatural),
         corrected_coordinates = ifelse((correction | upgraded_to_species | downgraded_to_genus | split | lumped), 
                                        FALSE, corrected_coordinates),
         new_coordinates = ifelse((correction | upgraded_to_species | downgraded_to_genus | split | lumped), 
                                  FALSE, new_coordinates))

changes <-
  changes %>%
  mutate(corrected_in = correction,
         corrected_out = correction) %>%
  select(-correction)

changes <-
  changes %>%
  gather(change, has_change, -id, -species_2007, -species_2017, -latitude_2007, -latitude_2017, -longitude_2007, -longitude_2017) %>%
  filter(has_change) %>%
  select(-has_change) %>%
  mutate(assessed_species_2007 = get_all_assessed_species(., name_changes, "old"),
         assessed_species_2017 = get_all_assessed_species(., name_changes, "new"))

# Calculate EOO attributable to each specimen change --------------------------
eoo_changes <-
  changes %>%
  left_join(old_assessments, by = c("assessed_species_2007" = "species")) %>%
  left_join(new_assessments, by = c("assessed_species_2017" = "species"), 
            suffix = c("_2007", "_2017")) %>%
  select(-dataset_2007, -dataset_2017) %>%
  mutate(eoo_change = get_specimen_eoo_changes(.)) %>%
  filter(!(change %in% c("nomenclatural", "no_change")))

# Find all the species that changes EOO ---------------------------------------
species_changes <-
  specimen_changes %>%
  filter(species_2007 == species_2017 | change == "nomenclatural") %>%
  count(species_2007, species_2017) %>%
  left_join(old_assessments, by = c("species_2007" = "species")) %>%
  left_join(new_assessments, by = c("species_2017" = "species"), suffix = c("_2007", "_2017")) %>%
  select(-dataset_2017, -dataset_2007, -n) %>%
  filter(!is.na(eoo_2007))

# Save everything to csv files ------------------------------------------------

conservation_assessments %>%
  write_csv(here("output", "conservation_assessments.csv"))

eoo_changes %>%
  write_csv(here("output", "specimen_eoo_changes.csv"))

specimen_changes %>%
  write_csv(here("output", "specimen_changes.csv"))

species_changes %>%
  write_csv(here("output", "species_status_changes.csv"))