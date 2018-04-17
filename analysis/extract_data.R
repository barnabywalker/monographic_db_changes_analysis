################################################################################
# Script to extract data from data base and save in right format for analysis. #
################################################################################

library(here)

source("R/data_processing_functions.R")

library(tidyverse)
library(uuid)

# dump a copy of the database -------------------------------------------------

dump_data(table = "Central_Myrtaceae2007")
dump_data(table = "Central_Myrtaceae2017_updatedCC20171002")


# load specimen data and clean ------------------------------------------------

# some post-hoc recoding of names specific to each dataset is needed

old_data <- 
  load_data(table = "Central_Myrtaceae2007", from_database = TRUE) %>%
  mutate(Checklist = str_trim(Checklist)) %>%
  clean_data() %>%
  mutate(Checklist = recode(Checklist,
                            `Myrcia fenzliana` = "Gomidesia lindeniana",
                            `Myrcia montana` = "Gomidesia montana sp. nov.",
                            `Gomidesia montana` = "Gomidesia montana sp. nov."))

new_data <- 
  load_data(table = "Central_Myrtaceae2017_updatedCC20171002", 
            from_database = TRUE) %>%
  mutate(Checklist = str_trim(Checklist)) %>%
  clean_data() %>%
  mutate(Checklist = recode(Checklist, `Myrcia pubipetala` = "Myrcia isaiana"))

# a subset of some old names need correcting, but not in the database

cordiifolia_specimens <-
  old_data %>%
  inner_join(new_data, by="Specimen_ID", suffix=c("_2007", "_2017")) %>%
  filter(Checklist_2007 == "Gomidesia cordiifolia", 
         Checklist_2017 == "Myrcia hebepetala") %>%
  pull(Specimen_ID)

grandifolia_specimens <-
  old_data %>%
  inner_join(new_data, by="Specimen_ID", suffix=c("_2007", "_2017")) %>%
  filter(Checklist_2007 == "Gomidesia grandifolia", 
         Checklist_2017 == "Myrcia amplexicaulis") %>%
  pull(Specimen_ID)

old_data <-
  old_data %>%
  mutate(Checklist = case_when(Specimen_ID %in% cordiifolia_specimens ~ "Gomidesia affinis",
                               Specimen_ID %in% grandifolia_specimens ~ "Gomidesia crocea",
                               TRUE ~ Checklist))

# rename and select only useful columns ---------------------------------------

old_data <-
  old_data %>%
  select(Specimen_ID, Checklist, LatDD, LongDD) %>%
  rename(id = Specimen_ID,
         species = Checklist,
         latitude = LatDD,
         longitude = LongDD)

new_data <-
  new_data %>%
  select(Specimen_ID, Checklist, LatDD, LongDD) %>%
  rename(id = Specimen_ID,
         species = Checklist,
         latitude = LatDD,
         longitude = LongDD)

specimen_data <- 
  rbind(old_data, new_data) %>%
  mutate(dataset = c(rep("2007", nrow(old_data)), rep("2017", nrow(new_data))),
         determined_to_species = str_detect(species, "[ ]+"))


# load name and clean changes -------------------------------------------------

name_changes <- 
  read_csv(here("data", "20180403_classified_change_types.csv")) %>%
  mutate(type = recode(type,
                       `Purely nomenclatural` = "nomenclatural",
                       `Sunk to older name` = "lumped",
                       `Sunk in` = "lumped",
                       `Segregated from complex` = "split",
                       `Correction` = "correction",
                       `working name published` = "nomenclatural",
                       `Working name sunk into existing species` = "lumped")) %>%
  filter(!(type %in% c("aff", "cf.", "ENL", "REMOVE _FIXED IN DB")),
         species_2007 %in% specimen_data$species,
         species_2017 %in% specimen_data$species)

# save data to files for analysis ---------------------------------------------

specimen_data %>%
  write_csv(here("output", "specimen_data.csv"))

name_changes %>%
  write_csv(here("output", "name_change_data.csv"))
# recode data for publication -------------------------------------------------

id_values <- unique(specimen_data$id)
species_values <- unique(specimen_data$species)

new_ids <- sample.int(length(id_values), length(id_values), replace=FALSE)
new_species <- replicate(length(species_values), UUIDgenerate())

names(new_ids) <- id_values
names(new_species) <- species_values

specimen_data$id <- new_ids[as.character(specimen_data$id)]
specimen_data$species <- new_species[specimen_data$species]

name_changes$species_2007 <- new_species[name_changes$species_2007]
name_changes$species_2017 <- new_species[name_changes$species_2017]

# save recoded data and name mappings -----------------------------------------
id_map <-tibble(id=names(new_ids),
                id_recoded=new_ids)

species_map <- tibble(species=names(new_species),
                      species_recoded=new_species)

specimen_data %>%
  write_csv(here("output", "specimen_data_for_publication.csv"))

name_changes %>%
  write_csv(here("output", "name_change_data_for_publication.csv"))

id_map %>%
  write_csv(here("output", "specimen_id_map.csv"))

species_map %>%
  write_csv(here("output", "specimen_species_name_map.csv"))

