###############################################################################
# Functions to transform the basic specimen and name change data, including   #
# generating conservation assessments and measuring change in EOO.            #
###############################################################################

library(rCAT)

assess_frame <- function(species_df, cellsize) {
  #' Generate conservation assessments for species in a data frame.
  #' 
  #' Calculates the EOO, and AOO for each of a species in a data frame
  #' and generates conservation assessments based on these, using the
  #' rCAT package.
  #' 
  #' @param species_df Data frame of specimens for one or more species.
  #' 
  #' @param cellsize The cell size, in metres, to use for the AOO
  #' calculations.
  #' 
  #' @return A data frame of species with their EOO, AOO, and
  #' conservation assessments.
  
  species_list <- unique(species_df$species)
  results_df <- data.frame(species=character(), 
                           occurrences=integer(),
                           eoo=character(), 
                           eoo_area=double(), 
                           aoo=character(), 
                           aoo_area=double(),
                           stringsAsFactors = FALSE)
  
  for (species in species_list) {
    species_coordinates <- (species_df[species_df$species==species, 
                                       c("latitude", "longitude")])
    colnames(species_coordinates) <- c("lat", "long")
    species_xy <- simProjWiz(species_coordinates, trueCOGll(species_coordinates))
    
    occurences <- nrow(species_xy)
    eoo_area <- EOOarea(species_xy) / -1000000
    aoo_area <- AOOsimp(species_xy, cellsize) * (cellsize / 1000)^2
    eoo_rating <- EOORating(eoo_area)
    aoo_rating <- AOORating(aoo_area)
    
    results_df[nrow(results_df)+1, ] <- list(species, occurences, eoo_rating,
                                             as.numeric(eoo_area), 
                                             aoo_rating, as.numeric(aoo_area))
  }
  results_df$occurrences <- as.numeric(results_df$occurrences)
  return(results_df)
}

classify_changes <- function(specimen_data, name_changes) {
  #' Classify changes to specimen data
  #' 
  #' @description
  #' This classifies changes to specimen data by first
  #' splitting the specimen data frame into new and old data,
  #' and then getting name pairs for each name change type
  #' from the name changes data frame. The old data is joined
  #' to the new data by specimen id and the changes to each
  #' specimen is classified according to the name change type,
  #' whether the coordinates have changed, or if the specimen
  #' is new to the database.
  #' 
  #' @param specimen_data A data frame of data about all
  #' specimens, including the specimen id, species name, 
  #' coordinates, and which version of the database it is from.
  #' 
  #' @param name_changes A data frame of name changes between
  #' the two databases, including the name in 2007, the name
  #' in 2017, and the type of change
   
  old_data <- 
    specimen_data %>%
    filter(dataset == 2007) %>%
    select(-dataset)
  
  new_data <- 
    specimen_data %>%
    filter(dataset == 2017) %>%
    select(-dataset)
  
  nomenclatural_pairs <- 
    name_changes %>%
    filter(type == "nomenclatural") %>%
    mutate(pairs = str_c(species_2007, species_2017)) %>%
    pull(pairs)
  
  lumped_pairs <- 
    name_changes %>%
    filter(type == "lumped") %>%
    mutate(pairs = str_c(species_2007, species_2017)) %>%
    pull(pairs)
  
  split_pairs <- 
    name_changes %>%
    filter(type == "split") %>%
    mutate(pairs = str_c(species_2007, species_2017)) %>%
    pull(pairs)
  
  correction_pairs <- 
    name_changes %>%
    filter(type == "correction") %>%
    mutate(pairs = str_c(species_2007, species_2017)) %>%
    pull(pairs)
  
  old_data %>%
    full_join(new_data, by = "id", suffix = c("_2007", "_2017")) %>%
    mutate(new_specimen = is.na(species_2007) & !is.na(species_2017),
           new_coordinates = (is.na(latitude_2007) | is.na(longitude_2007)) & !is.na(latitude_2017) & !is.na(longitude_2017) & !is.na(species_2007),
           corrected_coordinates = (latitude_2007 != latitude_2017 | longitude_2007 != longitude_2017) & !is.na(latitude_2007) & !is.na(longitude_2007) & !is.na(latitude_2017) & !is.na(longitude_2017),
           upgraded_to_species = !determined_to_species_2007 & determined_to_species_2017 & !is.na(species_2007),
           downgraded_to_genus = determined_to_species_2007 & !determined_to_species_2017 & !is.na(species_2017),
           lumped = str_c(species_2007, species_2017) %in% lumped_pairs,
           split = str_c(species_2007, species_2017) %in% split_pairs,
           nomenclatural = str_c(species_2007, species_2017) %in% nomenclatural_pairs,
           correction = str_c(species_2007, species_2017) %in% correction_pairs,
           no_change = !(new_specimen | new_coordinates | corrected_coordinates | upgraded_to_species | downgraded_to_genus | lumped | split| nomenclatural | correction)) %>%
    select(-starts_with("determined_to_species")) %>%
    gather(change, has_change, -id, -species_2007, -species_2017, -latitude_2007, -latitude_2017, -longitude_2007, -longitude_2017, -collection_year_2007, -collection_year_2017) %>%
    filter(has_change) %>%
    select(-has_change)
  
}

get_all_assessed_species <- function(records, name_changes, old_or_new = "old") {
  #' Get the species name for which a specimen was used in a conservation assessment
  #' for at a specific time.
  #' 
  #' A conservation assessment for a particular species may be listed under
  #' another species name, depending on what changes have occured to a
  #' particular species name. This finds which species name the specimen was
  #' associated to for a new or old conservation assessment.
  #' 
  #' @param records A data frame of records with an old and new species name for a
  #' list of specimens, as well as the type of name change that has occurred.
  #' 
  #' @param name_changes A data frame of name changes classified by type.
  #' 
  #' @param old_or_new A string, "old" or "new", specifying which name to get the
  #' conservation assessment name for.
  #' 
  #' @return A vector of species names that link to a conservation assessment.
  
  get_assessed_species <- switch(old_or_new,
                                 old = get_assessed_species_old,
                                 new = get_assessed_species_new)
  
  
  species <- c()
  for (row in seq(1, nrow(records))) {
    species <- c(species, get_assessed_species(records[row, ], name_changes))
  }
  
  return(species)
}

get_assessed_species_new <- function(record, name_changes) {
  #' Get the assessed species name for 2007 for a specimen.
  #' 
  #' @param record Data for the specimen to find the new
  #' assessed species name for.
  #' 
  #' @param name_changes A data frame of name changes that
  #' occurred between the two database snapshots, and their type.
  #' 
  #' @return The new assessed species name.
  
  switch(record$change,
         nomenclatural = record$species_2017,
         no_change = record$species_2017,
         new_coordinates = record$species_2017,
         corrected_coordinates = record$species_2017,
         new_specimen = record$species_2017,
         lumped = record$species_2017,
         upgraded_to_species = record$species_2017,
         corrected_in = record$species_2017,
         split = find_new_name(record$species_2007, name_changes),
         downgraded_to_genus = find_new_name(record$species_2007, name_changes),
         corrected_out = find_new_name(record$species_2007, name_changes),
         NA)
}

get_assessed_species_old <- function(record, name_changes) {
  #' Get the assessed species name for 2007 for a specimen.
  #' 
  #' @param record Data for the specimen to find the old
  #' assessed species name for.
  #' 
  #' @param name_changes A data frame of name changes that
  #' occurred between the two database snapshots, and their type.
  #' 
  #' @return The old assessed species name.
  
  switch(record$change,
         nomenclatural = record$species_2007,
         no_change = record$species_2007,
         new_coordinates = record$species_2007,
         corrected_coordinates = record$species_2007,
         new_specimen = find_old_name(record$species_2017, name_changes),
         lumped = find_old_name(record$species_2017, name_changes),
         upgraded_to_species = find_old_name(record$species_2017, name_changes),
         corrected_in = find_old_name(record$species_2017, name_changes),
         split = record$species_2007,
         downgraded_to_genus = record$species_2007,
         corrected_out = record$species_2007,
         NA)
}

find_old_name <- function(new_name, name_changes) {
  #' Find the old name of a species.
  #' 
  #' A specimen may have had a nomenclatural
  #' change, but be listed under another change name
  #' as purely nomenclatural changes do not affect EOO. This
  #' function finds the old name for a species, even
  #' if it has had a nomenclatural change.
  #' 
  #' @param new_name The species name the specimen has in
  #' the 2017 snapshot of the database.
  #' 
  #' @param name_changes A data frame of name changes between
  #' the two database snapshots, classified by type.
  #' 
  #' @return A string, the old species name.
  
  nomenclatural_change <-
    name_changes %>%
    filter(species_2017 == new_name, 
           type == "nomenclatural")
  
  if (nrow(nomenclatural_change) < 1) {
    return(new_name)
  }
  
  return(nomenclatural_change$species_2007)
}

find_new_name <- function(old_name, name_changes) {
  #' Find the new name of a species.
  #' 
  #' A specimen may have had a nomenclatural
  #' change, but be listed under another change name
  #' as purely nomenclatural changes do not affect EOO. This
  #' function finds the new name for a species, even
  #' if it has had a nomenclatural change.
  #' 
  #' @param new_name The species name the specimen has in
  #' the 2007 snapshot of the database.
  #' 
  #' @param name_changes A data frame of name changes between
  #' the two database snapshots, classified by type.
  #' 
  #' @return A string, the new species name.
  
  nomenclatural_change <-
    name_changes %>%
    filter(species_2007 == old_name,
           type == "nomenclatural")
  
  if (nrow(nomenclatural_change) < 1) {
    return(old_name)
  }
  
  return(nomenclatural_change$species_2017)
}


get_specimen_eoo_changes <- function(records) {
  #' Find the EOO change attributable to each specimen.
  #' 
  #' Iterates over each row of a data frame containing
  #' information about changes to specimen data between
  #' the two snapshots of the database, and calculates the
  #' change in EOO for the associated species when just
  #' that specimen change is made to the data for that
  #' species.
  #' 
  #' @param records A data frame where each row is a
  #' specimen, and its data in the 2007 and 2017 database
  #' snapshots.
  #' 
  #' @return A vector containing the EOO change for each
  #' specimen in records.
  
  results <- c()
  for (row in seq(1, nrow(records))) {
    results <- c(results, find_eoo_change(records[row, ], records))
  }
  
  return(results)
}


find_eoo_change <- function(entry, records) {
  #' Find the EOO change attributable one specimen.
  #' 
  #' Returns NA for anything with a change
  #' that a priori has no effect (no change, 
  #' purely a nomenclatural change) and anything
  #' where the conservation status changes to or
  #' from DD. For everything else, calculates the
  #' EOO by making just this specimen change to
  #' the collection of old specimens associated to
  #' that species, and then takes the difference
  #' from the old EOO.
  #' 
  #' @param entry Data for the specimen being examined.
  #' 
  #' @param records A data frame of all specimen information
  #' in the database
  #' 
  #' @return The change in EOO for this specimen.
  
  if (nrow(entry) < 1) {
    return(NA)
  }
  
  if (entry$change == "no_change" | entry$change == "nomenclatural") {
    return(NA)
  }
  
  if (is.na(entry$latitude_2017) | is.na(entry$longitude_2017)) {
    return(NA)
  }

  old_records <-
    entry %>%
    get_old_records(records)
  
  if (nrow(old_records) < 1) {
    return(NA)
  }
  
  updated_records <-
    entry %>%
    rename(lat = latitude_2017, long = longitude_2017) %>%
    select(id, lat, long) %>%
    make_record_change(entry$change, old_records)
  
  
  if(is.null(nrow(updated_records))) {
    print(updated_records)
    print(entry)
  }
  if (nrow(updated_records) < 1) {
    return(NA)
  }
  
  updated_eoo <-
    updated_records %>%
    calculate_eoo()
  
  old_eoo <-
    old_records %>%
    calculate_eoo()
  
  return(updated_eoo - old_eoo)
}

get_old_records <- function(entry, records) {
  #' Get all old records associated to the species
  #' of a particular specimen change.
  #' 
  #' If the specimen change type is in a particular set,
  #' the old records are all records for the species name
  #' from 2007. Otherwise it is all records for the
  #' species name from 2017.
  #' 
  #' @param entry Data for the specimen being examined.
  #' 
  #' @param records A data frame of all specimen information
  #' in the database
  #' 
  #' @return All records associated to the species.
  
  species <- entry$assessed_species_2007
  
  records %>%
    filter(species_2007 == species, !is.na(latitude_2007), !is.na(longitude_2007)) %>%
    rename(lat = latitude_2007, long = longitude_2007) %>%
    select(id, lat, long)
}

make_record_change <- function(changed_record, type, old_records) {
  #' Update a set of specimen records with a particular specimen change.
  #' 
  #' Depending on the specimen change, this either removes a specimen
  #' from a collection of records, adds one in, or updates the coordinates
  #' of a particular record.
  #' 
  #' @param changed_record Data for the specimen record in question.
  #' 
  #' @param type A string, the type of change that has occured to the specimen
  #' 
  #' @param old_records A data frame of all specimen information
  #' for the associated species in the old database snapshot.
  #' 
  #' @return A data frame of specimen information the same as old_records, but with
  #' the change being examined made to it.
  
  updated_records <- switch(type,
                            new_specimen = rbind(old_records, changed_record),
                            new_coordinates = rbind(old_records, changed_record),
                            upgraded_to_species = rbind(old_records, changed_record),
                            lumped = rbind(old_records, changed_record),
                            corrected_in = rbind(old_records, changed_record),
                            corrected_coordinates = mutate(old_records, 
                                                           lat = ifelse(id == changed_record$id, 
                                                                        changed_record$lat, lat),
                                                           long = ifelse(id == changed_record$id, 
                                                                         changed_record$long, long)),
                            downgraded_to_genus = filter(old_records, id != changed_record$id),
                            split = filter(old_records, id != changed_record$id),
                            corrected_out = filter(old_records, id != changed_record$id))
  
  return(updated_records)
}

calculate_eoo <- function(specimens) {
  #' Calculate the eoo for a given set of specimens.
  #' 
  #' @param specimens A data frame of specimen
  #' latitudes and longitudes.
  #' 
  #' @return The EOO of the specimens.
  
  eoo <- 
    specimens %>%
    select(lat, long) %>%
    simProjWiz(., trueCOGll(.)) %>%
    EOOarea() %>%
    "/"(-1000000)
  
  return(eoo)
}
