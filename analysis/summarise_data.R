###############################################################################
# Script to summarise the data in the 

library(here)
source(here("R", "analysis_functions.R"))
library(tidyverse)
library(knitr)
library(writexl)

# load in the data ------------------------------------------------------------

# if running as a stand-alone script, you will need to specify filepaths here

specimen_data <- read_csv(specimen_data_file)
name_changes <- read_csv(name_change_data_file)

# loading the outputs of the transform_data script here

categories <- c("DD", "LC", "NT", "VU", "EN", "CR")
assessments <- 
  here("output", "conservation_assessments.csv") %>%
  read_csv() %>%
  mutate(eoo = factor(eoo, levels=categories, ordered=TRUE),
         dataset = factor(dataset, levels=c("2007", "2017", ordered=TRUE)))

eoo_changes <- read_csv(here("output", "specimen_eoo_changes.csv"))
specimen_changes <- read_csv(here("output", "specimen_changes.csv"))
species_status_changes <- read_csv(here("output", "species_status_changes.csv"))


# conservation status summaries -----------------------------------------------
conservation_status_summary <-
  assessments %>%
  count(eoo, dataset)

status_change_summary <-
  species_status_changes %>%
  count(eoo_2007, eoo_2017) %>%
  complete(eoo_2007, eoo_2017, fill = list(n = 0)) %>%
  rename(preliminary = eoo_2007,
         final = eoo_2017)

# raster cell summaries -------------------------------------------------------
raster_ylim <- c(-40, 35)
raster_xlim <- c(-115, -25)

cell_size_km <- 100
project <- TRUE

old_specimen_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2007) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'specimen_count', 
                    project = project, 
                    xlim=raster_xlim, 
                    ylim=raster_ylim)

new_specimen_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2017) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'specimen_count', 
                    project = project, 
                    xlim=raster_xlim, 
                    ylim=raster_ylim)

old_species_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2007) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'species_count', project = project, xlim=raster_xlim, ylim=raster_ylim)

new_species_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2017) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'species_count', project = project, xlim=raster_xlim, ylim=raster_ylim)

# create overall summary table ------------------------------------------------

spec_count <- 
  specimen_data %>%
  group_by(dataset) %>%
  summarise(n_specimens = n())

species_count <- 
  specimen_data %>%
  filter(determined_to_species) %>%
  group_by(dataset) %>%
  summarise(n_species = n_distinct(species))

with_coordinates <- 
  specimen_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  group_by(dataset) %>%
  summarise(n_coordinates = n())

spec_per_spec <- 
  specimen_data %>%
  filter(determined_to_species) %>%
  group_by(dataset, species) %>%
  summarise(n = n()) %>%
  group_by(dataset) %>%
  summarise(per_species = median(n), mn_per_species = mean(n))

cell_count <- 
  as.data.frame(old_specimen_cells, xy = TRUE) %>%
  mutate(dataset = "2007") %>%
  rbind(as.data.frame(new_specimen_cells, xy = TRUE) %>%
          mutate(dataset = "2017")
  ) %>%
  group_by(dataset) %>%
  summarise(n_cells = sum(!is.na(layer)))

cell_spec <- 
  as.data.frame(old_specimen_cells, xy = TRUE) %>%
  mutate(dataset = "2007") %>%
  rbind(as.data.frame(new_species_cells, xy = TRUE) %>%
          mutate(dataset = "2017")
  ) %>%
  group_by(dataset) %>%
  summarise(av_specimens = median(layer, na.rm = TRUE), mn_specimens = mean(layer, na.rm = TRUE))


cell_species <- 
  as.data.frame(old_species_cells, xy = TRUE) %>%
  mutate(dataset = "2007") %>%
  rbind(as.data.frame(new_species_cells, xy = TRUE) %>%
          mutate(dataset = "2017")
  ) %>%
  group_by(dataset) %>%
  summarise(av_species = median(layer, na.rm = TRUE), mn_species = mean(layer, na.rm = TRUE))

cell_specimens_same <- 
  as.data.frame(old_specimen_cells, xy=TRUE) %>%
  inner_join(as.data.frame(new_specimen_cells, xy=TRUE), 
             by = c("x", "y"),
             suffix = c("_2007", "_2017")) %>%
  filter(!is.na(layer_2007), !is.na(layer_2017)) %>%
  gather(dataset, specimens, -x, -y) %>% 
  mutate(dataset = str_sub(dataset, 7, 10)) %>%
  group_by(dataset) %>%
  summarise(mn_specimens_same = mean(specimens, na.rm=TRUE), av_specimens_same = median(specimens, na.rm=TRUE))

cell_species_same <- 
  as.data.frame(old_species_cells, xy=TRUE) %>%
  inner_join(as.data.frame(new_species_cells, xy=TRUE), 
             by = c("x", "y"),
             suffix = c("_2007", "_2017")) %>%
  filter(!is.na(layer_2007), !is.na(layer_2017)) %>%
  gather(dataset, species, -x, -y) %>% 
  mutate(dataset = str_sub(dataset, 7, 10)) %>%
  group_by(dataset) %>%
  summarise(mn_species_same = mean(species, na.rm=TRUE), av_species_same = median(species, na.rm=TRUE))

only_genus_count <-
  specimen_data %>%
  group_by(dataset) %>%
  summarise(n = sum(!determined_to_species))

statistics_table <- 
  Reduce(function(...) merge(..., all = TRUE, by = c("dataset")), 
       list(spec_count, with_coordinates, species_count, spec_per_spec, 
            cell_count, cell_spec, cell_species, cell_specimens_same, 
            cell_species_same, only_genus_count)) %>%
  select(dataset, n_specimens, n_coordinates, n_species, per_species, n_cells,
         av_specimens, av_species, mn_species, mn_specimens, mn_species_same,
         av_species_same, mn_specimens_same, av_specimens_same, mn_per_species, n) %>%
  rename(Specimens = n_specimens,
         Species = n_species,
         `Specimens per Species (median)` = per_species,
         `Specimens per Species (mean)` = mn_per_species,
         `Cells covered` = n_cells,
         `Specimens per Cell (median)` = av_specimens,
         `Specimens per Cell (mean)` = mn_specimens,
         `Specimens per Cell (median, 2007 cells)` = av_specimens_same,
         `Specimens per Cell (mean, 2007 cells)` = mn_specimens_same,
         `Species per Cell (median)` = av_species,
         `Species per Cell (mean)` = mn_species,
         `Species per Cell (median, 2007 cells)` = av_species_same,
         `Species per Cell (mean, 2007 cells)` = mn_species_same,
         `Specimens with Geolocation` = n_coordinates,
         `Specimens only determined to genus` = n) %>%
  gather(metric, value, -dataset) %>%
  spread(dataset, value)


# summarise the specimen changes by species with assessments in both periods --
change_types <- c("corrected_coordinates", "new_coordinates", "corrected_in",
                  "corrected_out", "downgraded_to_genus", "upgraded_to_species", 
                  "split", "lumped",
                  "new_specimen", "nomenclatural", "no_change")

eoo_change_summaries <-
  eoo_changes %>%
  count(assessed_species_2007, change) %>%
  spread(change, n, fill = 0)

no_change_summaries <-
  specimen_changes %>%
  filter(change %in% c("nomenclatural", "no_change")) %>%
  count(species_2007, species_2017, change) %>%
  spread(change, n, fill = 0)

changes_with_assessments <-
  species_status_changes %>%
  left_join(eoo_change_summaries, by = c("species_2007" = "assessed_species_2007")) %>%
  left_join(no_change_summaries, by = c("species_2007", "species_2017")) %>%
  mutate_at(vars(one_of(change_types)), funs(ifelse(is.na(.), 0, .))) %>%
  select(-eoo_area_2007, -eoo_area_2017) %>%
  gather(change, n, -species_2007, -species_2017, -eoo_2007, -eoo_2017, -occurrences_2007, -occurrences_2017) %>%
  filter(!is.na(eoo_2017))

changes_by_assessment_change <-
  changes_with_assessments %>%
  group_by(eoo_2007, eoo_2017, change) %>%
  summarise(n_specimens = sum(n, na.rm=TRUE),
            n_species = sum(n > 0, na.rm=TRUE),
            total_species = n()) %>%
  group_by(eoo_2007, eoo_2017) %>%
  mutate(total_specimens = sum(n_specimens, na.rm=TRUE))

most_common_changes_by_assessment_change <-
  changes_by_assessment_change %>%
  slice(which.max(n_specimens))

most_common_changes_narrowed <-
  changes_by_assessment_change %>%
  filter(!(change %in% c("nomenclatural", "no_change"))) %>%
  slice(which.max(n_specimens))

# test significance of changes in the number species with IUCN statuses -------

# for all conservation statuses

modelled_status_summary <- 
  assessments %>%
  model_proportion_difference(dataset, eoo)

# for threatened and not threatened

modelled_threatened_summary <-
  assessments %>%
  mutate(threatened = eoo > "NT") %>%
  model_binomial_difference(dataset, threatened, 10000)

# summarise the EOO change by specimen change type ----------------------------
eoo_change_types <- c("corrected_coordinates", "new_coordinates", 
                      "corrected_out", "downgraded_to_genus", "corrected_in", "upgraded_to_species", 
                      "split", "lumped", 
                      "new_specimen")

filtered_changes <- 
  eoo_changes %>%
  filter(!(change %in% c("nomenclatural", "no_change"))) %>%
  mutate(eoo_2017 = factor(eoo_2017, levels=categories, ordered=TRUE),
         eoo_2007 = factor(eoo_2007, levels=categories, ordered=TRUE),
         change = factor(change, levels=eoo_change_types, ordered=TRUE),
         moved = eoo_2017 != eoo_2007) %>%
  filter(!(eoo_2007 %in% "DD"), 
         !(eoo_2017 %in% "DD"), 
         !is.na(change), 
         !is.na(eoo_change),
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species, ]$species |
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species, ]$species,
         !is.na(moved))

mean_changes <-
  filtered_changes %>%
  group_by(change) %>%
  summarise(mean = mean(eoo_change), 
            std_err = sd(eoo_change) / sqrt(n()),
            n_specimen = n(),
            n_species = length(unique(assessed_species_2007))) %>%
  mutate(proportion_specimens = n_specimen / sum(n_specimen),
         proportion_species = n_species / sum(n_species))
  

proportions <-
  filtered_changes %>%
  group_by(change, moved) %>%
  summarise(n_specimen = n(), 
            n_species = length(unique(assessed_species_2007)),
            mean = mean(eoo_change),
            std_err = sd(eoo_change) / sqrt(n())) %>%
  group_by(moved) %>%
  mutate(proportion_specimens = n_specimen / sum(n_specimen), proportion_species = n_species / sum(n_species))

# test significance of changes in the number of each specimen change type -----

# looking for a significant difference between numbers for species that have
# and haven't changed conservation status

modelled_change_proportions <-
  filtered_changes %>%
  model_proportion_difference(moved, change)


# save the data ---------------------------------------------------------------

write_xlsx(list(overall_summary = statistics_table,
                conservation_status_summary = conservation_status_summary %>% spread(dataset, n),
                conservation_status_changes = status_change_summary,
                conservation_status_signif = modelled_status_summary,
                conservation_threat_signif = modelled_threatened_summary,
                most_common_changes = most_common_changes_by_assessment_change,
                most_common_changes_filtered = most_common_changes_narrowed,
                change_counts_specimens = changes_by_assessment_change %>% select(-n_species, -total_species) %>% spread(change, n_specimens, fill = 0),
                change_counts_species = changes_by_assessment_change %>% select(-n_specimens, -total_specimens) %>% spread(change, n_species, fill = 0), 
                all_changes = changes_with_assessments,
                changes_numbers_signif = modelled_change_proportions,
                eoo_change_summary = mean_changes,
                eoo_change_proportions = proportions),
           here("output", "myrt_cons_biol_data_summary.xlsx"))




