# script to make figures
library(here)

source(here("R", "plotting_functions.R"))
source(here("R", "analysis_functions.R"))

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(circlize)
library(viridis)
library(ggsignif)

# changeable setting ----------------------------------------------------------
if (!exists("save_figures")) {
  save_figures <- FALSE
}

# load data -------------------------------------------------------------------
specimen_data <- read_csv(specimen_data_file)

categories <- c("DD", "LC", "NT", "VU", "EN", "CR")

assessments <- 
  read_csv(here("output", "conservation_assessments.csv")) %>%
  mutate(eoo = factor(eoo, levels=categories, ordered=TRUE),
         dataset = factor(dataset, levels=c("2007", "2017", ordered=TRUE)))

eoo_changes <- read_csv(here("output", "specimen_eoo_changes.csv"))
specimen_changes <- read_csv(here("output", "specimen_changes.csv"))
species_status_changes <- read_csv(here("output", "species_status_changes.csv"))

# conservation status bar chart -----------------------------------------------
modelled_status_summary <- 
  assessments %>%
  model_proportion_difference(dataset, eoo)

assessment_counts <-
  assessments %>%
  count(eoo, dataset)

annotations <- 
  modelled_status_summary %>%
  filter(significant) %>%
  left_join(assessment_counts, by = "eoo") %>%
  group_by(eoo) %>%
  summarise(y = max(n) + 15) %>%
  mutate(significant = TRUE,
         x = as.numeric(eoo) - 0.3,
         xend = as.numeric(eoo) + 0.3)

status_bar_chart <- 
  assessment_counts %>%
  ggbarplot(x = "eoo", y = "n", 
            fill = "dataset", color = "white", 
            position = position_dodge(width=0.75), label = "n", 
            ggtheme = theme_pubclean(), lab.size = 3.375,
            xlab = "Assessment category", ylab = "Number of species",
            lab.vjust = 0.5, lab.hjust = -0.1) +
  geom_segment(data = annotations, aes(x = x, xend = xend, y = y, yend = y)) +
  geom_text(data = annotations, aes(x = (x + xend) / 2, y = y + 4), label = "*", vjust=0.7) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1", name = "")

threatened_bar_chart <-
  assessments %>%
  mutate(threatened = ifelse(eoo > "NT", "threatened", "not threatened")) %>%
  count(threatened, dataset) %>%
  ggbarplot(x = "threatened", y = "n",
            fill = "dataset", color = "white",
            position = position_dodge((width=0.75)), label="n",
            ggtheme = theme_pubclean(), lab.size=3.375, lab.vjust = 0.5, lab.hjust = -0.1) +
  scale_fill_brewer(palette = "Set1", name="") +
  coord_flip() +
  labs(x="Number of species", y="")


if (save_figures) {
  ggsave(here("figures", "status_bar_chart.svg"), height = 4.2, width = 8, units="in", plot = status_bar_chart)
  ggsave(here("figures", "threatened_bar_chart.svg"), height = 4.2, width = 6.5, units = "in", plot = threatened_bar_chart)
} else {
  status_bar_chart
  threatened_bar_chart
}

# chord diagrams of conservation assessments that change ----------------------
status_changes <- 
  species_status_changes %>%
  count(eoo_2007, eoo_2017) %>%
  complete(eoo_2007, eoo_2017, fill = list(n = 0)) %>%
  rename(preliminary = eoo_2007, 
         final = eoo_2017)

if (save_figures) {
  svg(here("figures", "status_change_preliminary.svg"))
}

status_changes %>%
  plot_status_chord(label_outside = "NT")

if (save_figures) {
  dev.off()
}


# make raster plots -----------------------------------------------------------

raster_ylim <- c(-40, 35)
raster_xlim <- c(-115, -25)

map_ylim <- c(-35, 30)
map_xlim <- c(-108, -28)

cell_size_km <- 100
project <- TRUE

## specimen raster maps --------------------------------------------------------

breaks <- c(1, 11, 31, 46, 101, 266)
labels <- c("1 - 10", "11 - 30", "31 - 45", "46 - 100", "101 - 265")

old_specimen_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2007) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'specimen_count', 
                    project = project, 
                    xlim=raster_xlim, 
                    ylim=raster_ylim)

old_specimen_raster <- 
  old_specimen_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

new_specimen_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2017) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'specimen_count', 
                    project = project, 
                    xlim=raster_xlim, 
                    ylim=raster_ylim)

new_specimen_raster <- 
  new_specimen_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

levels <- levels(old_specimen_raster)[[1]][, "ID"]

old_specimen_map <- 
  old_specimen_raster %>%
  ggplot_raster_map(ylim = map_ylim, 
                    xlim = map_xlim, 
                    categorical = TRUE, 
                    levels = levels, 
                    labels = labels,
                    north.symbol = 10, 
                    north.scale = 0.075,
                    scale.size = 500,
                    scale.text.size = 2.2) +
  scale_fill_viridis(discrete = TRUE, name = "") + 
  theme(legend.position = c(0.13, 0.28))

new_specimen_map <- 
  new_specimen_raster %>%
  ggplot_raster_map(ylim = map_ylim, 
                    xlim = map_xlim, 
                    categorical = TRUE, 
                    levels = levels, 
                    labels = labels,
                    north.symbol = 10, 
                    north.scale = 0.075,
                    scale.size = 500,
                    scale.text.size = 2.2) +
  scale_fill_viridis(discrete = TRUE, name = "") + 
  theme(legend.position = c(0.13, 0.28))

diff_specimen_cells <- subtract_rasters(old_specimen_cells, new_specimen_cells)

breaks <- c(-7, 0, 1, 11, 51, 116)
labels <- c("< 0", "0", "1 - 10", "11 - 50", "51 - 115")

diff_specimen_raster <-
  diff_specimen_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

levels <- levels(diff_specimen_raster)[[1]][, "ID"]

colours <- c("#4bb062", '#dddcdc', '#c0dceb', '#68abd0', '#2870b1')

diff_specimen_map <-
  diff_specimen_raster %>%
  ggplot_raster_map(ylim = map_ylim, 
                    xlim = map_xlim, 
                    categorical = TRUE, 
                    levels = levels, 
                    labels = labels,
                    north.symbol = 10, 
                    north.scale = 0.075,
                    scale.size = 500,
                    scale.text.size = 2.2) +
  scale_fill_manual(values = colours, labels = labels, name = "") +
  theme(legend.position = c(0.13, 0.28))

if (save_figures) {
  ggsave(here("figures", "specimen_map_2007.svg"), plot = old_specimen_map)
  ggsave(here("figures", "specimen_map_2017.svg"), plot = new_specimen_map)
  ggsave(here("figures", "specimen_map_difference.svg"), plot = diff_specimen_map)
} else {
  old_specimen_map
  new_specimen_map
  diff_specimen_map
}

## species raster maps ---------------------------------------------------------------

breaks <- c(1, 11, 20, 40, 70)
labels <- c("1 - 10", "11 - 20", "21 - 40", "41 - 70")

old_species_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2007) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
               method = 'species_count', project = project, xlim=raster_xlim, ylim=raster_ylim)

old_species_raster <- 
  old_species_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify


new_species_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2017) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
               method = 'species_count', project = project, xlim=raster_xlim, ylim=raster_ylim)

new_species_raster <- 
  new_species_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

levels <- levels(old_species_raster)[[1]][, "ID"]

old_species_map <- 
  old_species_raster %>%
  ggplot_raster_map(ylim = map_ylim, 
                    xlim = map_xlim, 
                    categorical = TRUE, 
                    levels = levels, 
                    labels = labels,
                    north.symbol = 10, 
                    north.scale = 0.075,
                    scale.size = 500,
                    scale.text.size = 2.2) +
  scale_fill_viridis(discrete = TRUE, name = "") + 
  theme(legend.position = c(0.13, 0.28))

new_species_map <- 
  new_species_raster %>%
  ggplot_raster_map(ylim = map_ylim, 
                    xlim = map_xlim, 
                    categorical = TRUE, 
                    levels = levels, 
                    labels = labels,
                    north.symbol = 10, 
                    north.scale = 0.075,
                    scale.size = 500,
                    scale.text.size = 2.2) +
  scale_fill_viridis(discrete = TRUE, name = "") + 
  theme(legend.position = c(0.13, 0.28))

diff_species_cells <- subtract_rasters(old_species_cells, new_species_cells)

breaks <- c(-6, 0, 1, 2, 11, 25)
labels <- c("< 0", "0", "1", "2 - 10", "11 - 24")

breaks <- c(-6, 0, 1, 2, 6, 11, 25)
labels <- c("< 0", "0", "1", "2 - 5", "6 - 10", "11 - 24")

diff_species_raster <-
  diff_species_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

levels <- levels(diff_species_raster)[[1]][, "ID"]

colours <- c("#4bb062", '#dddcdc', '#d1e5f0', '#90c4dd', '#4393c3', '#2065ab')

diff_species_map <-
  diff_species_raster %>%
  ggplot_raster_map(ylim = map_ylim, 
                    xlim = map_xlim, 
                    categorical = TRUE, 
                    levels = levels, 
                    labels = labels,
                    north.symbol = 10, 
                    north.scale = 0.075,
                    scale.size = 500,
                    scale.text.size = 2.2) +
  scale_fill_manual(values = colours, labels = labels, name = "") +
  theme(legend.position = c(0.13, 0.28))

if (save_figures) {
  ggsave(here("figures", "species_map_2007.svg"), old_species_map)
  ggsave(here("figures", "species_map_2017.svg"), new_species_map)
  ggsave(here("figures", "species_map_difference_red.svg"), plot = diff_species_map)
} else {
  old_species_map
  new_species_map
  diff_species_map
}

# plot counts of specimen changes ---------------------------------------------
specimen_types <- c("corrected_coordinates", "new_coordinates", 
                    "correction", "downgraded_to_genus", "upgraded_to_species", 
                    "split", "lumped", 
                    "new_specimen", "nomenclatural", "no_change")

type_colours <- c('#898989', '#d4d4d4', '#1f78b4', '#6a3d9a', '#cab2d6', '#ff7f00', '#fdbf6f', '#e31a1c', '#33a02c', '#b2df8a')

labels <- str_to_sentence(str_replace_all(rev(specimen_types), "_", " "))

all_specimen_changes_pie <-
  specimen_changes %>%
  filter(!is.na(change)) %>%
  plot_changes_pie(order = rev(specimen_types), labels=labels, colours=type_colours) +
  labs(title = "All Specimen Changes")+
  theme(axis.text.x = element_blank())
  

eoo_change_types <- c("corrected_coordinates", "new_coordinates", 
                      "corrected_out", "corrected_in", "downgraded_to_genus", "upgraded_to_species", 
                      "split", "lumped", 
                      "new_specimen")

eoo_type_colours <- c('#1f78b4', '#6a3d9a', '#cab2d6', '#ff7f00', '#fdbf6f', '#e31a1c', '#fb9a99', '#33a02c', '#b2df8a')

labels <- str_to_sentence(str_replace_all(rev(eoo_change_types), "_", " ")) 

all_eoo_change_specimens <-
  eoo_changes %>%
  mutate(moved = eoo_2007 != eoo_2017) %>%
  filter(!(eoo_2007 %in% "DD"), 
         !(eoo_2017 %in% "DD"), 
         !is.na(change), 
         !is.na(eoo_change),
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species,]$species | 
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species,]$species,
         !is.na(moved)) %>%
  plot_changes_pie(order = rev(eoo_change_types), labels=labels, colours=eoo_type_colours) +
  labs(title = "Specimen Changes Used For EOO Change Calculation") +
  theme(axis.text.x = element_blank())

changed_eoo_specimens <-
  eoo_changes %>%
  mutate(moved = eoo_2007 != eoo_2017) %>%
  filter(!(eoo_2007 %in% "DD"), 
         !(eoo_2017 %in% "DD"), 
         !is.na(change), 
         !is.na(eoo_change),
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species,]$species | 
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species,]$species,
         !is.na(moved)) %>%
  filter(moved == TRUE) %>%
  plot_changes_pie(order = rev(eoo_change_types), labels=labels, colours = eoo_type_colours) +
  labs(title = "Specimen Changes For Species that Changed Conservation Status") +
  theme(axis.text.x = element_blank())


if (save_figures) {
  ggsave(here("figures", "pie_all_specimen_changes.svg"), all_specimen_changes_pie)
  ggsave(here("figures", "pie_eoo_calculation_specimen_changes.svg"), all_eoo_change_specimens)
  ggsave(here("figures", "pie_eoo_change_specimen_changes.svg"), changed_eoo_specimens)
} else {
  all_specimen_changes_pie
  all_eoo_change_specimens
  changed_eoo_specimens
}

# plot changes in eoo attributable to specimen changes ------------------------

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
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species,]$species | 
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species,]$species,
         !is.na(moved))

mean_changes <-
  filtered_changes %>%
  group_by(change) %>%
  summarise(mean = mean(eoo_change), 
            std_err = sd(eoo_change) / sqrt(n())) %>%
  mutate(too_few = change %in% c("downgraded_to_genus", "split"))

mean_change_per_type <- 
  ggplot(data = mean_changes, mapping = aes(x = change, y = mean)) +
  geom_hline(yintercept = 0, linetype = 4, colour = "grey50") +
  geom_errorbar(mapping = aes(ymin = mean - std_err*1.57, 
                              ymax = mean + std_err*1.57,
                              colour = too_few), 
                size = 0.3, width = 0.2, 
                position = position_dodge(width = 0.75)) +
  geom_point(mapping = aes(colour = too_few),
             position = position_dodge(width = 0.75)) +
  coord_flip() +
  theme_pubclean() + 
  scale_x_discrete(labels = str_replace(eoo_change_types, "_", "\n")) +
  scale_colour_manual(values = c("black", "grey50")) +
  guides(colour = FALSE) +
  labs(x = "", y = expression("EOO change /" ~ km^2))

mean_change_per_type_split <-
  filtered_changes %>%
  group_by(change, moved) %>%
  summarise(mean = mean(eoo_change), 
            std_err = sd(eoo_change) / sqrt(n())) %>%
  mutate(moved = ifelse(moved, "Status change", "No status change")) %>%
  ungroup() %>%
  rbind(mean_changes %>% select(-too_few) %>% mutate(moved = "Overall")) %>%
  mutate(moved = factor(moved, levels=c("Overall", "No status change", "Status change"), ordered=TRUE)) %>%
  ggplot(mapping = aes(x = change, y = mean, colour = change)) +
  geom_hline(yintercept = 0, linetype = 4, colour = "grey50") +
  geom_errorbar(mapping = aes(ymin = mean - std_err*1.57, 
                              ymax = mean + std_err*1.57), 
                size = 0.3, width = 0.2, 
                position = position_dodge(width = 0.75)) +
  geom_point(position = position_dodge(width = 0.75)) +
  coord_flip() +
  facet_grid(moved ~ .) +
  theme_pubclean() + 
  scale_x_discrete(labels = str_to_sentence(str_replace(str_replace(eoo_change_types, "_", "\n"), "_", " "))) +
  scale_color_manual(labels = str_replace_all(eoo_change_types, "_", " "), values = rev(eoo_type_colours)) +
  guides(colour = FALSE) +
  labs(x = "", y = expression("EOO change /" ~ km^2)) +
  theme(axis.text.y = element_text(color=rev(eoo_type_colours)), panel.border = element_rect(fill=NA, colour="grey80"))

proportions <-
  filtered_changes %>%
  count(change, moved) %>%
  group_by(moved) %>%
  mutate(proportion = n / sum(n))

modelled_change_proportions <-
  filtered_changes %>%
  model_proportion_difference(moved, change)

annotations <- 
  modelled_change_proportions %>%
  filter(significant) %>%
  left_join(proportions, by = "change") %>%
  group_by(change) %>%
  summarise(y = max(proportion) + 0.01) %>%
  mutate(significant = TRUE,
         x = as.numeric(change) - 0.3,
         xend = as.numeric(change) + 0.3)

proportions_by_type <- 
  proportions %>%
  ggbarplot(x = "change",
            y = "proportion",
            fill = "moved",
            color = "white",
            position = position_dodge(width = 0.75),
            ggtheme = theme_pubclean(),
            xlab = "",
            ylab = "Proportion of\nspecimen changes") +
  geom_segment(data = annotations, aes(x = x, xend = xend, y = y, yend = y)) +
  geom_text(data = annotations, aes(x = (x + xend) / 2, y = y + 0.01), label = "*") +
  coord_flip() +
  scale_fill_manual(values = c("#dddcdc", "#e41a1c"), name = "", labels = c("No change in\nspecies status", "Change in\nspecies status")) +
  scale_x_discrete(labels = str_to_sentence(str_replace(str_replace(eoo_change_types, "_", "\n"), "_", " ")))

combined_eoo_change_figs <- ggarrange(proportions_by_type, 
                                      mean_change_per_type + rremove("y.text"), 
                                      common.legend = T, 
                                      align = "h", 
                                      labels = c("(a)", "(b)"))

if (save_figures) {
  ggsave(here("figures", "eoo_changes_mean.svg"), mean_change_per_type)
  ggsave(here("figures", "eoo_changes_mean_split.svg"), mean_change_per_type_split, height=9)
  ggsave(here("figures", "eoo_changes_proportions.svg"), proportions_by_type)
  ggsave(here("figures", "eoo_changes_combined_plots.svg"), combined_eoo_change_figs, width=9)
} else {
  mean_change_per_type
  mean_change_per_type_split
  proportions_by_type
  combined_eoo_change_figs
}
