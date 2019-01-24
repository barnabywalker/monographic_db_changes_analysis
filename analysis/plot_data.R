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

# overall constants -----------------------------------------------------------
eoo_change_types <- c("corrected_coordinates", "new_coordinates", 
                      "corrected_out", "corrected_in", "downgraded_to_genus", "upgraded_to_species", 
                      "split", "lumped", 
                      "new_specimen")

# figure 1 --------------------------------------------------------------------

# This is the figure mapping the number and changes of specimens and species

# overall settings
raster_ylim <- c(-40, 35)
raster_xlim <- c(-115, -25)

map_ylim <- c(-35, 30)
map_xlim <- c(-108, -28)

cell_size_km <- 100
project <- TRUE

#### specimen rasters

# settings for binning the values in the raster cells
breaks <- c(1, 11, 31, 46, 101, 266)
labels <- c("1 - 10", "11 - 30", "31 - 45", "46 - 100", "101 - 265")

# converting 2007 specimen counts to rasters
old_specimen_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2007) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'specimen_count', 
                    project = project, 
                    xlim=raster_xlim, 
                    ylim=raster_ylim)

# binning the 2007 specimen counts
old_specimen_raster <- 
  old_specimen_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

# converting 2017 specimen counts to rasters
new_specimen_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2017) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'specimen_count', 
                    project = project, 
                    xlim=raster_xlim, 
                    ylim=raster_ylim)

# binning the 2017 specimen counts
new_specimen_raster <- 
  new_specimen_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

# getting the levels from the rasters for plotting
levels <- levels(old_specimen_raster)[[1]][, "ID"]

# 2007 specimen map (fig 1a)
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

# 2017 specimen map (fig 1b)
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

# subtracting the two rasters above for the difference map
diff_specimen_cells <- subtract_rasters(old_specimen_cells, new_specimen_cells)

# need different settings to bin these
breaks <- c(-7, 0, 1, 11, 51, 116)
labels <- c("< 0", "0", "1 - 10", "11 - 50", "51 - 115")

# binning the specimen difference values
diff_specimen_raster <-
  diff_specimen_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

# set the levels and colours for plotting
levels <- levels(diff_specimen_raster)[[1]][, "ID"]
colours <- c("#4bb062", '#dddcdc', '#c0dceb', '#68abd0', '#2870b1')

# specimen difference map (fig 1c)
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

# save the specimen maps
if (save_figures) {
  ggsave(here("figures", "figure1a.svg"), plot = old_specimen_map)
  ggsave(here("figures", "figure1b.svg"), plot = new_specimen_map)
  ggsave(here("figures", "figure1c.svg"), plot = diff_specimen_map)
} else {
  old_specimen_map
  new_specimen_map
  diff_specimen_map
}

#### species rasters

# settings for binning the values in the raster cells
breaks <- c(1, 11, 20, 40, 70)
labels <- c("1 - 10", "11 - 20", "21 - 40", "41 - 70")

# converting 2007 species counts to rasters
old_species_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2007) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'species_count', project = project, xlim=raster_xlim, ylim=raster_ylim)

# binning the 2007 species counts
old_species_raster <- 
  old_species_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify

# converting 2017 species counts to rasters
new_species_cells <- 
  specimen_data %>%
  filter(!is.na(latitude), !is.na(longitude), dataset == 2017) %>%
  select(species, longitude, latitude) %>%
  convert_to_raster(cell_size = cell_size_km, 
                    method = 'species_count', project = project, xlim=raster_xlim, ylim=raster_ylim)

# binning the 2017 species counts
new_species_raster <- 
  new_species_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

# getting the levels from the rasters for plotting
levels <- levels(old_species_raster)[[1]][, "ID"]

# 2007 species map (fig 1d)
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

# 2017 species map (fig 1e)
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

# subtracting the two rasters above for the difference map
diff_species_cells <- subtract_rasters(old_species_cells, new_species_cells)

# settings for binning the differences
breaks <- c(-6, 0, 1, 2, 6, 11, 25)
labels <- c("< 0", "0", "1", "2 - 5", "6 - 10", "11 - 24")

# bin the differences
diff_species_raster <-
  diff_species_cells %>%
  cut(breaks = breaks, right = FALSE) %>%
  ratify()

# set the levels and colours for plotting
levels <- levels(diff_species_raster)[[1]][, "ID"]
colours <- c("#4bb062", '#dddcdc', '#d1e5f0', '#90c4dd', '#4393c3', '#2065ab')

# difference species map (fig 1f)
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

# save species maps
if (save_figures) {
  ggsave(here("figures", "figure1d.svg"), old_species_map)
  ggsave(here("figures", "figure1e.svg"), new_species_map)
  ggsave(here("figures", "figure1f.svg"), diff_species_map)
} else {
  old_species_map
  new_species_map
  diff_species_map
}

# figure 2 --------------------------------------------------------------------

# figure illustrating the change in conservation assessments
# figure 2a is a bar chart of species in each conservation assessment

# test significance of differences in number of species in assessment categories
modelled_status_summary <- 
  assessments %>%
  model_proportion_difference(dataset, eoo)

# test significant of differences in number of threatened species
modelled_threat_summary <-
  assessments %>%
  mutate(threatened = case_when(eoo == "DD" ~ "DD",
                                eoo > "NT" ~ "threatened",
                                TRUE ~ "not threatened"),
         threatened = factor(threatened)) %>%
  model_proportion_difference(dataset, threatened)

# count number of species in each assessment category for each snapshot
assessment_counts <-
  assessments %>%
  count(eoo, dataset)

# generate annotations for categories with significant changes
annotations <- 
  modelled_status_summary %>%
  filter(significant) %>%
  left_join(assessment_counts, by = "eoo") %>%
  group_by(eoo) %>%
  summarise(y = max(n) + 15) %>%
  mutate(significant = TRUE,
         x = as.numeric(eoo) - 0.3,
         xend = as.numeric(eoo) + 0.3,
         grouping = "in")

# generate annotation for significant change in threatened species
if (filter(modelled_threat_summary, threatened == "threatened")$significant == TRUE) {
  annotations <-
    annotations %>%
    bind_rows(tibble(eoo="threatened", 
                     x=3.7, 
                     xend=6.3, 
                     y=max(assessment_counts[assessment_counts$eoo > "NT", "n"]) + 30,
                     grouping="out"))
}
# generate annotation for significant change in not threatened species
if (filter(modelled_threat_summary, threatened == "not threatened")$significant == TRUE) {
  annotations <-
    annotations %>%
    bind_rows(tibble(eoo="not threatened", 
                     x=1.7, 
                     xend=3.3, 
                     y=max(assessment_counts[assessment_counts$eoo %in% c("NT", "LC"), "n"]) + 30,
                     grouping="out"))
}

# plot bar chart with annotations (fig 2a)
status_bar_chart <- 
  assessment_counts %>%
  ggbarplot(x = "eoo", y = "n", 
            fill = "dataset", color = "white", 
            position = position_dodge(width=0.75), label = "n", 
            ggtheme = theme_pubclean(), lab.size = 3.375,
            xlab = "Assessment category", ylab = "Number of species",
            lab.vjust = 0.5, lab.hjust = -0.1) +
  geom_segment(data = annotations, aes(x = x, xend = xend, y = y, yend = y, linetype=grouping)) +
  geom_text(data = annotations, aes(x = (x + xend) / 2, y = y + 4), label = "*", vjust=0.7) +
  coord_flip() +
  guides(linetype=FALSE) +
  scale_fill_brewer(palette = "Set1", name = "")

# save bar chart
if (save_figures) {
  ggsave(here("figures", "figure2a.svg"), height = 4.2, width = 8, units="in", plot = status_bar_chart)
} else {
  status_bar_chart
}

# figure 2b is a chord diagram showing movement of species between assessment categories

# count the number of species for each possible change
status_changes <- 
  species_status_changes %>%
  count(eoo_2007, eoo_2017) %>%
  complete(eoo_2007, eoo_2017, fill = list(n = 0)) %>%
  rename(preliminary = eoo_2007, 
         final = eoo_2017)

# open the device for saving the figure
if (save_figures) {
  svg(here("figures", "figure2b.svg"))
}

# plot the chord diagram
status_changes %>%
  plot_status_chord(label_outside = "NT")

# close the device for saving the figure
if (save_figures) {
  dev.off()
}

# figure 3 --------------------------------------------------------------------
# figure 3 is a scatter plot of the species EOO in 2007 against the EOO in 2017

eoo_comparison_plot <-
  species_status_changes %>%
  ggplot(mapping=aes(x=log(eoo_area_2007 + 1), y=log(eoo_area_2017 + 1))) +
  geom_point(colour="steelblue", alpha=0.5) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey50", size=1) +
  coord_equal() +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x=expression("log(" ~ EOO["2007"] ~ ")"), y=expression("log(" ~ EOO["2017"] ~ ")")) +
  theme(panel.border=element_blank())

if (save_figures) {
  ggsave(here("figures", "figure3.svg"), eoo_comparison_plot, width=5.39, height=5.39)
} else {
  eoo_comparison_plot
}

# figure 4 --------------------------------------------------------------------
# figure 4 shows the effect of each type of specimen change on species EOO

# only looking at specimen changes that could effect the EOO
filtered_changes <- 
  eoo_changes %>%
  # remove changes that can't effect EOO
  filter(!(change %in% c("nomenclatural", "no_change"))) %>%
  # create/order variables for easier plotting
  mutate(eoo_2017 = factor(eoo_2017, levels=categories, ordered=TRUE),
         eoo_2007 = factor(eoo_2007, levels=categories, ordered=TRUE),
         change = factor(change, levels=eoo_change_types, ordered=TRUE),
         moved = eoo_2017 != eoo_2007) %>%
  # remove species without an EOO in both time periods
  filter(!(eoo_2007 %in% "DD"), 
         !(eoo_2017 %in% "DD"), 
         !is.na(change), 
         !is.na(eoo_change),
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species, ]$species |
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species, ]$species,
         !is.na(moved))

# figure 4a is a bar chart of the proportion of each change

# get the proportion of specimen change types
change_proportions <-
  filtered_changes %>%
  count(change, moved) %>%
  group_by(moved) %>%
  mutate(proportion = n / sum(n))

# test significance of changes in the proportions
modelled_change_proportions <-
  filtered_changes %>%
  model_proportion_difference(moved, change)

# make annotations for significant changes
annotations <- 
  modelled_change_proportions %>%
  filter(significant) %>%
  left_join(change_proportions, by = "change") %>%
  group_by(change) %>%
  summarise(y = max(proportion) + 0.01) %>%
  mutate(significant = TRUE,
         x = as.numeric(change) - 0.3,
         xend = as.numeric(change) + 0.3)

# plot the bar chart of specimen change proportions (fig 4a)
proportions_by_type <- 
  change_proportions %>%
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

# figure 4b is a dot plot of the average change in EOO for each change type

# summarise change in EOO 
mean_changes <-
  filtered_changes %>%
  group_by(change) %>%
  summarise(mean = mean(eoo_change), 
            std_err = sd(eoo_change) / sqrt(n())) %>%
  mutate(too_few = change %in% c("downgraded_to_genus", "split"))

# make the dot plot (fig 4b)
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

# combine the two figures (fig 4)
combined_eoo_change_figs <- ggarrange(proportions_by_type, 
                                      mean_change_per_type + rremove("y.text"), 
                                      common.legend = T, 
                                      align = "h", 
                                      labels = c("(a)", "(b)"))

if (save_figures) {
  ggsave(here("figures", "figure4.svg"), combined_eoo_change_figs, width=9)
} else {
  combined_eoo_change_figs
}


# appendix s2 -----------------------------------------------------------------
# This is a histogram of the number of the collection date of new specimens

specimen_accumulation_plot <-
  specimen_changes %>%
  filter(change == "new_specimen") %>%
  filter(!is.na(collection_year_2017)) %>%
  mutate(collected = ifelse(collection_year_2017 <= 2007, "before 2007", "after 2007")) %>%
  ggplot(mapping=aes(x=collection_year_2017)) +
  geom_histogram(mapping=aes(fill=collected), binwidth=6, alpha=0.8, colour="black") +
  scale_fill_manual(values=c(`before 2007`="#636363", `after 2007`="#de2d26"), name="") +
  theme_pubclean() +
  labs(x="Collection year", y="Number of new specimens")

if (save_figures) {
  ggsave(here("figures", "appendix_s1.svg"), specimen_accumulation_plot, height=5, width=7)
} else {
  specimen_accumulation_plot
}

# appendix s3 -----------------------------------------------------------------
# Pie charts summarising the number of specimens with each change

# the first is the number of specimens overall with each type

# settings for the plot
specimen_types <- c("corrected_coordinates", "new_coordinates", 
                    "correction", "downgraded_to_genus", "upgraded_to_species", 
                    "split", "lumped", 
                    "new_specimen", "nomenclatural", "no_change")

type_colours <- c('#898989', '#d4d4d4', '#1f78b4', '#6a3d9a', '#cab2d6', 
                  '#ff7f00', '#fdbf6f', '#e31a1c', '#33a02c', '#b2df8a')

labels <- str_to_sentence(str_replace_all(rev(specimen_types), "_", " "))

# make the pie chart (appendix s3a)
all_specimen_changes_pie <-
  specimen_changes %>%
  filter(!is.na(change)) %>%
  plot_changes_pie(order = rev(specimen_types), labels=labels, colours=type_colours) +
  labs(title = "All Specimen Changes")+
  theme(axis.text.x = element_blank())

# the second is the number of specimens for each change that altered the EOO (appendix s3b)



eoo_type_colours <- c('#1f78b4', '#6a3d9a', '#cab2d6', '#ff7f00', '#fdbf6f', '#e31a1c', '#fb9a99', '#33a02c', '#b2df8a')

labels <- str_to_sentence(str_replace_all(rev(eoo_change_types), "_", " ")) 

all_eoo_change_specimens <-
  eoo_changes %>%
  mutate(moved = eoo_2007 != eoo_2017) %>%
  filter(!(eoo_2007 %in% "DD"), 
         !(eoo_2017 %in% "DD"), 
         !is.na(change), 
         !is.na(eoo_change),
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species, ]$species |
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species, ]$species,
         !is.na(moved)) %>%
  plot_changes_pie(order = rev(eoo_change_types), labels=labels, colours=eoo_type_colours) +
  labs(title = "Specimen Changes Used For EOO Change Calculation") +
  theme(axis.text.x = element_blank())

# the third is the number of specimens for each change that altered the assessmnet category (appendix s3c)

changed_eoo_specimens <-
  eoo_changes %>%
  mutate(moved = eoo_2007 != eoo_2017) %>%
  filter(!(eoo_2007 %in% "DD"), 
         !(eoo_2017 %in% "DD"), 
         !is.na(change), 
         !is.na(eoo_change),
         assessed_species_2007 %in% specimen_data[specimen_data$determined_to_species, ]$species |
           assessed_species_2017 %in% specimen_data[specimen_data$determined_to_species, ]$species,
         !is.na(moved)) %>%
  filter(moved == TRUE) %>%
  plot_changes_pie(order = rev(eoo_change_types), labels=labels, colours = eoo_type_colours) +
  labs(title = "Specimen Changes For Species that Changed Conservation Status") +
  theme(axis.text.x = element_blank())


if (save_figures) {
  ggsave(here("figures", "appendix_s3a.svg"), all_specimen_changes_pie)
  ggsave(here("figures", "appendix_s3b.svg"), all_eoo_change_specimens)
  ggsave(here("figures", "appendix_s3c.svg"), changed_eoo_specimens)
} else {
  all_specimen_changes_pie
  all_eoo_change_specimens
  changed_eoo_specimens
}

# appendix s4 -----------------------------------------------------------------
# a bar chart of the Sorensen dissimilarity in assessment categories

# calculate Sorensen dissimilarity
assessment_turnover <-
  species_status_changes %>%
  mutate(species = 1:nrow(.)) %>%
  filter(!is.na(eoo_2007), !is.na(eoo_2017)) %>%
  select(species, eoo_2007, eoo_2017) %>%
  gather(dataset, eoo, -species) %>%
  mutate(dataset = str_replace(dataset, "eoo_", "")) %>%
  group_by(eoo) %>%
  nest() %>%
  mutate(turnover=map_dbl(data, calculate_turnover)) %>%
  select(-data)

# make bar chart (appendix s4)
turnover_plot <- 
  ggplot(data=assessment_turnover, mapping=aes(x=eoo, y=turnover)) +
  geom_col(fill="steelblue") +
  geom_text(mapping=aes(label=sprintf("%.2f", turnover)), vjust=-0.25) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x="Assessment category", y="Dissimilarity") +
  theme_pubclean()

if (save_figures) {
  ggsave(here("figures/appendix_s4.svg"), turnover_plot, width=7, height=4)
} else {
  turnover_plot
}

# appendix s5 -----------------------------------------------------------------
# dot plots summarising the change in EOO for the specimen categories in S3

mean_change_per_type_split <-
  filtered_changes %>%
  group_by(change, moved) %>%
  summarise(mean = mean(eoo_change), 
            std_err = sd(eoo_change) / sqrt(n())) %>%
  mutate(moved = ifelse(moved, "Status change", "No status change")) %>%
  ungroup() %>%
  bind_rows(mean_changes %>% mutate(moved = "Overall") %>% select(-too_few)) %>%
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

if (save_figures) {
  ggsave(here("figures/appendix_s5.svg"), mean_change_per_type_split, width=7, height=9)
} else {
  mean_change_per_type_split
}

# record accumulation ---------------------------------------------------------

specimen_accumulation_plot <-
  specimen_data %>%
  mutate(collected = ifelse(collection_year <= 2007, "before 2007", "after 2007")) %>%
  ggplot(mapping=aes(x=collection_year)) +
  geom_histogram(mapping=aes(fill=collected), binwidth=6, alpha=0.8, colour="black") +
  scale_fill_manual(values=c(`before 2007`="#636363", `after 2007`="#de2d26"), name="") +
  theme_pubclean() +
  labs(x="Collection year", y="Number of new specimens")

if (save_figures) {
  ggsave(here("figures", "appendix_s1.svg"), specimen_accumulation_plot, height=5, width=7)
} else {
  specimen_accumulation_plot
}

# specimen turnover -----------------------------------------------------------

assessment_turnover <-
  species_status_changes %>%
  mutate(species = 1:nrow(.)) %>%
  filter(!is.na(eoo_2007), !is.na(eoo_2017)) %>%
  select(species, eoo_2007, eoo_2017) %>%
  gather(dataset, eoo, -species) %>%
  mutate(dataset = str_replace(dataset, "eoo_", "")) %>%
  group_by(eoo) %>%
  nest() %>%
  mutate(turnover=map_dbl(data, calculate_turnover)) %>%
  select(-data)

turnover_plot <- 
  ggplot(data=assessment_turnover, mapping=aes(x=eoo, y=turnover)) +
  geom_col(fill="steelblue") +
  geom_text(mapping=aes(label=sprintf("%.2f", turnover)), vjust=-0.25) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x="Assessment category", y="Dissimilarity") +
  theme_pubclean()

if (save_figures) {
  ggsave(here("figures/appendix_s4.svg"), turnover_plot, width=7, height=4)
} else {
  turnover_plot
}

