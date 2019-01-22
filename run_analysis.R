# Script to run analysis of myrcia data for paper
library(here)

# settings that can be changed
save_figures = TRUE
specimen_data_file <- here("data", "specimen_data.csv")
name_change_data_file <- here("data", "name_change_data.csv")

# run data transformation script
source(here("analysis", "transform_data.R"))

# make plots
source(here("analysis", "plot_data.R"))

# generate statistics table
source(here("analysis", "summarise_data.R"))
