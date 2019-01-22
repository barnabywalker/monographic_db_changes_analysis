###############################################################################
#      Script to run the total analysis of the Myrcia specimen database       #
#                                                                             #
# Runs the analysis from end to end on the specimen data extracted from our   #
# Myrcia database. Will reproduce results and figures included in the paper.  #
###############################################################################
library(here)

# changeable variables --------------------------------------------------------

# Change these values to point to where you have the data files stored
specimen_data_file <- here("data", "specimen_data.csv")
name_change_data_file <- here("data", "name_change_data.csv")

# Change this to switch on or off figure saving
save_figures = TRUE

# running the script ----------------------------------------------------

# run data transformation script
source(here("analysis", "transform_data.R"))

# make plots
source(here("analysis", "plot_data.R"))

# generate statistics table
source(here("analysis", "summarise_data.R"))
