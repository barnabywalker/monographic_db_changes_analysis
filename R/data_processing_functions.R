###############################################################################
# Functions for loading and cleaning data from the database.                  #
###############################################################################

library(here)
library(RODBC)
library(infuser)
library(biogeo)

clean_data <- function(data) {
  #` Clean data from the Myrcia database
  #`
  #` Specific sequence of steps to clean Myrcia database data,
  #' starting by renaming columns, then removing unnecessary words
  #' from conservation status, then calculating coordinates in
  #' decimal degrees from DMS.
  #' 
  #' @param data A data frame of data from the Myrcia database.
  #' 
  #' @return A data frame of cleaned data from the Myrcia database.
  
  data <- data %>% 
    rename_all(. %>% gsub("Conservation ", "", .)) %>%
    rename_all(. %>% gsub(" ", "_", .))
  
  if (!("AOO_Rating" %in% colnames(data))) {
    data$AOO_Rating <- str_match(data$AOO_Rating, "\\(([A-Z]+)\\)")[,2]
    data$EOO_Rating <- str_match(data$EOO_Rating, "\\(([A-Z]+)\\)")[,2]
    data$Final_Assessment <- str_match(data$Final_Assessment, "\\(([A-Z]+)\\)")[,2]
  }
  
  data <-
    data %>%
    mutate(Longit_Secs = as.numeric(Longit_Secs),
           Latit_Secs = as.numeric(Latit_Secs)) %>%
    mutate(Latit_Secs = ifelse(is.na(Latit_Secs), 0, Latit_Secs),
           Longit_Secs = ifelse(is.na(Longit_Secs), 0, Longit_Secs)) %>%
    mutate(Longitude = "W") %>%
    mutate(LongDD = dms2dd(Longitud_Degree, Longitu_Minutes, Longit_Secs, Longitude),
           LatDD = dms2dd(Latitud_Degree, Latitu_Minutes, Latit_Secs, Latitude)) %>%
    mutate(LongDD = ifelse(LongDD == 0 & LatDD == 0, NA, LongDD),
           LatDD = ifelse(LongDD == 0 & LatDD == 0, NA, LatDD))

  data <- 
    data %>%
    mutate(Checklist = ifelse(Checklist == "Gomidesia cordiifolia" & Species == "affinis", "Gomidesia affinis", Checklist)) %>%
    mutate(Checklist = ifelse(Checklist == "Gomidesia grandifolia" & Species == "crocea", "Gomidesia crocea", Checklist))
  
  return(data)
}

load_data <- function(table, from_database = FALSE, directory = here("data", "dumps"), query = NULL, filename = NULL) {
  #` Load data from the specified table in the Myrcia database.
  #`
  #' Load data from the specified table, either from the Myrcia database
  #' or from a csv file of the most recent dump from that table. An SQL
  #' query can be specified for more control over the data loaded.
  #' 
  #' @param table A string, the name of the table to load from.
  #' @param from_database A boolean, whether to load from the database or the latest dump.
  #' @param directory The directory to load the latest dump from.
  #' @param query A string of a custom SQL query to use.
  #' @param filename The name of the dump file to load from.
  #' 
  #' @return A data frame of the loaded data.
  
  if (from_database) {
    db_path <- "T:/GIS/Eve/Myrtaceae_Main_2018.accdb"
    driver <- "Microsoft Access Driver (*.mdb, *.accdb)"
    conn_str <- sprintf("Driver={%s};DBQ=%s", driver, db_path)
    cnxn <- odbcDriverConnect(conn_str)
    
    if (is.null(query)) {
      query <- "SELECT *
      FROM {{table}}
      WHERE Checklist IS NOT NULL
      AND `Not to include` = 0
      AND
      (
      Checklist LIKE 'Myrcia %' 
      OR Checklist LIKE 'Marlierea %' 
      OR Checklist LIKE 'Gomidesia %' 
      OR Checklist LIKE 'Calyptranthes %' 
      OR Checklist LIKE 'Mitranthes %'
      OR Checklist LIKE 'Myrcia' 
      OR Checklist LIKE 'Marlierea' 
      OR Checklist LIKE 'Gomidesia' 
      OR Checklist LIKE 'Calyptranthes' 
      OR Checklist LIKE 'Mitranthes'
      )"
      }
    
    query <- infuse(query, table = table)
    data <- sqlQuery(cnxn, gsub("\n", " ", query), stringsAsFactors = FALSE)
    odbcClose(cnxn)
    
    } else {
      
      if (is.null(filename)) {
        files <- list.files(directory, pattern = paste(table, ".+\\.csv", sep=""))
        if (length(files) > 0) {
          filename <- files[length(files)]
        } else {
          stop("No data file exists for that table, query from database or call 'dumpData' to create data file")  
        }
      } else {
        if (!file.exists(paste(directory, filename))) {
          stop("No file found with that name")
        }
      }
      
      data <- read_csv(paste(directory, files[length(files)], sep=""))
    }
  
  return(data)
  }

dump_data <- function(table, directory = here("data", "dumps"), filename = NULL, query = NULL) {
  #' Dumps data from a table in the Myrcia database to a file.
  #' 
  #' @param table A string, the name of the table to load from.
  #' @param directory The directory to dump the data to.
  #' @param query A string of a custom SQL query to use.
  #' @param filename The name of the file to dump to.
  
  dir.create(directory)
  data <- load_data(table, from_database = TRUE, query = query)
  if (is.null(filename)) {
    date <- format(Sys.Date(), "%Y%m%d")
    filename <- file.path(directory, paste(table, "_", date, ".csv", sep=""))
  } else {
    filename <- file.path(directory, filename)
  }
  write_csv(data, filename)
  sprintf("Wrote date from SQL table '%s' to file '%s'", table, filename)
}

