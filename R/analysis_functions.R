###############################################################################
# Functions used to analyse the data from the specimen database, and the      #
# generated data like conservation assessments.                               #
###############################################################################

library(rCAT)
library(gtools)
library(raster)
library(HDInterval)
library(magrittr)

convert_to_raster <- function(species_occurrences, 
                              cell_size, 
                              method=c("species_count", "specimen_count"), 
                              project=FALSE, 
                              crs=CRS("+proj=wag4 +units=m"), 
                              xlim=c(-180, 180), 
                              ylim=c(-90, 90)) {
  #' Convert a data frame of species occurrences to a raster.
  #' 
  #' Converts a data frame of species occurrences to a Raster object
  #' using the specified cell size (in km) either aggregating by species
  #' (method = "species_count") or specimen (method = "specimen_count").
  #' Also converts to a specified projection, default being an equal area
  #' one.
  #' 
  #' @param species_occurrences A data frame of occurrences for one or more species.
  #' 
  #' @param cell_size The size of cell to use in the Raster, in km.
  #' 
  #' @param method The method for aggregating the occurrences, by species ("species_count")
  #' or by specimen ("specimen_count").
  #' 
  #' @param project Whether to reproject the occurrence data or not.
  #' 
  #' @param crs The project to reproject the data to.
  #' 
  #' @param xlim A vector containing the x-direction limits (longitude) of the Raster,
  #' in decimal degrees.
  #' 
  #' @param ylim A vector containing the y-direction limits (latitude) of the Raster,
  #' in decimal degrees.
  #' 
  #' @return A Raster Object with counts of species occurrence as cell values.
  
  default_crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  method <- match.arg(method)
  
  # generate raster to right size
  if (project) {
    rast <- raster(crs=crs)
    
    # reproject given limits
    zoom_area <- c(xlim, ylim) %>%
      extent %>%
      as("SpatialPolygons")
    
    proj4string(zoom_area) <- default_crs
    extent(rast) <- extent(spTransform(zoom_area, rast@crs))
    cell_size <- cell_size * 1000
  } else {
    rast <- raster()
    extent(rast) <- c(xlim, ylim)
    res(rast) <- 1
    cell_conversion <- sqrt(mean(values(area(rast))))
    cell_size <- cell_size / cell_conversion
  }
  res(rast) <- cell_size
  
  split_occurrences <- split(species_occurrences, species_occurrences$species, drop = TRUE)
  
  occurrences_projection <- lapply(lapply(split_occurrences, '[', c("longitude","latitude")), SpatialPoints, proj4string = default_crs)
  if (project){
    occurrences_projection <- lapply(occurrences_projection, spTransform, crs)
  }
  
  cell_values <- switch(method, 
                        species_count = aggregate_by_species_count(rast, occurrences_projection), 
                        specimen_count = aggregate_by_specimen_count(rast, occurrences_projection))
  
  return(cell_values)
}

aggregate_by_species_count <- function(rast, occurrences) {
  #' Aggregate occurrences for raster cells by species.
  #' 
  #' @param rast The raster to aggregate the occurrences into.
  #' 
  #' @param occurrences A list of occurrences.
  #' 
  #' @return A raster object with species count values in the cells.
  
  cell_positions <- list()
  for (i in 1:length(occurrences)) {
    cell_positions[[i]] <- unique(cellFromXY(rast, xy = occurrences[[i]]))
  }
  
  cell_counts <- rle(sort(unlist(cell_positions)))
  for (i in 1:max(cell_counts$lengths)) {
    rast[cell_counts$values[cell_counts$lengths == i]] <- i
  }
  
  return(rast)
}

aggregate_by_specimen_count <- function(rast, occurrences) {
  #' Aggregate occurrences for raster cells by specimen.
  #' 
  #' @param rast The raster to aggregate the occurrences into.
  #' 
  #' @param occurrences A list of occurrences.
  #' 
  #' @return A raster object with specimen count values in the cells.
  
  cell_positions <- list()
  for (i in 1:length(occurrences)) {
    cell_positions[[i]] <- cellFromXY(rast, xy = occurrences[[i]])
  }
  
  cell_counts <- rle(sort(unlist(cell_positions)))
  for (i in 1:max(cell_counts$lengths)) {
    rast[cell_counts$values[cell_counts$lengths == i]] <- i
  }
  
  return(rast)
}

subtract_rasters <- function(raster1, raster2) {
  #' Subtract two rasters from each other.
  #' 
  #' Normal raster subtraction results in NA if one
  #' of the cell values is NA. This subtracts two rasters
  #' cell by cell, but taking NA values as zero when
  #' something is subtracted from them.
  #' 
  #' @param raster1 The raster to subtract.
  #' 
  #' @param raster2 The raster to subtract from.
  #' 
  #' @return A raster containing the difference between the two.
  
  na_idx1 <- Which(is.na(raster1))
  na_idx2 <- Which(is.na(raster2))
  raster1[na_idx1] <- 0
  raster2[na_idx2] <- 0
  diff_raster <- raster2 - raster1
  diff_raster[na_idx1 & na_idx2] <- NA
  return(diff_raster)
}



model_binomial_difference <- function(data, filter_var, category_var, n_samples, comparison_value=0) {
  #' Model the difference in the probability of getting one of two outcomes between
  #' two datasets.
  #' 
  #' Models the given data as a binomial process and simulates the posterior distribution of
  #' the difference between two datasets. This uses a beta
  #' prior with equal probabilities, and updates it with the two category counts from
  #' the data. It then draws the specified number of samples from each posterior and
  #' subtracts the two distributions to get the posterior distribution of the difference.
  #' This is summaries by the mean value of the posterior difference, and the 95 % credible
  #' interval on the difference. The change is judged significant if 0 does not appear
  #' in the 95% CI.
  #' 
  #' @param data The observed values to be modelled.
  #' 
  #' @param filter_var The variable defining the two datasets to find a difference between.
  #' 
  #' @param category_var The variable that has the two categories to use in the binomial model,
  #' c.f. "successes" and "failures".
  #' 
  #' @param n_samples The number of samples to draw from the posterior distribution.
  #' 
  #' @return A data frame summarising the difference, with the mean posterior difference,
  #' the 95% CI, and if the difference is significant or not.
  
  filter_var <- enquo(filter_var)
  category_var <- enquo(category_var)
  
  filter_values <- 
    data %>%
    pull(!!filter_var) %>%
    unique()
  
  summary <-
    data %>%
    group_by(!!filter_var, !!category_var) %>%
    summarise(n = n()) %>%
    ungroup()
  
  
  set_a <-
    summary %>%
    filter((!!filter_var) == filter_values[1]) %>%
    spread(!!category_var, n)

  set_b <-
    summary %>%
    filter((!!filter_var) == filter_values[2]) %>%
    spread(!!category_var, n)
  
  posterior_samples_a <- rbeta(n_samples, 1+set_a$`TRUE`, 1+set_a$`FALSE`)
  posterior_samples_b <- rbeta(n_samples, 1+set_b$`TRUE`, 1+set_b$`FALSE`)
  
  posterior_summary <-
    tibble(`2007` = posterior_samples_a, `2017` = posterior_samples_b) %>%
    mutate(sample = 1:nrow(.)) %>%
    gather(dataset, p, -sample) %>%
    group_by(sample) %>%
    summarise(difference = diff(p)) %>%
    summarise(mean_diff = mean(difference),
              low_ci = hdi(difference, 0.95)[1],
              high_ci = hdi(difference, 0.95)[2]) %>%
    mutate(significant = (comparison_value < low_ci) | (comparison_value > high_ci))
 
  return(posterior_summary) 
}

model_proportion_difference <- function(data, 
                                        filter_var, 
                                        category_var, 
                                        n_samples=10000,
                                        comparison_value=0) {
  #' Model the difference in proportions of data in two sets
  #' 
  #' Use separate multinomial models of two sets of data to
  #' first draw samples from Dirichlet posteriors, and then
  #' summarise the differences in distributions.
  #' 
  #' Model the difference in the probability of getting each of number of outcomes between
  #' two datasets.
  #' 
  #' Models the given data as a multinomial process and simulates the posterior distribution of
  #' the difference between two datasets. This uses a dirichlet
  #' prior with equal probabilities, and updates it with the category counts from
  #' the data. It then draws the specified number of samples from each posterior and
  #' subtracts the two distributions to get the posterior distribution of the difference.
  #' This is summaries by the mean value of the posterior difference, and the 95 % credible
  #' interval on the difference. The change is judged significant if 0 does not appear
  #' in the 95% CI.
  #' 
  #' @param data The observed values to be modelled.
  #' 
  #' @param filter_var The variable defining the two datasets to find a difference between.
  #' 
  #' @param category_var The variable containing the categories to model.
  #' 
  #' @param n_samples The number of samples to draw from the posterior distribution.
  #' 
  #' @return A data frame summarising the difference, with the mean posterior difference,
  #' the 95% CI, and if the difference is significant or not.
  
  filter_var <- enquo(filter_var)
  category_var <- enquo(category_var)
  
  filter_name <- quo_name(filter_var)
  category_name <- quo_name(category_var)
  
  filter_values <- 
    data %>%
    pull(!!filter_var) %>%
    unique()
  category_values <- 
    data %>%
    pull(!!category_var) %>%
    levels()
  
  posterior_samples_a <- 
    data %>%
    filter((!!filter_var) == filter_values[1]) %>%
    group_by(!!category_var) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    draw_dirichlet_posterior(categories = category_values, n_samples) %>%
    mutate(!!filter_name := filter_values[1],
           !!category_name := factor(category, levels = category_values, ordered = TRUE))
  
  posterior_samples_b <- 
    data %>%
    filter((!!filter_var) == filter_values[2]) %>%
    group_by(!!category_var) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    draw_dirichlet_posterior(categories = category_values, n_samples) %>%
    mutate(!!filter_name := filter_values[2],
           !!category_name := factor(category, levels = category_values, ordered = TRUE))
  
  posterior_summary <-
    rbind(posterior_samples_a, posterior_samples_b) %>%
    group_by(sample, !!category_var) %>%
    summarise(difference = diff(p)) %>%
    group_by(!!category_var) %>%
    summarise(mean_diff = mean(difference),
              low_ci = hdi(difference, 0.95)[1],
              high_ci = hdi(difference, 0.95)[2]) %>%
    mutate(significant = (comparison_value < low_ci) | (comparison_value > high_ci))
  
  return(posterior_summary)
  
}


draw_dirichlet_posterior <- function(data, categories, n_samples) {
  #' Draw samples from a Dirichlet posterior.
  #' 
  #' First creates the Dirichlet prior with equal probability of
  #' getting each category, then updates it with the data to 
  #' get the posterior, then draws samples and shapes into a
  #' data frame.
  #' 
  #' @param data The data to model.
  #' 
  #' @param categories A vector of the categories in the data.
  #' 
  #' @param n_samples The number of samples to draw from the
  #' posterior.
  #' 
  #' @return A data frame of the samples from the posterior.
  
  posterior_samples <- data %>%
    mutate(n = n + 1) %>%
    pull(n) %>%
    rdirichlet(n_samples, .) %>%
    set_colnames(categories) %>%
    as.tibble() %>%
    mutate(sample = seq(1, nrow(.))) %>%
    gather(category, p, -sample)
  
  return(posterior_samples)
}

calculate_turnover <- function(data) {
  #' Calculate turnover in sets of species between 2017 and 2007.
  #' 
  #' Calculates turnover as the Sorensen distance for sets of species,
  #' based on the assumption that there is one set of species for 2007
  #' and one for 2017.
  #' 
  #' @param data A dataframe with a column for species and one for
  #' the dataset (if the species belongs to 2007 or 2017).
  #' @return a number between 0 and 1, the dissimilarity.
  
  species_before <- filter(data, dataset == 2007)$species
  species_after <- filter(data, dataset == 2017)$species
  
  unique_before <- length(setdiff(species_before, species_after))
  unique_after <- length(setdiff(species_after, species_before))
  in_both <- length(intersect(species_after, species_before))
  
  (unique_before + unique_after) / (2*in_both + unique_before + unique_after)
}
