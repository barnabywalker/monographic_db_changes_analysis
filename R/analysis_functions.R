###############################################################################
# Functions used to analyse the data from the specimen database, and the      #
# generated data like conservation assessments.                               #
###############################################################################

library(rCAT)
library(gtools)
library(raster)

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



model_binomial_difference <- function(data, filter_var, category_var, n_samples) {
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
              low_ci = quantile(difference, 0.025, names=FALSE),
              high_ci = quantile(difference, 0.975, names=FALSE)) %>%
    mutate(significant = (0 < low_ci) | (0 > high_ci))
 
  return(posterior_summary) 
}

model_proportion_difference <- function(data, 
                                        filter_var, 
                                        category_var, 
                                        n_samples = 10000) {
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
              low_ci = quantile(difference, 0.025, names=FALSE),
              high_ci = quantile(difference, 0.975, names=FALSE)) %>%
    mutate(significant = (0 < low_ci) | (0 > high_ci))
  
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


model_status_numbers <- function(assessments) {
  #' model proportions of data with particular statuses in two data sets
  #' 
  #' TODO: maybe delete?
  
  posterior_samples_2007 <- assessments %>%
    filter(dataset == 2007) %>%
    count(eoo) %>%
    draw_dirichlet_posterior(categories = levels(assessments$eoo), 10000) %>%
    mutate(dataset = 2007,
           status = factor(category, levels=levels(assessments$eoo), ordered=TRUE))
  
  posterior_samples_2017 <- assessments %>%
    filter(dataset == 2017) %>%
    count(eoo) %>%
    draw_dirichlet_posterior(categories = levels(assessments$eoo), 10000) %>%
    mutate(dataset = 2017,
           status = factor(category, levels=levels(assessments$eoo), ordered=TRUE))
  
  posterior_summary <- rbind(posterior_samples_2007, posterior_samples_2017) %>%
    spread(dataset, p) %>%
    mutate(difference = `2017` - `2007`) %>%
    group_by(status) %>%
    summarise(mean_diff = mean(difference),
              low_ci = quantile(difference, 0.025, names=FALSE), 
              high_ci = quantile(difference, 0.975, names=FALSE)) %>%
    mutate(significant = (0 < low_ci) | (0 > high_ci),
           status = factor(status, levels=levels(assessments$eoo), ordered=TRUE))
  
  return(posterior_summary)
}

wilcoxon_pvalues <- function(contributions) {
  #' calculate p-values using Wilcoxon's test
  
  contributions %>%
    group_by(change) %>%
    summarise(pvalue = wilcox.test(eoo_change ~ moved)$p.value)
}

permutation_pvalues <- function(contributions, statistic = "mean") {
  #' calculate p-values by permutation
  
  obs_diffs <-
    contributions %>%
    group_by(change, moved) %>%
    summarise(stat = switch(statistic,
                            mean = mean(eoo_change),
                            median = median(eoo_change))) %>%
    spread(moved, stat) %>%
    mutate(obs_diff = `FALSE` - `TRUE`) %>%
    select(change, obs_diff)
  
  shuffled_contributions <-
    do(5000) *
    (contributions %>%
       group_by(change) %>%
       mutate(moved = shuffle(moved)) %>%
       group_by(change, moved) %>%
       summarise(stat = switch(statistic,
                               mean = mean(eoo_change),
                               median = median(eoo_change))))
  
  null_diffs <-
    shuffled_contributions %>%
    group_by(.index, change) %>%
    summarise(diffs = diff(stat))
  
  null_diffs %>%
    left_join(obs_diffs, by = c("change" = "change")) %>%
    mutate(null_gt_obs = abs(diffs) >= abs(obs_diff)) %>%
    group_by(change) %>%
    summarise(pvalue = mean(null_gt_obs))
  
}

test_differences_significance <- function(contributions, 
                                          test=c("wilcoxon", "permutation"), 
                                          stat="mean", level=0.05) {
  #' test for a significant difference between two data sets using the
  #' specified test
   
  pvalues <- switch(test,
                    wilcoxon = wilcoxon_pvalues(contributions),
                    permutation = permutation_pvalues(contributions, stat = "mean"))
  
  pvalues %>%
    mutate(significant = pvalue <= level)
}

permutation_test <- function(data, group_var, measure_var, statistic=c("mean", "median")) {
  #' test for significance using permutation
  
  if (is_character(group_var)) {
    group_var <- as.name(group_var)
  } else {
    group_var <- enquo(group_var)
  }
  
  if (is_character(measure_var)) {
    measure_var <- as.name(measure_var)
  } else {
    measure_var <- enquo(measure_var)  
  }
  
  obs_diffs <-
    data %>%
    group_by(!!group_var) %>%
    summarise(stat = switch(statistic,
                            mean = mean(!!measure_var),
                            median = median(!!measure_var))) %>%
    pull(stat) %>%
    diff()
  
  shuffled_contributions <-
    do(5000) *
    (data %>%
       mutate(label = shuffle((!!group_var))) %>%
       group_by(label) %>%
       summarise(stat = switch(statistic,
                               mean = mean(!!measure_var),
                               median = median(!!measure_var))))
  
  null_dist <- 
    shuffled_contributions %>%
    group_by(.index) %>%
    summarise(diff = diff(stat))
  
  null_dist %>%
    mutate(null_gt_obs = abs(diff) >= abs(obs_diff)) %>%
    pull(null_gt_obs) %>%
    mean()
}

permute_change_pair <- function(pair, data, statistic) {
  #' carry out a permutation significance test for a pair in data
  
  data <- 
    data %>% 
    filter(change %in% c(as.character(pair$type1), as.character(pair$type2)))
  
  permutation_test(data, "change", "eoo_change", statistic = statistic)
}

test_pairwise_permutations <- function(contributions, statistic = c("mean", "median")) {
  #' test for significant differences pairwise, by permutation
  
  combinations <-
    crossing(contributions$change, contributions$change) %>%
    set_colnames(c("type1", "type2")) %>%
    filter(type1 != type2)
  
  combination_list <- split(combinations, seq(nrow(combinations)))
  p_values <- map(combination_list, function(x) permute_change_pair(x, contributions, statistic))
  
  tibble(type1 = combinations$type1, type2 = combinations$type2, pvalue = p_values)
  
}
