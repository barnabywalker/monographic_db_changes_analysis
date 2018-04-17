###############################################################################
# Convenience functions to generate some of the plots used.                   #
###############################################################################

library(circlize)
library(maptools)
library(ggsn)

ggplot_raster_map <- function(raster_cells, 
                              raster_cells2 = NULL, 
                              xlim=c(-130, -20), 
                              ylim=c(-35, 35), 
                              categorical=FALSE, 
                              levels = NULL, 
                              labels = NULL,
                              scale.size = NULL,
                              north.symbol = NULL,
                              north.scale = 0.2,
                              scale.text.size = 5) {
  #' Plot a map of one or two rasters using ggplot.
  #' 
  #' @param raster_cells The raster to plot.
  #' @param raster_cells2 An optional second raster to plot.
  #' @param xlim The longitudinal limits in decimal degrees.
  #' @param ylim The latitudinal limits in decimal degrees.
  #' @param categorical Whether the raster values are categorical or not.
  #' @param levels The levels of a categorical raster.
  #' @param labels The labels of the levels.
  #' @param scale.size The size of map scale to add.
  #' @param north.symbol The type of north symbol to add.
  #' @param north.scale The size of the north symbol.
  #' @param scale.text.size The size of the scale text.
  #' 
  #' @return A ggplot object of the raster over a map.
  
  raster_polygons <- rasterToPolygons(raster_cells)
  default_crs <- raster()@crs
  data_crs <- raster_cells@crs
  
  # load world map
  data(wrld_simpl)
  
  if (data_crs@projargs != default_crs@projargs) {
    # calculate extent of area of interest in crs
    zoom_area <- extent(c(xlim, ylim))
    e <- as(zoom_area, "SpatialPolygons")
    proj4string(e) <- default_crs
    zoom_area <- extent(spTransform(e, data_crs))
    
    world_map <- spTransform(wrld_simpl, data_crs)
    xlim = c(zoom_area@xmin, zoom_area@xmax)
    ylim = c(zoom_area@ymin, zoom_area@ymax)
  } else {
    world_map <- wrld_simpl
  }
  
  if (is.null(labels)) {
    labels <- levels
  }
  
  
  # convert polygons to dataframe for ggplot to handle
  raster_polygons@data$id <- rownames(raster_polygons@data)
  polygon_points <- fortify(raster_polygons, region = "id")
  polygon_df <- merge(polygon_points, raster_polygons@data, by = "id")
  
  if (!is.null(raster_cells2)) {
    raster_polygons <- rasterToPolygons(raster_cells2)
    raster_polygons@data$id <- rownames(raster_polygons@data)
    polygon_df %>% 
      mutate(dataset = "2007") %>%
      rbind(merge(raster_polygons %>%
                    fortify(region = "id"), 
                  raster_polygons@data, by = "id") %>%
              mutate(dataset = "2007"))
  }
  
  # plot map
  p <- ggplot() + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), alpha=0.5, fill=NA, colour="black")
  # add cells
  if (categorical) {
    if (!is.null(levels)) {
      p <- p + geom_polygon(data=polygon_df, aes(x=long, y=lat, fill=factor(layer, levels = levels, labels = labels), group=group))
    } else {
      p <- p + geom_polygon(data=polygon_df, aes(x=long, y=lat, fill=factor(layer), group=group))
    }
  } else {
    p <- p + geom_polygon(data=polygon_df, aes(x=long, y=lat, fill=layer, group=group), colour="black")
  }
  # format plot
  p <- p + coord_equal(xlim=xlim, ylim=ylim) + theme_bw()
  p <- p + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
                 axis.text.x=element_blank(), axis.text.y=element_blank(), 
                 legend.position=c(0.05,0.23))
  p <- p + labs(x="", y="")
  
  if (!is.null(raster_cells2)) {
    p <- p + facet_wrap(~ dataset)
  }
  
  if (!is.null(scale.size)) {
    p <- p + ggsn::scalebar(x.min=xlim[1], x.max=xlim[2], 
                            y.min=ylim[1], y.max=ylim[2], 
                            dist=scale.size, st.size=scale.text.size, 
                            height=0.01)
  }
  
  if (!is.null(north.symbol)) {
    p <- p + north(x.min=xlim[1], x.max=xlim[2],
                   y.min=ylim[1], y.max=ylim[2], 
                   symbol=north.symbol, scale=north.scale)
  }
  
  return(p)
  
}


plot_status_chord <- function(status_changes, colour_palette=NULL, highlight_status=NULL, label_outside=c()) {
  #' Plot the change in conservation statuses as a chord diagram.
  #' 
  #' @param status_changes A data frame containing the count of each conservation status change.
  #' @param colour_palette The (optional) colour palette to use.
  #' @param highlight_status The name of a conservation status to highlight.
  #' @param label_outside A vector of conservation statuses that need their label plotting outside their segment.
  #' 
  #' @return A base R plot of a chord diagram.
  
  transparency <- ifelse(status_changes$preliminary == status_changes$final, 1, 0.5)
  if (is.null(colour_palette)){
    colour_palette <- brewer.pal(7, "YlOrRd")[2:7]
  }
  
  if (!is.null(highlight_status)) {
    colours <- c("grey50" %>%
                   c %>%
                   rep(6),
                 rev(colour_palette))
    
    status_idx <- which(highlight_status == levels(status_changes$preliminary))
    colours[status_idx] <- colour_palette[status_idx]
    
  } else {
    colours <- c(colour_palette, rev(colour_palette))
  }
  
  to_plot <- status_changes %>%
    mutate(preliminary = paste(as.character(preliminary), "_P", sep=""),
           final = paste(as.character(final), "_F", sep=""))
  
  row_sum <- to_plot %>% group_by(preliminary) %>% summarise(n = sum(n)) %>% pull(n) %>% sum()
  col_sum <- to_plot %>% group_by(final) %>% summarise(n = sum(n)) %>% pull(n) %>% sum()
  
  small_gap <- 1
  big_gap <- 20
  
  nr <- 6
  nc <- 6
  
  n_sector <- nr + nc
  row_sector_degree <- (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) + small_gap*(nr-1)
  
  start_degree = 90 - (180 - row_sector_degree)/2 + 180
  
  gaps = c(rep(small_gap, nr - 1), big_gap, rep(small_gap, nc - 1), big_gap)
  circos.par(gap.after = gaps, start.degree = start_degree)
  
  to_plot %>%
    chordDiagram(transparency = transparency, grid.col = colours, annotationTrack = "grid",
                 order = c("DD_P", "LC_P", "NT_P", "VU_P", "EN_P", "CR_P", 
                           "CR_F", "EN_F", "VU_F", "NT_F", "LC_F", "DD_F"),
                 annotationTrackHeight = uh(5, "mm"),
                 link.sort = TRUE,
                 link.decreasing = TRUE)
  
  for(si in get.all.sector.index()) {
    if (substr(si, 1, 2) %in% label_outside){
      xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
      ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
      circos.text(mean(xlim), max(ylim), substr(si, 1, 2), sector.index = si, track.index = 1, 
                  facing = "bending.inside", niceFacing = TRUE, col = "black", adj = c(0.45, -0.5))
    } else {
      xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
      ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
      circos.text(mean(xlim), mean(ylim), substr(si, 1, 2), sector.index = si, track.index = 1, 
                  facing = "bending.inside", niceFacing = TRUE)
    }
  }
  circos.clear()
}

plot_changes_pie <- function(changes_df, order=NA, labels=NA, colours=NA) {
  #' Plot specimen changes by type as a pie chart.
  #' 
  #' @param changes_df A data frame of specimen changes.
  #' @param order The order to plot the changes in.
  #' @param labels The labels to use for each order.
  #' @param colours A vector of colours for the segments.
  
  changes_count <- 
    changes_df %>%
    count(change)
  
  if (!all(is.na(order))) {
    changes_count <- changes_count %>%
      mutate(change=factor(change, levels=order, ordered=TRUE)) %>%
      complete(change, fill=list(n=0)) %>%
      arrange(rev(change))
  }
  
  changes_count <- 
    changes_count %>%
    mutate(breaks=cumsum(n) - n/2)
  
  amounts <- 
    changes_count %>%
    filter(n > 0) %>%
    pull(n)
  
  breaks <- cumsum(amounts) - amounts/2
  
  p <- 
    changes_count %>%
    ggplot(aes(x=1, y=n, fill=change)) +
    geom_bar(colour="black", stat="identity") +
    coord_polar(theta="y") +
    guides(fill=guide_legend(override_aes=list(colour=NA))) +
    theme_pubclean() + 
    theme(axis.ticks=element_blank(), 
          axis.title=element_blank(), 
          axis.text.y=element_blank(),  
          legend.position = "right") +
    scale_y_continuous(breaks=breaks, labels=amounts) + 
    labs(fill="")
  
  if (!is.na(labels)) {
    labels <-
      labels %>%
      paste("(") %>%
      paste(rev(amounts), sep="") %>%
      paste(")", sep="")
    
    p <- p + scale_fill_manual(labels = labels, values = colours)
  }
  
  return(p)
}

str_to_sentence <- function(string) {
  #' Converts strings to sentence case.
  #' 
  #' A convenience function to capitalise the first
  #' word in a string, and have everything else lower case.
  #' 
  #' @param string The string to convert to sentence case.
  #' 
  #' @return The sentence case string.
   
  string <- str_to_lower(string)
  substr(string, 1, 1) <- str_to_upper(substr(string, 1, 1))
  return (string)
}

plot_status_matrix <- function(status_changes, levels) {
  #' plot a heat map of conservation statuses across two data sets
  
  all_statuses <- status_changes %>% 
    count(eoo_2007, eoo_2017) %>% 
    mutate(eoo_2007 = factor(eoo_2007, levels=levels), 
           eoo_2017 = factor(eoo_2017, levels=levels)) %>% 
    complete(eoo_2007, eoo_2017)
  
  p <- ggplot(all_statuses, aes(x = eoo_2007, y = eoo_2017, fill = n))
  p <- p + geom_tile()
  
  p <- p + scale_fill_viridis(direction=1, name = "", na_value = "white") 
  p <- p + geom_text(aes(label=n), hjust=0.5, vjust=0.5, colour="white", na_rm = TRUE) 
  p <- p + theme_minimal() + theme(panel_grid.major = element_blank(), panel_grid.minor = element_blank()) 
  p <- p + labs(x="2007 Preliminary Assessment", y="2017 Preliminary Assessment")
  
  return(p)
}
