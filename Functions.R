# Functions of various purposes
# Written by Margaret Bolton mb804@exeter.ac.uk

# restrict (location.data, rastr)
## rasterizes location data from a data frame and returns TRUE if more than one cell is occupied
## location.data = dataframe with Longitude and Latitude
## rastr = empty raster of chosen resolution, e.g. raster(resolution = (60/60))

# range_distances (dat, range.pol)
## extracts maximum, minimum, and median latitude from range polygons
## extracts range area and range span (vertical distance between maximum and minimum latitudes)
## calculates vertical distance from sample location to maximum latitude

# gbif_plotter (synonym_row, dat, data_type)
## methods for plotting gbif_data

library(geosphere)
library(ggplot2)
library(ggspatial)
library(ggtext)
library(RColorBrewer)
library(raster)
library(tidyverse)

restrict <- function(location.data, rastr){
  c.rast <- raster::rasterize(location.data[[1]], rastr, fun = "count")
  return(length(Which(c.rast, cells = TRUE))>1) 
}


range_distances <- function(dat, range.pol){
  range.pol <- range.pol[range.pol@data$binomial %in% unique(dat$HostCorrectedName), ] # restrict range.pol to match hosts in dat
  
  dat.sets <- split(dat, f = dat$HostCorrectedName)
  
  dat.out <- lapply(dat.sets, function(k){
    range.sub <- range.pol[range.pol@data$binomial == unique(k$HostCorrectedName), ] # subset range.pol to host of interest
    range.abs.vals <- abs(raster::geom(range.sub)[, "y"]) # extract coordinates from range.sub and convert to absolute values
    
    range.max <- max(range.abs.vals)
    range.min <- min(range.abs.vals)
    range.span <- geosphere::distGeo(c(0, range.max), c(0, range.min))
    range.area <- sum(geosphere::areaPolygon(range.sub))
    range.median <- median(range.max, range.min)
    
    range.traits <- data.frame("HostCorrectedName" = unique(k$HostCorrectedName),
                               "RangeMax" = range.max,
                               "RangeMin" = range.min,
                               "RangeSpan" = range.span,
                               "RangeArea" = range.area,
                               "RangeMedian" = range.median)
    
    k$zeros <- 0 # for use as a longitude calculating distances from sample points in k
    
    k$EquatorwardsDist <- geosphere::distGeo(k[,c("zeros","Latitude")], c(0,range.min))
    k$EquatorwardsProp <- k$EquatorwardsDist/range.span
    k$MedianDist <- geosphere::distGeo(k[,c("zeros","Latitude")], c(0,range.median))
    k$MedianProp <- k$MedianDist/(range.span/2)
    
    k$zeros <- NULL
    
    return(list(k, range.traits))
  })
  
  out <- do.call(rbind, lapply(dat.out, function(dat) dat[[1]]))
  range.traits <- do.call(rbind, lapply(dat.out, function(dat) dat[[2]]))
  
  return(list("DistanceMetrics" = out, "RangeTraits" = range.traits))
}

gmpd_plotter <- function(host, dat, polys){
  species_dat <- filter(dat, HostCorrectedName == host)
  polys$legend <- as.factor(polys$legend)
  polys$id <- rownames(polys@data)
  species_IUCN <- polys[polys$binomial == host, ]
  species_IUCN <-  base::merge(species_IUCN@data[c("legend", "id")], fortify(species_IUCN), by = "id")
  
  # to make active levels bold in legend
  curr <- unique(species_IUCN$legend)
  new <- c(paste0("**",curr, "**"), levels(species_IUCN$legend)[!levels(species_IUCN$legend) %in% unique(species_IUCN$legend)])
  curr <- c(as.character(unique(species_IUCN$legend)), levels(species_IUCN$legend)[!levels(species_IUCN$legend) %in% unique(species_IUCN$legend)])
  species_IUCN$legend <- dplyr::recode(species_IUCN$legend, !!!deframe(data.frame(curr, new))) 
  
  # from https://sashamaps.net/docs/resources/20-colors/
  g.cols <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', 
              '#f58231', '#911eb4', '#42d4f4', '#f032e6', 
              '#bfef45', '#fabed4', '#469990', '#dcbeff', 
              '#9A6324', '#fffac8', '#800000', '#aaffc3', 
              '#808000', '#ffd8b1', '#000075', '#a9a9a9')
  
  ggplot() + coord_fixed() +
    borders("world", colour = "gray50", fill = "gray50") +
    
    geom_polypath(data = species_IUCN, 
                 aes(x = long, y = lat, group = group, 
                     colour = legend, fill = legend)) +
    
    scale_color_manual(values = g.cols, drop = FALSE) +
    scale_fill_manual(values = g.cols, drop = FALSE) +
    
    theme(legend.text = element_markdown(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.3, "lines"),
          legend.position = "bottom") +
    
    geom_point(data = species_dat,
               aes(x = Longitude, y = Latitude),
               colour = "navy", shape = 1, alpha = 0.5) +
    
    ggtitle(host)
  
}

gbif_plotter <- function(synonym_row, dat, data_type, polys){
  if (data_type == "base"){
    species_dat <- filter(dat, species == synonym_row["GBIFName"])
    
    # sorting out polygons to have data for legend
    polys$legend <- as.factor(polys$legend)
    polys$id <- rownames(polys@data)
    species_IUCN <- polys[polys$binomial == synonym_row["IUCNName"], ]
    species_IUCN <-  base::merge(species_IUCN@data[c("legend", "id")], fortify(species_IUCN), by = "id")
    
    # to make active levels bold in legend
    curr <- unique(species_IUCN$legend)
    new <- c(paste0("**",curr, "**"), levels(species_IUCN$legend)[!levels(species_IUCN$legend) %in% unique(species_IUCN$legend)])
    curr <- c(as.character(unique(species_IUCN$legend)), levels(species_IUCN$legend)[!levels(species_IUCN$legend) %in% unique(species_IUCN$legend)])
    species_IUCN$legend <- dplyr::recode(species_IUCN$legend, !!!deframe(data.frame(curr, new))) 
    
    # from https://sashamaps.net/docs/resources/20-colors/
    g.cols <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', 
                '#f58231', '#911eb4', '#42d4f4', '#f032e6', 
                '#bfef45', '#fabed4', '#469990', '#dcbeff', 
                '#9A6324', '#fffac8', '#800000', '#aaffc3', 
                '#808000', '#ffd8b1', '#000075', '#a9a9a9')
    
    ggplot() + coord_fixed() +
      borders("world", colour = "gray50", fill = "gray50") +
      
      geom_polypath(data = species_IUCN, 
                    aes(x = long, y = lat, group = group, 
                        colour = legend, fill = legend)) +
      
      scale_color_manual(values = g.cols, drop = FALSE) +
      scale_fill_manual(values = g.cols, drop = FALSE) +
      
      theme(legend.text = element_markdown(size = 6),
            legend.title = element_text(size = 6),
            legend.key.size = unit(0.3, "lines"),
            legend.position = "bottom") +
      
      geom_point(data = species_dat,
                 aes(x = decimalLongitude, y = decimalLatitude),
                 colour = "navy", shape = 1, alpha = 0.5) +
      
      ggtitle(paste0(synonym_row["GBIFName"], " (", synonym_row["IUCNName"], ")"))
    
  } else if (data_type == "bor"){
    species_dat <- filter(dat, species == synonym_row["GBIFName"])
    species_IUCN <- fortify(IUCN_Data_List[[synonym_row["IUCNName"]]])
    
    ggplot() + coord_fixed() +
      borders("world", colour = "gray50", fill = "gray50") +
      
      geom_polypath(data = species_IUCN, 
                   aes(x = long, y = lat, group = group),
                   colour = "skyblue", fill = "skyblue") +
      
      geom_point(data = species_dat,
                 aes(x = decimalLongitude, y = decimalLatitude, colour = basisOfRecord),
                 shape = 1) +
      scale_colour_brewer(palette = "Paired") +
      
      ggtitle(paste0(synonym_row["GBIFName"], " (", synonym_row["IUCNName"], ")"))
    
  } else if (data_type == "outlier"){
    species_dat <- dat[[synonym_row["GBIFName"]]] %>%
      mutate(outlier = case_when(is.na(outlier) ~ "untested",
                                 outlier        ~ "accepted",
                                 !outlier       ~ "rejected"))
    species_IUCN <- fortify(IUCN_Data_List[[synonym_row["IUCNName"]]])
    
    ggplot() + coord_fixed() +
      borders("world", colour = "gray50", fill = "gray50") +
      
      geom_point(data = species_dat,
                 aes(x = decimalLongitude, y = decimalLatitude, colour = outlier),
                 shape = 1, size = 1) +
      scale_colour_manual(values = c("firebrick", "chartreuse1")) +
      
      geom_polypath(data = species_IUCN, 
                   aes(x = long, y = lat, group = group),
                   colour = "skyblue", fill = NA) +
      
      ggtitle(paste0(synonym_row["GBIFName"], " (", synonym_row["IUCNName"], ")")) +
      
      if(!as.logical(synonym_row["cc_outl"]) & !as.logical(synonym_row["cc_iucn"])){
        labs(caption = "Coordinates untested")
      } else if (as.logical(synonym_row["cc_outl"])){
        labs(caption = paste0("Coordinates tested using 'outlier' method from 'cc_outl' with a multiple of ",
                              synonym_row["mltpl"]))
      } else if (as.logical(synonym_row["cc_iucn"])){
        labs(caption = paste0("Coordinates tested using iucn polygon with a buffer of ",
                              synonym_row["buffer"], " decimal degrees"))
      }
    
  }
}

species_cleaner <- function(synonym_row, dat){
  print(synonym_row["GBIFName"])
  species_dat <- filter(dat, species == synonym_row["GBIFName"])
  species_dat <- rename(species_dat, binomial = species)
  
  if(synonym_row["cc_outl"]){
    species_dat$outlier <- cc_outl(x = species_dat,
                                   lon = "decimalLongitude",
                                   lat = "decimalLatitude",
                                   method = "quantile",
                                   species = "binomial",
                                   mltpl = as.numeric(synonym_row["mltpl"]),
                                   value = "flagged")
    
  } else if(synonym_row["cc_iucn"]){
    species_IUCN <- IUCN_Data_List[[synonym_row["IUCNName"]]]
    species_IUCN$binomial <- synonym_row["GBIFName"]
    species_dat$outlier <- cc_iucn(x = species_dat,
                                   range = species_IUCN,
                                   lon = "decimalLongitude",
                                   lat = "decimalLatitude",
                                   species = "binomial",
                                   buffer = as.numeric(synonym_row["buffer"]),
                                   value = "flagged")
  } else {
    species_dat$outlier <- NA
  }
  
  species_dat <- rename(species_dat, species = binomial)
  return(species_dat)
}
