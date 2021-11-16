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
library(rnaturalearth)
library(sp)
library(tidyverse)

restrict <- function(location.data, rastr){
  c.rast <- raster::rasterize(location.data[[1]], rastr, fun = "count")
  return(length(Which(c.rast, cells = TRUE))>1) 
}


range_distances <- function(dat, range.pol, range.dat, method){
  if(method == "iucn"){
    range.pol <- range.pol[range.pol@data$binomial %in% unique(dat$HostCorrectedName), ] # restrict range.pol to match hosts in dat
    
    dat.sets <- split(dat, f = dat$HostCorrectedName)
    
    dat.out <- lapply(dat.sets, function(k){
      range.pol.sub <- range.pol[range.pol@data$binomial == unique(k$HostCorrectedName), ] # subset range.pol to current host 
      range.abs.vals <- abs(raster::geom(range.pol.sub)[, "y"]) # extract latitudes from range.pol.sub and convert to absolute values
      
      range.max <- max(range.abs.vals)
      range.min <- min(range.abs.vals)
      range.span <- geosphere::distGeo(c(0, range.max), c(0, range.min))
      range.area <- sum(geosphere::areaPolygon(range.pol.sub))
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
    
  } else if(method == "gbif"){
    range.pol <- range.pol[range.pol@data$binomial %in% unique(dat$HostCorrectedName), ] # restrict range.pol to match hosts in dat
    range.dat <- filter(range.dat, species %in% unique(dat$HostCorrectedName))
    
    dat.sets <- split(dat, f = dat$HostCorrectedName)
    
    dat.out <- lapply(dat.sets, function(k){
      range.pol.sub <- range.pol[range.pol@data$binomial == unique(k$HostCorrectedName), ] # subset range.pol to current host
      range.dat.sub <- filter(range.dat, species == unique(k$HostCorrectedName))
      range.abs.vals <- abs(range.dat.sub$decimalLatitude)
      
      range.max <- max(range.abs.vals)
      range.min <- min(range.abs.vals)
      range.span <- geosphere::distGeo(c(0, range.max), c(0, range.min))
      range.area <- sum(geosphere::areaPolygon(range.pol.sub))
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
    
  }
  
  return(list("DistanceMetrics" = out, "RangeTraits" = range.traits))
  
}

drop_introduced <- function(host, range.polygon, native.df){
  #keep <- filter(native.df, HostCorrectedName == host & keep)$Status
  throw <- filter(native.df, HostCorrectedName == host & !keep)$Status
  range.polygon <- range.polygon[range.polygon$binomial == host, ]
  
  if(length(throw) > 0){
    range.polygon <- range.polygon - range.polygon[range.polygon$legend %in% throw, ]
  }
  
  return(range.polygon)
}

gmpd_plotter <- function(host, dat, range.polygon, legend.text = Legend_Text, plot_type){
  if(plot_type == "iucn"){
    species_dat <- filter(dat, HostCorrectedName == host)
    range.polygon$legend <- factor(range.polygon$legend, levels = legend.text)
    range.polygon$id <- rownames(range.polygon@data)
    sp.range.polygon <- range.polygon[range.polygon$binomial == host, ]
    sp.range.polygon <-  base::merge(sp.range.polygon@data[c("legend", "id")], fortify(sp.range.polygon), by = "id")
    
    # to make active levels bold in legend
    curr <- unique(sp.range.polygon$legend)
    new <- c(paste0("**",curr, "**"), levels(sp.range.polygon$legend)[!levels(sp.range.polygon$legend) %in% curr])
    curr <- c(as.character(curr), levels(sp.range.polygon$legend)[!levels(sp.range.polygon$legend) %in% curr])
    sp.range.polygon$legend <- dplyr::recode(sp.range.polygon$legend, !!!deframe(data.frame(curr, new))) 
    
    # from https://sashamaps.net/docs/resources/20-colors/
    g.cols <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', 
                '#f58231', '#911eb4', '#42d4f4', '#f032e6', 
                '#bfef45', '#fabed4', '#469990', '#dcbeff', 
                '#9A6324', '#fffac8', '#800000', '#aaffc3', 
                '#808000', '#ffd8b1', '#000075', '#a9a9a9')
    
    # creating bounding box to clip plots
    # Blocked out for now because coord_sf is buggy
    # bbox <- range.polygon@bbox
    # dat.bbox <- SpatialPoints(species_dat[c("Longitude", "Latitude")])@bbox
    # bbox <- matrix(c(min(c(bbox[1,1], dat.bbox[1,1])),
    #                   min(c(bbox[2,1], dat.bbox[2,1])),
    #                   max(c(bbox[1,2], dat.bbox[1,2])),
    #                   max(c(bbox[2,2], dat.bbox[2,2]))),
    #                   c(2,2)) # there must be a one-liner for this..
    # 
    # # expanding bboxes that are too small to give context
    # if(bbox[1,2] - bbox[1,1] < 10){
    #   bbox[1,1] <- bbox[1,1] -5
    #   bbox[1,2] <- bbox[1,2] +5
    # }
    # if(bbox[2,2] - bbox[2,1] < 10){
    #   bbox[2,1] <- bbox[2,1] -5
    #   bbox[2,2] <- bbox[2,2] +5
    # }
    
    ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
      # coord_sf(xlim = bbox[1, ], ylim = bbox[2, ], expand = TRUE) +
      xlab("Longitude") + ylab("Latitude") +
      
      geom_polypath(data = sp.range.polygon, 
                    aes(x = long, y = lat, group = group, 
                        colour = legend, fill = legend),
                    size = 0.15) +
      
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
  } else if(plot_type == "base") {
    species_dat <- filter(dat, HostCorrectedName == host)
    
    ## Currently there's a bug from the recent crs change
    ## https://rgdal.r-forge.r-project.org/articles/CRS_projections_transformations.html
    # bbox <- SpatialPoints(species_dat[c("Longitude", "Latitude")])@bbox
    # 
    # #expanding smaller bboxes to give context
    # if(bbox[1,2] - bbox[1,1] < 15){
    #   bbox[1,1] <- bbox[1,1] - round((15 - (bbox[1,2] - bbox[1,1]))/2, 2)
    #   bbox[1,2] <- bbox[1,2] + round((15 - (bbox[1,2] - bbox[1,1]))/2, 2)
    # }
    # if(bbox[2,2] - bbox[2,1] < 15){
    #   bbox[2,1] <- bbox[2,1] - round((15 - (bbox[2,2] - bbox[2,1]))/2, 2)
    #   bbox[2,2] <- bbox[2,2] + round((15 - (bbox[2,2] - bbox[2,1]))/2, 2)
    # }
    
    ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
      # coord_sf(xlim = bbox[1, ], ylim = bbox[2, ], expand = TRUE) +
      xlab("Longitude") + ylab("Latitude") +
      
      geom_point(data = species_dat,
                 aes(x = Longitude, y = Latitude),
                 colour = "navy", shape = 1, alpha = 0.5) +
      
      ggtitle(host)
  }
}

iucn_test <- function(host, dat = GMPD_Data, range.polygon = IUCN_Data_List, native.df = Native_Clean){
  species_dat <- dat %>%
    filter(HostCorrectedName == host) 
  native.df <- filter(native.df, HostCorrectedName == host)
  
  throw <- native.df$Status
  range.polygon <- range.polygon[[host]]
  
  if(length(throw) > 0){
    removed.polygon <- range.polygon[range.polygon$legend %in% throw, ]
  }
  
  species_dat$out <-  cc_iucn(x = rename(species_dat, binomial = HostCorrectedName),
                              range = removed.polygon,
                              lon = "Longitude",
                              lat = "Latitude",
                              species = "binomial",
                              buffer = unique(native.df$buffer),
                              value = "flagged")
  species_dat <- species_dat %>% 
    mutate(out = case_when(out  ~ "Discarded",
                           !out ~ "Kept")) 
  species_dat$out <- factor(species_dat$out, levels = c("Kept", "Discarded"))
  
  bbox <- range.polygon@bbox
  dat.bbox <- SpatialPoints(species_dat[c("Longitude", "Latitude")])@bbox
  bbox <- matrix(c(min(c(bbox[1,1], dat.bbox[1,1])),
                   min(c(bbox[2,1], dat.bbox[2,1])),
                   max(c(bbox[1,2], dat.bbox[1,2])),
                   max(c(bbox[2,2], dat.bbox[2,2]))),
                 c(2,2)) # there must be a one-liner for this..
  
  ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
    geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
    coord_sf(xlim = bbox[1, ], ylim = bbox[2, ], expand = TRUE) +
    xlab("Longitude") + ylab("Latitude") +
    
    geom_polypath(data = fortify(range.polygon), 
                  aes(x = long, y = lat, group = group),
                  colour = "palegreen3",
                  fill = "palegreen3") +
    
    geom_polypath(data = fortify(gBuffer(removed.polygon, byid = FALSE, width = 0.1)), 
                  aes(x = long, y = lat, group = group),
                  colour = NA,
                  fill = "firebrick",
                  alpha = 0.5) +
    
    geom_polypath(data = fortify(removed.polygon), 
                  aes(x = long, y = lat, group = group),
                  colour = "firebrick",
                  fill = "firebrick") +
    
    geom_point(data = species_dat,
               aes(x = Longitude, y = Latitude, colour = out)) +
    scale_colour_discrete(drop = FALSE) +
    
    theme(legend.position = "bottom") +
    ggtitle(host)
}

iucn_cleaning <- function(dat, native.df, range.polygon = IUCN_Data){
  hosts <- unique(native.df$HostCorrectedName)
  for(host in hosts){
    sp.range.polygon <- range.polygon[range.polygon$binomial == host, ]
    species.dat <- split(dat, dat$HostCorrectedName == host)
    sp.native.df <- filter(native.df, HostCorrectedName == host)
    
    throw <- sp.native.df$Status
    removed.polygon <- sp.range.polygon[sp.range.polygon$legend %in% throw, ]
    
    species.dat[["TRUE"]]$out <-  cc_iucn(x = rename(species.dat[["TRUE"]], binomial = HostCorrectedName),
                                range = removed.polygon,
                                lon = "Longitude",
                                lat = "Latitude",
                                species = "binomial",
                                buffer = unique(sp.native.df$buffer),
                                value = "flagged")
    
    species.dat[["TRUE"]] %>%
      filter(!out) %>%
      dplyr::select(-out)
    
    dat <- bind_rows(species.dat)
  }
  return(dat)
}

gbif_plotter <- function(synonym_row, dat, data_type, range.polygon, legend.text = Legend_Text){
  if (data_type == "base"){
    species_dat <- filter(dat, species == synonym_row["GBIFName"])
    
    # sorting out polygons to have data for legend
    range.polygon$legend <- factor(range.polygon$legend, levels = legend.text)
    range.polygon$id <- rownames(range.polygon@data)
    sp.range.polygon <- range.polygon[range.polygon$binomial == synonym_row["IUCNName"], ]
    sp.range.polygon <-  base::merge(sp.range.polygon@data[c("legend", "id")], fortify(sp.range.polygon), by = "id")
    
    # to make active levels bold in legend
    curr <- unique(sp.range.polygon$legend)
    new <- c(paste0("**", curr, "**"), levels(sp.range.polygon$legend)[!levels(sp.range.polygon$legend) %in% curr])
    curr <- c(as.character(curr), levels(sp.range.polygon$legend)[!levels(sp.range.polygon$legend) %in% curr])
    sp.range.polygon$legend <- dplyr::recode(sp.range.polygon$legend, !!!deframe(data.frame(curr, new))) 
    
    # from https://sashamaps.net/docs/resources/20-colors/
    g.cols <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', 
                '#f58231', '#911eb4', '#42d4f4', '#f032e6', 
                '#bfef45', '#fabed4', '#469990', '#dcbeff', 
                '#9A6324', '#fffac8', '#800000', '#aaffc3', 
                '#808000', '#ffd8b1', '#000075', '#a9a9a9')
    
    # bbox <- range.polygon@bbox
    # dat.bbox <- SpatialPoints(species_dat[c("Longitude", "Latitude")])@bbox
    # bbox <- matrix(c(min(c(bbox[1,1], dat.bbox[1,1])),
    #                  min(c(bbox[2,1], dat.bbox[2,1])),
    #                  max(c(bbox[1,2], dat.bbox[1,2])),
    #                  max(c(bbox[2,2], dat.bbox[2,2]))),
    #                c(2,2)) # there must be a one-liner for this..
    
    ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
      # coord_sf(xlim = bbox[1, ], ylim = bbox[2, ], expand = TRUE) +
      xlab("Longitude") + ylab("Latitude") +
      
      geom_polypath(data = sp.range.polygon, 
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
    dat$basisOfRecord <- as.factor(dat$basisOfRecord)
    species_dat <- filter(dat, species == synonym_row["GBIFName"])
    sp.range.polygon <- range.polygon[range.polygon$binomial == synonym_row["IUCNName"], ]
    
    # to make active levels bold in legend
    curr <- unique(species_dat$basisOfRecord)
    new <- c(paste0("**", curr, "**"), levels(species_dat$basisOfRecord)[!levels(species_dat$basisOfRecord) %in% curr])
    curr <- c(as.character(curr), levels(species_dat$basisOfRecord)[!levels(species_dat$basisOfRecord) %in% curr])
    species_dat$basisOfRecord <- dplyr::recode(species_dat$basisOfRecord, !!!deframe(data.frame(curr, new))) 
    
    ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
      
      geom_polypath(data = sp.range.polygon, 
                   aes(x = long, y = lat, group = group),
                   colour = "skyblue", fill = "skyblue") +
      
      geom_point(data = species_dat,
                 aes(x = decimalLongitude, y = decimalLatitude, colour = basisOfRecord),
                 shape = 1) +
      scale_colour_brewer(palette = "Paired", drop = FALSE) +
      
      theme(legend.text = element_markdown()) +
      
      ggtitle(paste0(synonym_row["GBIFName"], " (", synonym_row["IUCNName"], ")"))
    
  } else if (data_type == "outlier"){
    species_dat <- dat[[synonym_row["GBIFName"]]] %>%
      mutate(outlier = case_when(is.na(outlier) ~ "untested",
                                 outlier        ~ "accepted",
                                 !outlier       ~ "rejected"))
    sp.range.polygon <- range.polygon[range.polygon$binomial == synonym_row["IUCNName"], ]
    
    ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
      
      geom_point(data = species_dat,
                 aes(x = decimalLongitude, y = decimalLatitude, colour = outlier),
                 shape = 1, size = 1) +
      scale_colour_manual(values = c("firebrick", "chartreuse1")) +
      
      geom_polypath(data = sp.range.polygon, 
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
    sp.range.polygon <- IUCN_Data_List[[synonym_row["IUCNName"]]]
    sp.range.polygon$binomial <- synonym_row["GBIFName"]
    species_dat$outlier <- cc_iucn(x = species_dat,
                                   range = sp.range.polygon,
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
