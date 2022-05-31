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
library(tmap)

restrict_grid <- function(location.data, rastr){
  c.rast <- raster::rasterize(location.data[[1]], rastr, fun = "count")
  return(length(Which(c.rast, cells = TRUE))>1) 
}

restrict_deci <- function(dat, subsp = FALSE){
  hostlist <- unique(dat$HostCorrectedName)

  enough <- lapply(hostlist, function(host){
    sp.dat <- dat[dat$HostCorrectedName == host,]
    sp.dat <- sp.dat %>%
      mutate(across(c(Latitude, Longitude), ~ trunc(.x)))
    
    if(subsp){
      subgroups <- unique(sp.dat$subgroup)
      nrow.sg.dat <- lapply(subgroups, function(sg){
        sg.dat <- sp.dat[sp.dat$subgroup == sg,]
        
        out <- data.frame(HostCorrectedName = host,
                          subgroup = sg, 
                          enough = nrow(unique(sg.dat[c("Latitude", "Longitude")])) > 1)
      })
      nrow.sg.dat <- bind_rows(nrow.sg.dat)
      
    } else {
      out <- data.frame(HostCorrectedName = host,
                        enough = nrow(unique(sp.dat[c("Latitude", "Longitude")])) > 1)
    }
  })
  enough <- bind_rows(enough)
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
    range.dat <- filter(range.dat, species %in% unique(dat$HostCorrectedName)) # restrict range.dat to match hosts in dat
    
    dat.sets <- split(dat, f = dat$HostCorrectedName) # split input data by host
    
    dat.out <- lapply(dat.sets, function(k){
      range.dat.sub <- range.dat[range.dat$species == unique(k$HostCorrectedName), ] # subset range.dat to current host
      range.abs.vals <- abs(range.dat.sub@coords[,"Latitude"]) # vector of absolute latitudes in gbif
      
      range.max <- max(range.abs.vals)
      range.min <- min(range.abs.vals)
      range.span <- geosphere::distGeo(c(0, range.max), c(0, range.min))
      # range.area <- ## Currently no method for range area from gbif, but could use convex hull
      range.median <- median(range.max, range.min)
      range.traits <- data.frame("HostCorrectedName" = unique(k$HostCorrectedName),
                                 "RangeMax" = range.max,
                                 "RangeMin" = range.min,
                                 "RangeSpan" = range.span,
                                 #"RangeArea" = range.area,
                                 "RangeMedian" = range.median)
      
      k$zeros <- 0 # for use as a longitude calculating distances from sample points in k
      
      k$EquatorwardsDist <- geosphere::distGeo(k[,c("zeros","Latitude")], c(0,range.min)) # vertical distances to lowest latitude from input data sample locations
      k$EquatorwardsProp <- k$EquatorwardsDist/range.span # as a proportion of range span
      k$MedianDist <- geosphere::distGeo(k[,c("zeros","Latitude")], c(0,range.median)) # vertical distances to range median
      k$MedianProp <- k$MedianDist/(range.span/2) # as proportion of half range span
      
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
    sp.range.polygon <- base::merge(sp.range.polygon@data[c("legend", "id")], fortify(sp.range.polygon), by = "id")
    
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
    
  } else {
    dat[,data_type] <- as.factor(dat[, data_type])
    species_dat <- filter(dat, species == synonym_row["GBIFName"])
    sp.range.polygon <- range.polygon[range.polygon$binomial == synonym_row["IUCNName"], ]
    
    # to make active levels bold in legend
    curr <- unique(species_dat[, data_type])
    new <- c(paste0("**", curr, "**"), levels(species_dat[, data_type])[!levels(species_dat[, data_type]) %in% curr])
    curr <- c(as.character(curr), levels(species_dat[, data_type])[!levels(species_dat[, data_type]) %in% curr])
    species_dat[, data_type] <- dplyr::recode(species_dat[, data_type], !!!deframe(data.frame(curr, new))) 
    
    ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
      
      geom_polypath(data = sp.range.polygon, 
                    aes(x = long, y = lat, group = group),
                    colour = "skyblue", fill = "skyblue") +
      
      geom_point(data = species_dat,
                 aes(x = decimalLongitude, y = decimalLatitude, colour = data_type),
                 shape = 1) +
      scale_colour_brewer(palette = "Paired", drop = FALSE) +
      
      theme(legend.text = element_markdown()) +
      
      ggtitle(paste0(synonym_row["GBIFName"], " (", synonym_row["IUCNName"], ")"))
    
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

# Not sure whether to keep
quick_map <- function(c.species){
  sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == c.species, ]
  sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == c.species, c("Longitude", "Latitude")])
  tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 
}

pip_test <- function(host, dat, range.polygon, buff, subsp = TRUE){
  if("HostCorrectedName" %in% names(dat)){
    sp.dat <- dat[dat$HostCorrectedName == host,]
  } else if("species" %in% names(dat)){
    sp.dat <- dat[dat$species == host,] 
  } 
  
  if(subsp){
    sg.dat <- lapply(unique(buff$subgroup), function(sg){
      sg.polygon <- range.polygon[range.polygon$subgroup == sg,] 
      if(buff[buff$subgroup == sg, "buff"] != 0){
        sg.polygon <- terra::buffer(sg.polygon, buff[buff$subgroup == sg, "buff"]) 
      }
      
      sg.dat <- sp.dat[sg.polygon,] 
      if(nrow(sg.dat) > 0){
        sg.dat@data$subgroup <- sg
      }
      
      return(sg.dat)
    })
    
    sg.dat <- sg.dat[which(lapply(sg.dat, function(x)nrow(x) != 0) == TRUE)] 
    if(length(sg.dat) == 0){
      return(NULL)
    } else if(length(sg.dat) == 1){
      return(sg.dat[[1]])
    } else {
      suppressWarnings(sg.dat <- do.call(raster::bind, sg.dat))
      
      if(nrow(sg.dat) > nrow(sp.dat)){
        warning("point(s) assigned multiple subgroups")
      } else if(nrow(sg.dat) < nrow(sp.dat)){
        warning("not all points assigned to subgroup")
      }
      
      return(sg.dat)
      
    }

  } else{
    
    if(buff != 0){
      range.polygon <- terra::buffer(range.polygon, buff)
    }
    
    sp.dat <- sp.dat[range.polygon,]
    return(sp.dat)
  }
}

pip_test_host <- function(hostlist, dat, range.polygon, buff, subsp = TRUE){
  out <- lapply(hostlist, function(host){
    range.polygon <- range.polygon[range.polygon$binomial == host,]
    out <- pip.test(host, dat, range.polygon, buff, subsp)
  })
  out <- do.call(raster::bind, sg.dat)
}

complete_plot <- function(synonym.row, dat, range.dat, range.polygon){
  # data prep
  sp.dat <- dat[dat$HostCorrectedName == synonym.row["IUCNName"],]
  sp.range.dat <- range.dat[range.dat$species == synonym.row["GBIFName"],]
  sp.range.pol <- range.polygon[range.polygon$binomial == synonym.row["IUCNName"],]
  
  sp.bbox <- bboxes[[synonym.row["scaling"]]]
  subgroups <- unique(sp.range.pol@data$subgroup)
  
  extant.group <- poly_fill[poly_fill$fill_group == "Extant", "legend"]
  present.group <- poly_fill[poly_fill$fill_group == "Presence Likely", "legend"]
  absent.group <- poly_fill[poly_fill$fill_group == "Absence Likely", "legend"]
  extinct.group <- poly_fill[poly_fill$fill_group == "Extinct", "legend"]
  
  border.factor <- 1
  if(synonym.row["scaling"] %in% c("Africa", "South", "North", "Central", "Asia", "Europe")) {border.factor <- 0.5}
  
  # plot 1, iucn range
  par(bg = 'powderblue')
  terra::plot(ne_countries(scale = "medium", returnclass = "sp"), col = "beige", border = "burlywood",
       xlim = sp.bbox[1,],
       ylim = sp.bbox[2,])
  
  for(i in 1:length(subgroups)){
    sp.range.pol.sub <- sp.range.pol[sp.range.pol$subgroup == subgroups[i],]
    terra::plot(sp.range.pol.sub, col = sub.colours[i], border = "transparent", add = TRUE)
    
    temp.border <- spTransform(sp.range.pol.sub, CRS("+init=epsg:3857"))
    temp.border <- try(terra::buffer(terra::buffer(temp.border, 10000), -100000*border.factor), silent = TRUE)

    if(class(temp.border) != "try-error") {
      temp.border <- spTransform(temp.border, CRS("+init=epsg:4326"))
      
      if(area(temp.border) > 100000){
        
        if(any(sp.range.pol.sub$legend %in% c(present.group, absent.group, extinct.group))){
          terra::plot(raster::intersect(sp.range.pol.sub[sp.range.pol.sub$legend %in% c(present.group, absent.group, extinct.group),], temp.border), 
               col = "grey65", border = "transparent", add = TRUE)
        }
        
        if(any(sp.range.pol.sub$legend %in% c(absent.group, extinct.group))){
          terra::plot(raster::intersect(sp.range.pol.sub[sp.range.pol.sub$legend %in% c(absent.group, extinct.group),], temp.border), 
               col = "black", density = 15, angle = 45, border = "transparent", add = TRUE)
        }
        
        if(any(sp.range.pol.sub$legend %in% extinct.group)){
          terra:plot(raster::intersect(sp.range.pol.sub[sp.range.pol.sub$legend %in% extinct.group,], temp.border), 
               col = "black", density = 15, angle = 135, border = "transparent", add = TRUE)
        }
        
        if(any(sp.range.pol.sub$legend %in% c(extant.group))){
          terra::plot(raster::intersect(sp.range.pol.sub[sp.range.pol.sub$legend %in% extant.group,], temp.border), 
               col = "grey", border = "transparent", add = TRUE)
        }
        
      }
      
    }
    
    terra::plot(sp.dat[sp.dat$subgroup == subgroups[i] & sp.dat$RestrAll & sp.dat$RestrSub,], 
                pch = 24, bg = sub.colours[i], col = "grey30", lwd = 2, add = TRUE)

    terra::plot(sp.dat[sp.dat$subgroup == subgroups[i] & sp.dat$RestrAll & !sp.dat$RestrSub,], 
                pch = 25, bg = sub.colours[i], col = "grey30", lwd = 2, add = TRUE)
    
  }
  
  mtext(paste0(synonym.row["IUCNName"], " with iucn range\n"), side = 3, cex = 1.5)
  
  legend("bottomleft", 
         legend = subgroups, 
         title = "Subgroup",
         fill = sub.colours,
         border = sub.colours,
         cex = 1,
         bg = "transparent") 
  
  legend("topleft", 
         legend = c("", "", "", ""), 
         title = "",
         fill = c("grey", "grey65", "grey65", "grey65"),
         cex = 1, 
         bty = "n") 
  
  legend("topleft", 
         legend = unique(poly_fill$fill_group), 
         title = "Status",
         fill = c("transparent", "transparent", "black", "black"),
         density = c(0, 0, 15, 15),
         angle = c(0, 0, 45, 135),
         cex = 1,
         bty = "o",
         bg = "transparent") 
  
  legend("bottomright",
         legend = c("All and subgrouped", "All only"),
         title = "Restricted Datasets",
         pch = c(24, 25),
         cex = 1,
         pt.bg = sub.colours[1],
         pt.lwd = 2,
         bg = "transparent")
  
  # plot 2, gbif range
  par(bg = 'powderblue')
  terra::plot(ne_countries(scale = "medium", returnclass = "sp"), col = "beige", border = "burlywood",
              xlim = sp.bbox[1,],
              ylim = sp.bbox[2,])
  
  for(i in 1:length(subgroups)){
    terra::plot(sp.range.dat[sp.range.dat$subgroup == subgroups[i],], 
                pch = 1, col = paste0(sub.colours[i], "FF"), add = TRUE)
    
    terra::plot(sp.dat[sp.dat$subgroup == subgroups[i] & sp.dat$CleanAll & sp.dat$CleanSub,], 
                pch = 22, bg = sub.colours[i], col = "grey30", lwd = 2, add = TRUE)

    terra::plot(sp.dat[sp.dat$subgroup == subgroups[i] & sp.dat$CleanAll & !sp.dat$CleanSub,], 
                pch = 23, bg = sub.colours[i], col = "grey30", lwd = 2, add = TRUE)
    
  }
  
  mtext(paste0(synonym.row["IUCNName"], " with gbif range\n"), side = 3, cex = 1.5)
  
  legend("bottomleft", 
         legend = subgroups, 
         title = "Subgroup",
         fill = sub.colours,
         border = sub.colours,
         cex = 1,
         bg = "transparent") 

  legend("bottomright",
         legend = c("All and subgrouped", "All only"),
         title = "Cleaned Datasets",
         pch = c(22, 23),
         cex = 1,
         pt.bg = sub.colours[1],
         pt.lwd = 2,
         bg = "transparent")
  
}

bboxer <- function(...){
  x <- list(...)
  
  matrix(c(min(sapply(x, function(x) x[1,1])), max(sapply(x, function(x) x[1,2])),
           min(sapply(x, function(x) x[2,1])), max(sapply(x, function(x) x[2,2]))),
         byrow = TRUE, nrow = 2,
         dimnames = list(c("Longitude", "Latitude"), c("min", "max")))
}
