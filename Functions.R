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


gbif_plotter <- function(synonym_row){
  species_dat <- filter(GBIF_Data, species == synonym_row["GBIFName"])
  species_IUCN <- fortify(IUCN_Data_List[[synonym_row["IUCNName"]]])
  
  ggplot() + coord_fixed() +
    borders("world", colour = "gray50", fill = "gray50") +
    
    geom_polygon(data = species_IUCN, 
                 aes(x = long, y = lat, group = group),
                 colour = "black",
                 fill = NA) +
    
    geom_point(data = species_dat,
               aes(x = decimalLongitude, y = decimalLatitude),
               colour = "blue") +
    
    ggtitle(paste0(synonym_row["GBIFName"], " (", synonym_row["IUCNName"], ")"))
}

