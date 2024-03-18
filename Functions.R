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
# library(raster) # TODO: use terra instead
# library(rnaturalearth) # TODO: depends on sp
# library(sp)
library(tidyverse)
library(tmap)

# TODO: sort functions
# used in 01 ############################################################
pip_test <- function(point.data, range.polygon, buff = 0){
  # TODO: get rid of this silly function
  if(buff != 0){
    range.polygon <- sf::st_buffer(range.polygon, buff)
  }
  
  point.data <- point.data[range.polygon,]
  return(point.data)
}


plot_native <- function(hostname, buff = 0){
  if(!exists("World")) data("World")
  
  species_polygon <- IUCN_Data[IUCN_Data$sci_name == hostname, ]
  species_dots <- GMPD_Data[GMPD_Data$HostCorrectedName == hostname, ]
  
  bbox_polygon <- sf::st_as_sfc(sf::st_bbox(species_polygon))
  bbox_dots <- sf::st_as_sfc(sf::st_bbox(species_dots))
  bbox_plot <- sf::st_bbox(sf::st_union(bbox_polygon, bbox_dots))

  if(buff == 0){
    tm_shape(World) + tm_fill() +
      tm_shape(species_polygon) +
      tm_polygons("keep") +
      tm_shape(species_dots) +
      tm_dots() +
      tm_layout(title = hostname)
  } else {
    tm_shape(World) + tm_fill() +
      tm_shape(species_polygon) +
      tm_polygons("keep") +
      tm_shape(sf::st_buffer(species_polygon, buff)) +
      tm_polygons("keep", alpha = 0.3) +
      tm_shape(species_dots) +
      tm_dots() +
      tm_layout(title = hostname)
  }
}

restrict_decimal <- function(dat, subsp = FALSE){
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

# used in 02 ###################################################################
plot_native_gbif <- function(hostname, buff = 0){
  if(!exists("World")) data("World")
  
  species_polygon <- IUCN_Data[IUCN_Data$sci_name == hostname, ]
  species_dots <- GBIF_Data[GBIF_Data$species == hostname, ]
  
  bbox_polygon <- sf::st_as_sfc(sf::st_bbox(species_polygon))
  bbox_dots <- sf::st_as_sfc(sf::st_bbox(species_dots))
  bbox_plot <- sf::st_bbox(sf::st_union(bbox_polygon, bbox_dots))
  
  if(buff == 0){
    
    plot_out <- tm_shape(World, bbox = bbox_plot) + tm_fill() +
      tm_shape(species_dots) +
      tm_dots() +
      tm_shape(species_polygon) +
      tm_fill("keep", alpha = 0.4) +
      tm_layout(title = hostname,)
    
  } else {
    
    plot_out <- tm_shape(World, bbox = bbox_plot) + tm_fill() +
      tm_shape(species_polygon) +
      tm_fill("keep") +
      tm_shape(species_dots) +
      tm_dots() +
      tm_layout(title = hostname) +
      tm_shape(sf::st_buffer(species_polygon, buff)) +
      tm_fill("keep", alpha = 0.3, legend.show = FALSE)
    
  }
  return(plot_out)
}


# used in 03 ###################################################################

range_distances_host <- function(host, dat, range.object, method){
  
  dat <- dat %>% filter(HostCorrectedName == host)
  
  # subset range.object to current host 
  if("sci_name" %in% names(range.object)){
  range.object <- range.object %>% filter(sci_name == host)
  
  } else if("species" %in% names(range.object)){
  range.object <- range.object %>% filter(species == host)
  
  }
  
  if(!"sf" %in% class(range.object)){
    range.object <- sf::st_as_sf(range.object,
                                 coords = c("decimalLongitude", "decimalLatitude"),
                                 crs = Projection_String)
  }

  out <- range_distances_method(dat = dat, range.object = range.object)
  
  out[[2]]$HostCorrectedName <- host
  out[[2]]$RangeArea <- sum(sf::st_area(IUCN_Native_Data[IUCN_Native_Data$sci_name == host,]))
  
  out[[1]]$RangeMethod <- method
  out[[2]]$RangeMethod <- method
  out[[1]]$RangeTaxonLvl <- "species"
  out[[2]]$RangeTaxonLvl <- "species"
  
  return(out)
}

# range_distances_subsp <- function(host, dat, range.pol, method){
#   
#   dat <- dat %>% filter(HostCorrectedName == host)
#   subsp.list <- dat %>% pull(subgroup) %>% unique()
#   
#   if(method == "iucn"){
#     
#     # subset range.pol to current host 
#     range.pol <- range.pol %>% filter(sci_name == host)
#     
#     out <- lapply(subsp.list, function(subspecies){
#       
#       # subset range.pol and data to current subgroup
#       dat <- dat %>% filter(subgroup == subspecies)
#       range.pol.sub <- range.pol %>% filter(sci_name == subspecies)
#       
#       out.sub <- iucn_method(dat = dat, range.polygon = range.pol.sub)
#       out.sub$range.traits <- cbind(out$range.traits, "subgroup" = subspecies)
#       
#       return(out.sub)
#     })
#     
#     out <- do.call(rbind, out)
#     
#   } else if(method == "gbif"){
#     
#     # subset range.dat to current host 
#     range.dat <- range.dat %>% filter(sci_name == host)
#     
#     out <- lapply(subsp.list, function(subspecies){
#       
#       # subset range.dat and data to current subgroup
#       dat <- dat %>% filter(subgroup == subspecies)
#       range.dat.sub <- range.dat %>% filter(sci_name == subspecies)
#       
#       # drop hosts with insufficient data
#       if(nrow(range.dat.sub) == 0){
#         
#         out.sub <- NULL
#         writeLines(paste(host, subspecies, "dropped, no data", sep = " "))
#         
#       } else if(nrow(range.dat.sub) < 2){
#         
#         out.sub <- NULL
#         writeLines(paste(host, subspecies, "dropped, insufficient data", sep = " "))
#         
#       } else {
#         
#         out.sub <- gbif_method(dat = dat, range.points = range.dat)
#         out.sub$range.traits <- cbind(out$range.traits, "subgroup" = subspecies)
#       
#       }
#       
#       return(out.sub)
#       
#     })
#     
#     out <- do.call(rbind, out)
#     
#   }
#   return(list("DistanceMetrics" = out, "RangeTraits" = range.traits))
#   
# }

range_distances_method <- function(dat, range.object){
  
  # extract latitudes from range.object and convert to absolute values
  range.vals.exact  <- sf::st_coordinates(range.object)[,"Y"]
  range.vals.abs    <- abs(range.vals.exact) 

  # absolute values (latitude)
  range.max.abs         <- max(range.vals.abs)
  range.max.abs.point   <- c(0, range.max.abs) %>% sf::st_point() %>% sf::st_sfc(crs = Projection_String) 
  
  range.min.abs         <- min(range.vals.abs)
  range.min.abs.point   <- c(0, range.min.abs) %>% sf::st_point() %>% sf::st_sfc(crs = Projection_String)
  
  # using the absolute minimum rather than zero because there may be some which span the equator but do not cross it
  # only really matters when it comes to error in analysis but this is avoided because we're now using median instead of span
  range.span.abs        <- sf::st_distance(range.max.abs.point, range.min.abs.point) %>% as.vector()
  
  # exact values (median)
  range.max.exact       <- max(range.vals.exact)
  range.max.exact.point <- c(0, range.max.exact) %>% sf::st_point() %>% sf::st_sfc() %>% sf::st_set_crs(Projection_String)
  
  range.min.exact       <- min(range.vals.exact)
  range.min.exact.point <- c(0, range.min.exact) %>% sf::st_point() %>% sf::st_sfc() %>% sf::st_set_crs(Projection_String)
  
  range.span.exact      <- sf::st_distance(range.max.exact.point, range.min.exact.point) %>% as.vector()
  
  range.median          <- median(c(range.max.exact, range.min.exact))
  range.median.point    <- c(0, range.median) %>% sf::st_point() %>% sf::st_sfc(crs = Projection_String)
  
  # range traits, distances, proportions
  range.traits <- data.frame("RangeMaxAbs"   = range.max.abs,
                             "RangeMaxExact" = range.max.exact,
                             
                             "RangeMinAbs"   = range.min.abs,
                             "RangeMinExact" = range.min.exact,
                             
                             "RangeSpanAbs"   = range.span.abs,
                             "RangeSpanExact" = range.span.exact,
                             
                             "RangeMedian" = range.median)
  
  # use absolute and exact values to calculate metrics
  
  # drop samples outside range limits
  dat.out <- dat[!(dat$Latitude > range.max.exact | dat$Latitude < range.min.exact),]
  
  if(nrow(dat.out) == 0){
    dat.out <- rbind(NA, dat.out) 
    writeLines("All points outside range margins")
    
  } else {
    
    # use absolute and exact values to calculate metrics
    
    dat.out$zeros <- 0 # for use as a longitude calculating distances from sample points in dat.out
    dat.out$absLat <- abs(dat.out$Latitude)
    
    dat.out.spatial.abs   <- sf::st_as_sf(dat.out,
                                      coords = c("zeros", "absLat"),
                                      crs = Projection_String)
    dat.out.spatial.exact <- sf::st_as_sf(dat.out,
                                      coords = c("zeros", "Latitude"),
                                      crs = Projection_String)
    
    dat.out$EquatorwardsDist <- sf::st_distance(dat.out.spatial.abs, range.min.abs.point) %>% as.vector() # vertical distances to lowest abs latitude from sample locations
    dat.out$EquatorwardsProp <- dat.out$EquatorwardsDist/range.span.abs # as a proportion of abs range span
    
    dat.out$MedianDist       <- sf::st_distance(dat.out.spatial.exact, range.median.point) %>% as.vector() # vertical distances to exact range median
    dat.out$MedianProp       <- dat.out$MedianDist/(range.span.exact/2) # as proportion of half range span
    dat.out$AboveMedn        <- (dat.out$Latitude > range.median & range.median > 0) | (dat.out$Latitude < range.median & range.median < 0) # "above" median really means polewards from median
    
    dat.out$zeros <- NULL
    dat.out$absLat <- NULL
    
    writeLines(paste(nrow(dat) - nrow(dat.out), "points outside range margins"))
    
  }
  
  return(list(dat.out, range.traits))
  
}


# used in 04 ###################################################################

## Climate niche functions #####################################################
## Written by Olivier Broennimann and Blaise Petitpierre. Departement of Ecology and Evolution (DEE). 
## University of Lausanne. Switzerland. April 2012.

### Modified by Regan Early
### Adapted by Margaret Bolton, July 2023

# Functions to calculate environmental niche position

# TODO: rewrite descriptions
## grid.clim(climate.pca.scores, species.pca.scores, R, threshold.species, threshold.env) 
## use the scores of an ordination (or SDM predictions) and create a grid species.density of RxR pixels 
## (or a vector of R pixels when using scores of dimension 1 or SDM predictions) with occurrence densities
## Only scores of one, or two dimensions can be used 

## climate.pca.scores = scores for the whole study area, 
## species.pca.scores = scores for occurrences of the species in the ordination
## samples.pca.scores = subset of scores for occurrences of the species in the ordination 
## R                  = resolution of the grid to be outputted
## threshold.species  = quantile of species density at species occurences used as a threshold to exclude low species density values
## threshold.env      = quantile of environmental density at all study sites used as a threshold to exclude low environmental density values

kd_prep <- function(clim.raw = BIO_050612, 
                    spat.dat, 
                    samp.dat, 
                    pca.full = PCA_Full,
                    bioclim.full.df = Bioclim_DF,
                    density.resolution = 10/60){
  list.out <- list()
  
  # species prep
  hostlist <- unique(samp.dat$HostCorrectedName)
  species.out <- lapply(hostlist, species_kd_prep, 
                        spat.dat = spat.dat,
                        samp.dat = samp.dat,
                        clim.raw = clim.raw,
                        pca.full = pca.full,
                        density.resolution = density.resolution)
  names(species.out) <- hostlist
  species.out <- purrr::list_transpose(species.out, simplify = FALSE)
  
  list.out$pca.species <- species.out$pca.species
  list.out$pca.samples <- species.out$pca.samples
  
  # climate prep
  # output for gridclim
  list.out$pca.climate <- ade4::suprow(pca.full, bioclim.full.df[, Clim_Variables])$lisup     #The pca scores for all climate
  
  return(list.out)
}

species_kd_prep <- function(host,
                            spat.dat = spat.dat,
                            samp.dat = samp.dat,
                            clim.raw = clim.raw, 
                            pca.full = pca.full,
                            density.resolution){
  list.out <- list()
  
  # Narrow down data to selected species
  spat.dat <- spat.dat %>%
    dplyr::filter(species == host) %>%
    dplyr::select("decimalLongitude", "decimalLatitude")
  
  samp.dat <- samp.dat %>%
    filter(HostCorrectedName == host)
  
  # (partially) account for sampling bias by gridding species occurence data and extracting filled cells
  raster.count <- terra::rast(extent = terra::ext(c(-180,180,-90,90)),
                              resolution = density.resolution) 
  
  raster.count <- samp.dat %>% 
    dplyr::select(matches("Longitude|Latitude")) %>% 
    dplyr::rename(decimalLongitude = Longitude,
                  decimalLatitude = Latitude) %>% # TODO: use str_detect
    rbind(spat.dat) %>% 
    as.matrix() %>% 
    terra::rasterize(y = raster.count,
                     fun = sum)
  
  # TODO: could be more robust, but serves for now
  if(all(terra::res(raster.count) < terra::res(clim.raw))) {
    bioclim.occurrences.df <- terra::extract(x = clim.raw,
                                             y = terra::xyFromCell(raster.count, terra::cells(raster.count)))
  } else if(all(terra::res(raster.count) == terra::res(clim.raw))){
    bioclim.occurrences.df <- terra::extract(x = clim.raw,
                                             y = terra::cells(raster.count))
  } else if(all(terra::res(raster.count) > terra::res(clim.raw))){
    bioclim.occurrences.df <- terra::extract(x = terra::aggregate(x = clim.raw, 
                                                                  fact = terra::res(raster.count) / terra::res(clim.raw),
                                                                  fun = "mean"),
                                             y = terra::xyFromCell(raster.count, terra::cells(raster.count)))
  } else warning("no method for differing res factors")
  # extract bioclim variables for points where species occurs
  # bioclim.occurrences.df <- terra::extract(x = clim.raw,
  #                                          y = terra::xyFromCell(raster.count, terra::cells(raster.count)))
  
  # extract bioclim variables for points where parasites are sampled
  bioclim.parasites.df <- terra::extract(x = clim.raw,
                                         y = samp.dat[c("Longitude", "Latitude")]) %>% 
    cbind(samp.dat)
  
  # output for grid.clim
  list.out$pca.species <- drop_na(ade4::suprow(pca.full, bioclim.occurrences.df[, Clim_Variables])$lisup)    #The pca scores for current species
  list.out$pca.samples   <- cbind(ade4::suprow(pca.full, bioclim.parasites.df[, Clim_Variables])$lisup, bioclim.parasites.df)   # pca scores for current species from gmpd
  
  return(list.out)
}

clim_density <- function(climate.pca.scores,
                         species.pca.scores,
                         samples.pca.scores,
                         R = 100,
                         threshold.species = 0,
                         threshold.env = 0){
  
  list.out <- list()
  mask.xy <- expand.grid(x = seq(0.01, 1, by = 0.01), y = seq(0.01, 1, by = 0.01))
  
  # normalise pca values
  xmin <- min(climate.pca.scores["Axis1"])
  xmax <- max(climate.pca.scores["Axis1"])
  ymin <- min(climate.pca.scores["Axis2"])
  ymax <- max(climate.pca.scores["Axis2"])
  
  climate.pca.normalised <- data.frame(cbind((climate.pca.scores["Axis1"] - xmin)/abs(xmax - xmin), 
                                             (climate.pca.scores["Axis2"] - ymin)/abs(ymax - ymin)))
  
  # calculate values of H to match adehabitatHR::kernelUD use of the ad hoc method
  climate.H <- (sqrt(0.5*(var(climate.pca.normalised["Axis1"]) + var(climate.pca.normalised["Axis2"]))))*(nrow(climate.pca.normalised)^-(1/6))
  
  # calculate the density of occurrences in a grid of RxR pixels along the score gradients
  # using a gaussian kernel density function, with RxR bins.
  
  climate.density <- MASS::kde2d(x = climate.pca.normalised[,"Axis1"],
                                 y = climate.pca.normalised[,"Axis2"],
                                 n = R,
                                 h = c(climate.H*4, climate.H*4),
                                 lims = c(range(mask.xy$x), range(mask.xy$y)))
  
  # rescale density to the number of sites in climate.pca.scores
  # or the number of occurrences in species.pca.scores
  climate.density$uncorrected  <- climate.density$z*nrow(climate.pca.scores)/sum(climate.density$z)
  
  # Fill grid using pca scores 
  pca.breaks <- data.frame(Axis1 = seq(from = min(climate.pca.scores["Axis1"]), # breaks on score gradient 1
                                       to   = max(climate.pca.scores["Axis1"]), 
                                       length.out = R),
                           Axis2 = seq(from = min(climate.pca.scores["Axis2"]), # breaks on score gradient 2
                                       to   = max(climate.pca.scores["Axis2"]), 
                                       length.out = R))
  
  climate.pca.image <- points_to_image(climate.pca.scores, pca.breaks)
  
  # calculate threshold density
  climate.density.threshold <- quantile(as.vector(climate.density$uncorrected[which(climate.pca.image == 1)]), threshold.env)  
  
  # remove tiny values generated by kernel density
  climate.density$uncorrected[climate.density$uncorrected < climate.density.threshold] <- 0
  
  # climate density needed by species_clim_density
  hostlist <- names(species.pca.scores)
  species.out <- lapply(hostlist, 
                        species_clim_density,
                        climate.pca.scores = climate.pca.scores,
                        species.pca.scores = species.pca.scores,
                        samples.pca.scores = samples.pca.scores,
                        climate.density = climate.density)
  names(species.out) <- hostlist
  species.out <- purrr::list_transpose(species.out, simplify = FALSE)
  
  list.out <- species.out
  list.out$climate.density <- climate.density
  
  return(list.out)
}

species_clim_density <- function(host,
                                 climate.pca.scores = climate.pca.scores,
                                 species.pca.scores = species.pca.scores,
                                 samples.pca.scores = samples.pca.scores,
                                 R = 100,
                                 threshold.species = 0,
                                 climate.density = climate.density){
  list.out <- list()
  mask.xy <- expand.grid(x = seq(0.01, 1, by = 0.01), y = seq(0.01, 1, by = 0.01))
  
  # subset species inputs
  species.pca.scores <- species.pca.scores[[host]]
  samples.pca.scores <- samples.pca.scores[[host]]
  
  # normalise pca values
  xmin <- min(climate.pca.scores["Axis1"])
  xmax <- max(climate.pca.scores["Axis1"])
  ymin <- min(climate.pca.scores["Axis2"])
  ymax <- max(climate.pca.scores["Axis2"])
  
  species.pca.normalised <- data.frame(cbind((species.pca.scores["Axis1"] - xmin)/abs(xmax - xmin), 
                                             (species.pca.scores["Axis2"] - ymin)/abs(ymax - ymin))) 
  
  # calculate values of H to match adehabitatHR::kernelUD use of the ad hoc method
  species.H <- (sqrt(0.5*(var(species.pca.normalised["Axis1"]) + var(species.pca.normalised["Axis2"]))))*(nrow(species.pca.normalised)^-(1/6))
  
  # calculate the density of occurrences in a grid of RxR pixels along the score gradients
  # using a gaussian kernel density function, with RxR bins.
  
  species.density <- MASS::kde2d(x = species.pca.normalised[,"Axis1"],
                                 y = species.pca.normalised[,"Axis2"],
                                 n = R,
                                 h = c(species.H*4, species.H*4),
                                 lims = c(range(mask.xy$x), range(mask.xy$y)))
  
  # rescale density to the number of sites in climate.pca.scores
  # or the number of occurrences in species.pca.scores
  species.density$uncorrected <- species.density$z*nrow(species.pca.scores)/sum(species.density$z)
  
  # Fill grid using pca scores 
  pca.breaks <- data.frame(Axis1 = seq(from = min(climate.pca.scores["Axis1"]), # breaks on score gradient 1
                                       to   = max(climate.pca.scores["Axis1"]), 
                                       length.out = R),
                           Axis2 = seq(from = min(climate.pca.scores["Axis2"]), # breaks on score gradient 2
                                       to   = max(climate.pca.scores["Axis2"]), 
                                       length.out = R))
  
  species.pca.image <- points_to_image(species.pca.scores, pca.breaks)
  
  # calculate threshold density
  species.density.threshold <- quantile(as.vector(species.density$uncorrected[which(species.pca.image == 1)]), threshold.species)
  
  # remove tiny values generated by kernel density
  species.density$uncorrected[species.density$uncorrected < species.density.threshold] <- 0 
  
  # scale between 0:1 for comparison with other species
  species.density$uncorrected <- species.density$uncorrected/max(species.density$uncorrected)	
  
  # record density as presence absence
  species.density$presence <- species.density$uncorrected
  species.density$presence[species.density$presence > 0] <- 1
  
  # correct for environment prevalence
  # TODO: probably shouldn't be using the uncorrected like this!!
  species.density$corrected <- species.density$uncorrected/climate.density$uncorrected
  
  # remove n/0 situations
  species.density$corrected[is.na(species.density$corrected)] <- 0
  species.density$corrected[species.density$corrected == "Inf"] <- 0
  
  # rescale between [0:1] for comparison with other species (again)
  species.density$corrected <- species.density$corrected/max(species.density$corrected)	
  
  # get peak of kernel density
  species.density$peak_uncorrected <- which(species.density$uncorrected == 1, arr.ind = T)
  species.density$peak_corrected <- which(species.density$corrected == 1, arr.ind = T)
  
  # parasite data time
  # get positions within climate pca from values for parasite data
  samples.pca.loci <- points_to_indices(samples.pca.scores, pca.breaks)
  
  # calculate distances and angle, and extract density for each parasite sample 
  samples.out <- distangles(loci = samples.pca.loci, 
                            origin = species.density$peak_uncorrected) %>% 
    dplyr::bind_cols(UncorrectedDensity = species.density$uncorrected[samples.pca.loci], 
                     samples.pca.scores) %>% 
    dplyr::rename(UncorrectedDistance = distance,
                  UncorrectedAngle = angle)
  
  samples.out <- distangles(loci = samples.pca.loci, 
                            origin = species.density$peak_corrected) %>% 
    dplyr::bind_cols(CorrectedDensity = species.density$corrected[samples.pca.loci], 
                     samples.out, 
                     samples.pca.loci) %>% 
    dplyr::rename(CorrectedDistance = distance,
                  CorrectedAngle = angle) 
  
  # output
  list.out$samples.out <- samples.out
  list.out$species.density <- species.density
  
  return(list.out)
}

# rotate_matrix <- function(mat) t(mat[nrow(mat):1,,drop = FALSE])


points_to_image <- function(pts, extent){
  img <- matrix(0, nrow = nrow(extent), ncol = nrow(extent))
  interval1 <- findInterval(pts$Axis1, extent$Axis1)
  interval2 <- findInterval(pts$Axis2, extent$Axis2)
  xy <- cbind(interval1, interval2)
  img[xy] <- 1
  return(img)
}


points_to_indices <- function(pts, extent){
  interval1 <- findInterval(pts$Axis1, extent$Axis1)
  interval2 <- findInterval(pts$Axis2, extent$Axis2)
  img.vals <- cbind(interval1, interval2)
  return(img.vals)
}


dists <- function(loci, origin){
  D <- sqrt((loci[, 1] - origin[, 1])^2+(loci[, 2] - origin[, 2])^2)
  return(D)
}

distangles <- function(loci, origin){
  D <- dists(loci, origin)
  
  loci[,1] <- loci[,1] - origin[,1]
  loci[,2] <- loci[,2] - origin[,2]
  
  radians <- atan2(loci[,1], loci[,2])
  
  degrees <- radians * (180/pi)
  
  DA <- cbind(D, degrees)
  colnames(DA) <- c("distance", "angle")
  return(DA)
}

# used in 05 ###################################################################
basic_barplot <- function(dat = GMPD_Climate_Data, xvar = ParClass, yvar = Prevalence){
  
  dat <- dat %>% 
    dplyr::mutate(AbsLatitude = abs(Latitude))
  
  # summarise for error bars
  plot_data <- dat %>% 
    dplyr::group_by({{xvar}}) %>% 
    dplyr::summarise(
      HostsSampledMax = mean(HostsSampled) + sd(HostsSampled),
      LatitudeMax     = mean(Latitude) + sd(Latitude),
      AbsLatitudeMax  = mean(AbsLatitude) + sd(AbsLatitude),
      PrevalenceMax   = mean(Prevalence) + sd(Prevalence),
      
      HostsSampledMin = mean(HostsSampled) - sd(HostsSampled),
      LatitudeMin     = mean(Latitude) - sd(Latitude),
      AbsLatitudeMin  = mean(AbsLatitude) - sd(AbsLatitude),
      PrevalenceMin   = mean(Prevalence) - sd(Prevalence),
      
      HostsSampled    = mean(HostsSampled),
      Latitude        = mean(Latitude),
      AbsLatitude     = mean(AbsLatitude),
      Prevalence      = mean(Prevalence),
      
      TotalHostsSampled = sum(HostsSampled), 
      n = n()
    ) %>% 
    dplyr::mutate(xvar2 = paste0({{xvar}}, " (", n, ")"),
                  xvar3 = paste0({{xvar}}, " (", TotalHostsSampled, ")"))
  
  # data for points
  jitter_data <- plot_data %>% 
    dplyr::select({{xvar}}, xvar2) %>% 
    right_join(dat, by = join_by({{xvar}}))
  
  # Trying to get around "glue-tunelling"
  yvarmin <- paste0(as_label(enquo(yvar)), "Min")
  yvarmax <- paste0(as_label(enquo(yvar)), "Max")
  
  # shit plot  
  plot_out <- ggplot(plot_data, aes(x = xvar2, y = {{yvar}})) +
    geom_bar(stat = "identity", fill =  "#3C6E71") +
    
    geom_jitter(data = jitter_data,
                col = "#B7CECE",
                width = 0.25) +
    
    geom_errorbar(aes(ymin = !! sym(yvarmin),
                      ymax = !! sym(yvarmax)),
                  width = .2,
                  col = "#1C0F13") +
    
    xlab(as_label(enquo(xvar))) +
    
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot_out)
}

relative_likelihood <- function(model1, model2){
  aics <- c(AIC(model1), AIC(model2))
  
  return(exp((min(aics)-max(aics))/2))
}

# used in script 0? ############################################################

restrict_grid <- function(location.data, rastr){
  # c.rast <- raster::rasterize(location.data[[1]], rastr, fun = "count")
  c.rast <- terra::rasterize(y = location.data[[1]],
                   fun = length)
  
  return(length(Which(c.rast, cells = TRUE))>1) 
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





