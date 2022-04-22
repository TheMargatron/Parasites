# GBIF import and data clean
# Written by Margaret Bolton mb804(at)exeter.ac.uk 
# Data citation: "GBIF.org (6 October 2021) GBIF Occurrence Download  https://www.gbif.org/occurrence/download/0018590-210914110416597"

############################################## Libraries and data #############################################

library(CoordinateCleaner)
library(countrycode)
library(ggplot2)
library(here)
#library(lubridate)
#library(maps)
library(rgbif)
library(rgdal)              # read shapefiles
library(rgeos)
library(rnaturalearthdata)
library(taxize)
library(tidyverse)
library(tmap)
library(beepr)
source(here::here("Functions.R"))

sf::sf_use_s2(FALSE) # For "invalid spherical geometry" errors
tmap_mode("view")

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data_01.csv"), header = TRUE, stringsAsFactors = FALSE)

Legend_Text <- sort(unique(read.csv(here::here("Data/Data back ups/Native_DF_01.csv"), header = TRUE)$Status))

IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))
#IUCN_Orders <- readRDS(here::here("Data/Data back ups/IUCN_Orders_01"))

Hostlist <- sort(unique(GMPD_Data$HostCorrectedName))

# Getting taxon keys ##########################################################################################
# Adding synonymous species names not picked up by taxize
Host_Synonyms <- data.frame("IUCNName" = Hostlist, "GBIFName" = Hostlist, stringsAsFactors = FALSE)
Host_Synonyms <- Host_Synonyms %>%
  mutate(GBIFName = case_when(GBIFName == "Martes pennanti"  ~ "Pekania pennanti",
                              GBIFName == "Melogale subaurantiaca"  ~ "Melogale moschata",
                              GBIFName == "Neovison vison"   ~ "Mustela vison",
                              GBIFName == "Tragelaphus oryx" ~ "Taurotragus oryx",
                              TRUE                       ~ GBIFName))

Taxon_Keys <- taxize::get_gbifid_(Host_Synonyms$GBIFName, method = "backbone")
Taxon_Keys <- lapply(Host_Synonyms$GBIFName, function(name) {
  Taxon_Keys[[name]]["GBIFName"] <- name
  return(Taxon_Keys[[name]])
})

Taxon_Keys <- Taxon_Keys %>%
  bind_rows() %>%
  filter(class == "Mammalia") %>%
  filter(status == "ACCEPTED" & matchtype == "EXACT")

# Actual download #############################################################################################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

Download_Key <- occ_download(
  pred_in("taxonKey", Taxon_Keys$usagekey),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV"
)

saveRDS(Download_Key, here::here("Data/GBIF/Download_Key"))
Download_Key <- readRDS(here::here("Data/GBIF/Download_Key"))

Download_Get <- occ_download_get(Download_Key, path = here::here("Data/GBIF/"), overwrite = TRUE)

saveRDS(Download_Get, here::here("Data/GBIF/Download_Get"))
Download_Get <- readRDS(here::here("Data/GBIF/Download_Get"))

GBIF_Raw_Data <- occ_download_import(Download_Get, path = here::here("Data/GBIF/"))
nrow(GBIF_Raw_Data) #3277366
Host_Synonyms[!Host_Synonyms$GBIFName %in% unique(GBIF_Raw_Data$species),]
# Cervus canadensis is recorded as Cervus elaphus subsp canadensis so not actually missing data

# GBIF_Raw_Plots_00 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Raw_Data, data_type = "base", range.polygon = IUCN_Native_Data)
# names(GBIF_Raw_Plots_00) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_Raw_Plots_00.pdf'), width = 10, height = 7)
# GBIF_Raw_Plots_00
# dev.off()

rm(Download_Key, Download_Get)
# Cleaning data ###############################################################################################

GBIF_Data <- filter(GBIF_Raw_Data, countryCode != "" & countryCode != "XK" & countryCode != "ZZ")
GBIF_Data$countryCode <- countrycode(GBIF_Data$countryCode, origin = "iso2c", destination = "iso3c")
nrow(GBIF_Data) #3276768

GBIF_Data <- clean_coordinates(x = GBIF_Data,
                               lon = "decimalLongitude", 
                               lat = "decimalLatitude", 
                               countries = "countryCode",
                               tests = c("capitals", "centroids", "countries", "gbif", "institutions", "zeros"), 
                               capitals_rad = 10000,
                               centroids_rad = 1000,
                               centroids_detail = "country",
                               inst_rad = 100,
                               zeros_rad = 0.5,
                               value = "clean")

nrow(GBIF_Data) #3068592
gc()

# GBIF_bor_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "basisOfRecord", range.polygon = IUCN_Native_Data)
# names(GBIF_bor_Plots_01) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_bor_Plots_01.pdf'), width = 10, height = 7)
# GBIF_bor_Plots_01
# dev.off()
# 
# GBIF_base_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base", range.polygon = IUCN_Native_Data)
# names(GBIF_base_Plots_01) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_base_Plots_01.pdf'), width = 10, height = 7)
# GBIF_base_Plots_01
# dev.off()
#
# GBIF_epithet_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "infraspecificEpithet", range.polygon = IUCN_Native_Data)
# names(GBIF_epithet_Plots_01) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_epithet_Plots_01.pdf'), width = 10, height = 7)
# GBIF_epithet_Plots_01
# dev.off()


## Issues #### 

### Coordinate uncertainty ####

# Raster resolution used later in analysis is 2.5 arcminutes (roughly 0.042 degrees)
# This corresponds to approx 4625m at the equator
# From histograms its clear that 5000m is a commonly estimated uncertainty for coordinates
# Proportionally few records have an uncertainty over 5000m
# Need to explain the logic jump a bit but it seems to make sense to use 5000m as the max uncertainty
# Wouldn't lose too much data and it means the true coordinates at worst will only be one raster cell over
# Because of spatial autocorrelation this shouldn't have too much effect on analyses

GBIF_Data %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, no max")

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 10000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 10,000m")

### Coordinate precision ####

# Raster resolution used later in analysis is 2.5 arcminutes (roughly 0.042 degrees)
# I could use similar reasoning as I did for coordinate uncertainty and apply 0.05 as a cut off
# But I using 0.02 only removes a marginal amount more while ensuring greater raster cell accuracy

GBIF_Data %>%
  pull(coordinatePrecision) %>%
  hist(main = "Histogram of coordinate precision, no max")

GBIF_Data %>%
  filter(coordinatePrecision <= 0.2) %>%
  pull(coordinatePrecision) %>%
  hist(main = "Histogram of coordinate precision, max 0.2") 

nrow(GBIF_Data[GBIF_Data$coordinatePrecision <= 0.02,]) / nrow(GBIF_Data)
nrow(GBIF_Data[GBIF_Data$coordinatePrecision <= 0.05,]) / nrow(GBIF_Data)

### Event date ####
 ## Not sure whether to filter by date
GBIF_Data %>%
  filter(!is.na(eventDate)) %>%
  pull(eventDate) %>%
  as.Date() %>%
  hist(breaks = "years",
       main = "Histogram of event dates") 

### Basis of record ####

# https://data-blog.gbif.org/post/living-specimen-to-preserved-specimen-understanding-basis-of-record/

# Using Chrysocyon brachyurus as an illustrative example for specimen basis of record

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Chrysocyon brachyurus",]), 
               aes(x = long, y = lat, group = group),
               colour = "black",
               fill = NA) +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus"),
             aes(x = decimalLongitude, y = decimalLatitude, colour = basisOfRecord)) +
  
  ggtitle("Chrysocyon brachyurus basis of record")

### Locality ####

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Chrysocyon brachyurus",]), 
               aes(x = long, y = lat, group = group),
               colour = "black",
               fill = NA) +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus"),
             aes(x = decimalLongitude, y = decimalLatitude, 
                 colour = (species == "Chrysocyon brachyurus" & 
                   str_detect(locality, regex("zoo", ignore_case = TRUE))))) +
  
  scale_colour_discrete(name = "Zoo") +
  ggtitle("Chrysocyon brachyurus zoos")

## Filtering ####

GBIF_Data <- GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.02 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(basisOfRecord, "_SPECIMEN")) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE))) %>%
  filter(basisOfRecord != "MATERIAL_SAMPLE") %>%
  filter(!str_detect(basisOfRecord, "UNKNOWN"))

nrow(GBIF_Data) #2491767

# GBIF_Base_Plots_02 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base", range.polygon = IUCN_Native_Data)
# names(GBIF_Base_Plots_02) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_Base_Plots_02.pdf'), width = 10, height = 7)
# GBIF_Base_Plots_02
# dev.off()

## Subgroups and outliers ####
# May move prep steps to 02.1 if I don't make other sub-scripts
GBIF_Data <- GBIF_Data %>%
  mutate(infraspecificEpithet = case_when(infraspecificEpithet == "" ~ NA_character_,
                                          TRUE ~ infraspecificEpithet)) %>%
  mutate(subgroup = infraspecificEpithet)

GBIF_Data <- as.data.frame(GBIF_Data)

GBIF_Spatial <- SpatialPointsDataFrame(coords      = GBIF_Data[, c("decimalLongitude", "decimalLatitude")],
                                       data        = GBIF_Data[, names(GBIF_Data)[!names(GBIF_Data) %in% c("decimalLongitude", "decimalLatitude")]], 
                                       proj4string = CRS(proj4string(IUCN_Native_Data)))

source(here::here("02.1Subgrouping GBIF.R"))

nrow(GBIF_Subgroups) #

# Plotting ####
# prep
GMPD_Spatial <- SpatialPointsDataFrame(coords      = GMPD_Data[, c("Longitude", "Latitude")],
                                       data        = GMPD_Data[, names(GMPD_Data)[!names(GMPD_Data) %in% c("Longitude", "Latitude")]], 
                                       proj4string = CRS(proj4string(IUCN_Native_Data)))

# polygon legend groups:
poly_fill <- data.frame(legend = unique(IUCN_Native_Data@data$legend),
           fill_group = NA)
poly_fill <- poly_fill %>%
  mutate(fill_group = case_when(legend == "Extinct" ~ "Extinct",
                                legend == "Extinct & Reintroduced" ~ "Extinct",
                                str_detect(legend, "Extinct|Presence") ~ "Absence Likely",
                                legend == "Extant (resident)" ~ "Extant",
                                str_detect(legend, "Extant & R|Extant & I") ~ "Extant",
                                str_detect(legend, "Extant") ~ "Presence Likely"))

# borders and scale
sub.colours <- c('#ea3c67', '#3cb44b', '#ffbb19', # red, green, yellow
                 '#4388d8', '#f58231', '#42d4f4', # blue, orange, cyan
                 '#f032e6', '#469990', # magenta, teal
                 '#b372ff', '#c37e2e', '#f577a5') # lavendar, brown, pink

Host_Synonyms_temp <- Host_Synonyms %>%
  mutate(minlong = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
                                                                           bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Longitude", "min"]}),
         maxlong = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
                                                                           bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Longitude", "max"]}),
         minlat = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
                                                                          bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Latitude", "min"]}),
         maxlat = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
                                                                          bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Latitude", "max"]})) %>%
  mutate(scaling = case_when(IUCNName == "Vulpes vulpes" ~ "Eurasia",
                             maxlong - minlong > 300 ~ "Global",
                             maxlat < 38 & maxlong < 58 & minlong > -20 ~ "Africa",
                             str_detect(IUCNName, "Herpestes ichneumon|Hyaena hyaena|Panthera leo") ~ "Africa",
                             str_detect(IUCNName, "Puma concolor|Odocoileus virginianus") ~ "Global", 
                             minlat > 6 & maxlong < -12  ~ "North",
                             IUCNName == "Leopardus tigrinus" ~ "Central",
                             maxlat < 13 & maxlong < -32  ~ "South",
                             maxlong < -30 ~ "Central",
                             minlong > 91 ~ "Asia",
                             minlong > -11  & maxlong < 117 ~ "Europe",
                             TRUE ~ "Eurasia")) %>%
  dplyr::select(IUCNName, GBIFName, scaling)

bboxes <- list("Global" = matrix(c(-180, -90, 180, 90), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "Eurasia" = matrix(c(-180, -90, 180, 90), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "Africa" = matrix(c(-60, -41, 120, 49), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "South" = matrix(c(-180, -67, 0, 23), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "North" = matrix(c(-180, -8, 0, 82), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "Central" = matrix(c(-180, -38, 0, 52), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "Asia" = matrix(c(0, -8, 180, 82), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
               "Europe" = matrix(c(-37, -8, 143, 82), 
                                 ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))))


BF_temp <- apply(Host_Synonyms, MARGIN = 1, function(syn.row) {
  xy <- bboxer(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == syn.row["IUCNName"],]@bbox,
         GBIF_Spatial[GBIF_Spatial$species == syn.row["GBIFName"],]@bbox,
         IUCN_Native_Data[IUCN_Native_Data$binomial == syn.row["IUCNName"],]@bbox)

  xy <- xy[,2] - xy[,1]
  xy <- xy[1]/xy[2]
})

apply(Host_Synonyms[1,], MARGIN = 1, FUN = complete_plot, 
      dat = GMPD_Spatial, range.dat = GBIF_Spatial, range.polygon = IUCN_Native_Data)

# sort out buffer function
# make legends pretty (mostly position)
# fix legends
## extend x/ylim in one direction
# all daata at once, symbols differ
# put iucn and gbif next to each other
# send collated data


# look for mathematica trainign courses
# ask bram for mathematica code
# bes, bob o hara

