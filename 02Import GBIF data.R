# GBIF import and data clean
# Written by Margaret Bolton mb804(at)exeter.ac.uk 
# Data citation: "GBIF.org (02 July 2021) GBIF Occurrence Download  https://doi.org/10.15468/dl.9hewb3"

############################################## Libraries and data #############################################

library(CoordinateCleaner)
library(countrycode)
library(ggplot2)
library(here)
#library(lubridate)
#library(maps)
library(rgbif)
library(rgdal)              # read shapefiles
#library(spatialEco)
library(taxize)
library(tidyverse)

GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data.csv"), header = TRUE, stringsAsFactors = FALSE)
Hostlist <- unique(GMPD_Data$HostCorrectedName)
IUCN_Mammals <- readOGR(here::here("Data/IUCN"), "MAMMALS") #this takes a while
IUCN_Data_List <- readRDS(here::here("Data/Data back ups/IUCN_Data_List"))

############################################## Getting taxon keys #############################################
# Adding synonymous species names not picked up by taxize
Host_Synonyms <- data.frame("IUCNName" = Hostlist, "GBIFName" = Hostlist, stringsAsFactors = FALSE)
Host_Synonyms <- Host_Synonyms %>%
  mutate(GBIFName = case_when(GBIFName == "Neovison vison"   ~ "Mustela vison",
                              GBIFName == "Tragelaphus oryx" ~ "Taurotragus oryx",
                              GBIFName == "Martes pennanti"  ~ "Pekania pennanti",
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

############################################## Actual download ################################################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

# Download_Key <- occ_download(
#   pred_in("taxonKey", Taxon_Keys$usagekey),
#   pred("hasCoordinate", TRUE),
#   format = "SIMPLE_CSV"
# )
# 
# saveRDS(Download_Key, here::here("Data/GBIF/Download_Key"))
# Download_Key <- readRDS(here::here("Data/GBIF/Download_Key"))
# 
# Download_Get <- occ_download_get(Download_Key, path = here::here("Data/GBIF/"), overwrite = TRUE)
# 
# saveRDS(Download_Get, here::here("Data/GBIF/Download_Get"))
Download_Get <- readRDS(here::here("Data/GBIF/Download_Get"))

GBIF_Raw_Data <- occ_download_import(Download_Get, path = here::here("Data/GBIF/"))

############################################## Cleaning data ##################################################

GBIF_Data <- filter(GBIF_Raw_Data, countryCode != "" & countryCode != "XK" & countryCode != "ZZ")
GBIF_Data$countryCode <- countrycode(GBIF_Data$countryCode, origin = "iso2c", destination = "iso3c")

GBIF_Data <- clean_coordinates(x = GBIF_Data,
                               lon = "decimalLongitude", 
                               lat = "decimalLatitude", 
                               countries = "countryCode",
                               tests = c("capitals", "centroids", "countries", "gbif", "institutions", "seas", "zeros"), # think about whether I need duplicates or outlier tests
                               capitals_rad = 10000,
                               centroids_rad = 1000,
                               centroids_detail = "country",
                               inst_rad = 100,
                               zeros_rad = 0.5,
                               value = "clean")

## Issues ## 

c(2, 3, 4, 6, 8, 9, 11, 12, 13, 14, 31, 32, 33, 42)

## Coordinate uncertainty ##

# Raster resolution used later in analysis is 2.5 arcminutes (roughly 0.042 degrees)
# This corresponds to approx 4625m at the equator
# From histograms its clear that 5000m is a commonly estimated uncertainty for coordinates
# Proportionally few records have an uncertainty over 5000m
# Need to explain the logic jump a bit but it seems to make sense to use 5000m as the max uncertainty
# Wouldn't lose too much data and it means the true coordinates at worst will only be one raster cell over
# Because of spatial autocorrelation this shouldn't have too much effect on analyses

GBIF_Data %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, no max") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 1000000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 1,000,000m") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 100000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 100,000m") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 10000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 10,000m") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 7000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 7,000m") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 4500) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 4,500m") %>%
  print()

## Coordinate precision ##

# Raster resolution used later in analysis is 2.5 arcminutes (roughly 0.042 degrees)
# I could use similar reasoning as I did for coordinate uncertainty and apply 0.05 as a cut off

GBIF_Data %>%
  pull(coordinatePrecision) %>%
  hist(main = "Histogram of coordinate precision, no max") %>%
  print()

GBIF_Data %>%
  filter(coordinatePrecision <= 0.1) %>%
  pull(coordinatePrecision) %>%
  hist(main = "Histogram of coordinate precision, max 0.1") %>%
  print()

GBIF_Data %>%
  filter(coordinatePrecision <= 0.05) %>%
  pull(coordinatePrecision) %>%
  hist(main = "Histogram of coordinate precision, max 0.05") %>%
  print()

GBIF_Data %>%
  filter(coordinatePrecision <= 0.01) %>%
  pull(coordinatePrecision) %>%
  hist(main = "Histogram of coordinate precision, max 0.01") %>%
  print()

## Event date ##
 ## Not sure whether to filter by date
GBIF_Data %>%
  mutate(eventDate = case_when(eventDate == "" ~ NA_character_,
                               TRUE            ~ eventDate)) %>%
  pull(eventDate) %>%
  as.Date() %>%
  hist(breaks = "years",
       main = "Histogram of event dates")

## Basis of record ##

# https://data-blog.gbif.org/post/living-specimen-to-preserved-specimen-understanding-basis-of-record/

# base map for plotting species data

base_map <- ggplot() + coord_fixed() +
  borders("world", colour = "gray50", fill = "gray50") 

# Using Chrysocyon brachyurus as an example for specimens in basis of record

base_map +
  geom_polygon(data = fortify(IUCN_Data_List[["Chrysocyon brachyurus"]]), 
               aes(x = long, y = lat, group = group),
               colour = "black",
               fill = NA) +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "blue") +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus" & basisOfRecord == "FOSSIL_SPECIMEN"),
             aes(x = decimalLongitude, y = decimalLatitude),
             color = "red") +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus" & basisOfRecord == "PRESERVED_SPECIMEN"),
             aes(x = decimalLongitude, y = decimalLatitude),
             color = "green") +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus" & basisOfRecord == "LIVING_SPECIMEN"),
             aes(x = decimalLongitude, y = decimalLatitude),
             color = "orange")

## Locality ##

base_map +
  geom_polygon(data = fortify(IUCN_Data_List[["Chrysocyon brachyurus"]]), 
               aes(x = long, y = lat, group = group),
               colour = "black",
               fill = NA) +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "blue") +
  
  geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus" & 
                             str_detect(locality, regex("zoo", ignore_case = TRUE))),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "red")

# Filtering

GBIF_Data <- GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.01 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(basisOfRecord, "_SPECIMEN")) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE))) %>%
  filter(basisOfRecord != "MATERIAL_SAMPLE")

## Mapping ##

GBIF_Plots <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base")
names(GBIF_Plots) <- Host_Synonyms$IUCNName

## Running outliers for mapping ## 

# Distance method seems a bit bats. Look at Leopardus geoffroyi [[3]]
GBIF_Outliers_dist$cc_outl2 <- cc_outl(x = GBIF_Data,
                                      lon = "decimalLongitude",
                                      lat = "decimalLatitude",
                                      species = "species",
                                      method = "distance",
                                      tdi = 1000,
                                      value = "flagged")

GBIF_Outliers_quantile <- clean_coordinates(x = GBIF_Data,
                                        lon = "decimalLongitude",
                                        lat = "decimalLatitude",
                                        species = "species",
                                        tests = c("outliers"),
                                        outliers_method = "quantile",
                                        outliers_mtp = 5)

GBIF_Plots_Test <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Outliers_dist, data_type = "tested")
names(GBIF_Plots) <- Host_Synonyms$IUCNName

temp <- apply(Host_Synonyms, MARGIN = 1, species_outlier, dat = GBIF_Data)
temp <- bind_rows(temp)

temp2 <- apply(Host_Synonyms, MARGIN = 1, species_outlier, dat = GBIF_Data)
temp2 <- bind_rows(temp2)

### Temporary plotting to get outlier parameters

# iucn test
Species <- c("Canis lupus",
             "Mustela erminea",
             "Mustela nivalis",
             "Mustela lutreola")

species_dat <- GBIF_Data %>%
  mutate(species = case_when(species == "Mustela vison"    ~ "Neovison vison",
                             species == "Taurotragus oryx" ~ "Tragelaphus oryx",
                             species == "Pekania pennanti" ~ "Martes pennanti",
                             TRUE                          ~ species)) %>%
  filter(species == Species[1]) %>%
  rename(binomial = species)

species_dat$out <-  cc_iucn(x = species_dat,
                            range = IUCN_Data_List[[Species[1]]],
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            species = "binomial",
                            buffer = 1,
                            value = "flagged")

base_map +
  geom_polygon(data = fortify(IUCN_Data_List[[Species[1]]]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen3",
               fill = "palegreen3") +
#  split points into two so I can easily switch between them and ensures outliers are plotted on top
  geom_point(data = filter(species_dat, out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_point(data = filter(species_dat, !out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  ggtitle(Species[1])
beep(2)


#outlier test
Species <- "Connochaetes gnou"

species_dat <- filter(GBIF_Data, species == Species)

species_dat$out <- cc_outl(x = species_dat,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           method = "quantile",
                           mltpl = 5,
                           value = "flagged")

base_map +
  geom_polygon(data = fortify(IUCN_Data_List[[Species]]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen4",
               fill = "palegreen4") +
  #  split points into two so I can easily switch between them and ensures outliers are plotted on top
  geom_point(data = filter(species_dat, out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_point(data = filter(species_dat, !out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  ggtitle(Species)
beep(2)

## Awkward species ####

### Alces alces ####

### Cervus elaphus ####
GBIF_Plots[[5]]

base_map +
  geom_polygon(data = fortify(IUCN_Mammals[IUCN_Mammals$binomial == "Cervus elaphus", ]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen4",
               fill = "palegreen4") +
  
  geom_polygon(data = fortify(IUCN_Mammals[IUCN_Mammals$binomial == "Cervus canadensis", ]), 
               aes(x = long, y = lat, group = group),
               colour = "forestgreen",
               fill = "forestgreen") +
  
#  split points into two so I can easily switch between them and ensures outliers are plotted on top
  
  geom_point(data = filter(GBIF_Data, species == "Cervus elaphus"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  geom_point(data = filter(GMPD_Raw_Data , HostCorrectedName == "Cervus elaphus"),
             aes(x = Longitude, y = Latitude),
             colour = "navy") +
  
  ggtitle("Cervus elaphus")
beep(2)

### Meles meles ####
GBIF_Plots[[6]]

# Meles meles is the only one recorded in GMPD_Raw_Data
GMPD_Raw_Data %>% filter(str_detect(HostCorrectedName, "Meles")) %>% pull(HostCorrectedName) %>% unique()

# Recorded separately in GBIF
Meles <- taxize::get_gbifid_(c("Meles meles", "Meles anakuma", "Meles leucurus"), method = "backbone")

# M anakuma, M. Leucurus
base_map +
  geom_polygon(data = fortify(IUCN_Mammals[IUCN_Mammals$binomial == "Meles anakuma", ]), 
               aes(x = long, y = lat, group = group),
               colour = "navy",
               fill = NA) +
  
  geom_polygon(data = fortify(IUCN_Mammals[IUCN_Mammals$binomial == "Meles leucurus", ]), 
               aes(x = long, y = lat, group = group),
               colour = "firebrick",
               fill = NA) +
  
  geom_polygon(data = fortify(IUCN_Mammals[IUCN_Mammals$binomial == "Meles meles", ]), 
               aes(x = long, y = lat, group = group),
               colour = "forestgreen",
               fill = NA) +
  
  geom_point(data = filter(GBIF_Data, species == "Meles meles"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_point(data = filter(GMPD_Raw_Data , HostCorrectedName == "Meles meles"),
             aes(x = Longitude, y = Latitude),
             colour = "orange") +
  
  #  split points into two so I can easily switch between them and ensures outliers are plotted on top
  ggtitle("Meles meles")
beep(2)

### Canis lupus ####
GBIF_Plots[[10]]

subsp <- c("familiaris", "dingo")

species_dat <- GBIF_Data %>%
  filter(species == "Canis lupus") %>%
  filter(!str_detect(verbatimScientificName, "familiaris|dingo")) %>%   # removed domestic dogs
  mutate(subsp = case_when(str_detect(verbatimScientificName, "arctos")       ~ "arctos",
                           str_detect(verbatimScientificName, "crassodon")    ~ "crassodon",
                           str_detect(verbatimScientificName, "arabs")        ~ "arabs",
                           str_detect(verbatimScientificName, "pallipes")     ~ "pallipes",
                           str_detect(verbatimScientificName, "occidentalis") ~ "occidentalis",
                           str_detect(verbatimScientificName, "chanco")       ~ "chanco",
                           str_detect(verbatimScientificName, "albus")        ~ "albus",
                           str_detect(verbatimScientificName, "baileyi")      ~ "baileyi",
                           str_detect(verbatimScientificName, "italicus")     ~ "italicus",
                           str_detect(verbatimScientificName, "lycaon")       ~ "lycaon",
                           str_detect(verbatimScientificName, "signatus")     ~ "signatus",
                           str_detect(verbatimScientificName, "rufus")        ~ "rufus",
                           TRUE                                               ~ "lupus")) %>%
  mutate(subsp = as.factor(subsp))

my_colors <- c("chartreuse", "chartreuse", "chartreuse",
               "orange", "orange", "orange",
               "deeppink",
               "darkcyan", "darkcyan",
               "turquoise", "turquoise", "turquoise")
my_shapes <- c(0,1,2,0,1,2,2,1,2,0,1,2)

base_map +
  geom_polygon(data = fortify(IUCN_Mammals[IUCN_Mammals$binomial == "Canis lupus", ]), 
               aes(x = long, y = lat, group = group),
               colour = "white",
               fill = "white") +
  
  geom_point(data = filter(species_dat, subsp != "lupus"),
             aes(x = decimalLongitude, y = decimalLatitude, colour = subsp, shape = subsp)) +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  
  #  split points into two so I can easily switch between them and ensures outliers are plotted on top
  ggtitle("Canis lupus")
beep(2)

### Mustela erminea ####
GBIF_Plots[[12]]
# Mustela erminea can be confused with M. nivalis. 
# The M. erminea samples outside its range generally fall within the M nivalis range

base_map +
  geom_polygon(data = fortify(IUCN_Data_List[["Mustela erminea"]]), 
               aes(x = long, y = lat, group = group),
               colour = "navy",
               fill = NA) +
  
  geom_polygon(data = fortify(IUCN_Data_List[["Mustela nivalis"]]), 
               aes(x = long, y = lat, group = group),
               colour = "firebrick",
               fill = NA)

# Using the same buffer method as CoordinateCleaner for consistency
M_erminea_Buff <- rgeos::gBuffer(IUCN_Data_List[["Mustela erminea"]], byid = TRUE, width = 1)
M_nivalis_overlap <- IUCN_Data_List[["Mustela nivalis"]] - M_erminea_Buff 
M_nivalis_overlap$binomial <- factor("Mustela erminea")

species_dat <- GBIF_Data %>%
  filter(species == "Mustela erminea") %>%
  rename(binomial = species)

species_dat$out <-  cc_iucn(x = species_dat,
                            range = M_nivalis_overlap,
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            species = "binomial",
                            buffer = 0.5,
                            value = "flagged")

base_map +
  geom_polygon(data = fortify(IUCN_Data_List[["Mustela nivalis"]]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen4",
               fill = "palegreen4") +
  
  geom_polygon(data = fortify(IUCN_Data_List[["Mustela erminea"]]), 
               aes(x = long, y = lat, group = group),
               colour = "white",
               fill = "white") +
  #  split points into two so I can easily switch between them and ensures outliers are plotted on top
  geom_point(data = filter(species_dat, !out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_point(data = filter(species_dat, out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  ggtitle("Mustela erminea")
beep(2)

# Doesn't catch the ones in America but having checked online records they are both subsp. Richardsonii so okay to keep

