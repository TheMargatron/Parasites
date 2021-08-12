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
library(taxize)
library(tidyverse)

GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data.csv"), header = TRUE, stringsAsFactors = FALSE)
Hostlist <- unique(GMPD_Data$HostCorrectedName)
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

GBIF_Data_Temp <- GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.01 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(basisOfRecord, "_SPECIMEN")) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE)))

# Next bit is issue, then having a look at images on a map. 

## Mapping ##
Host_Synonyms

gbif_plotter <- function(Species){
  species_dat <- filter(GBIF_Data, species == Species)
  current_dir <- here::here(paste0("GBIF cleaning/", Species))
  species_IUCN <- fortify(IUCN_Data_List[[Species]])
  
  if(!dir.exists(current_dir)){
    dir.create(current_dir)
  }
  
  png(file = paste0(current_dir, "/", Species, ".png"), width = 500, height = 500, pointsize = 12)
  
  par(mfrow = c(2, 1))
  
  print(ggplot() + coord_fixed() +
          borders("world", colour = "gray50", fill = "gray50") +
          geom_polygon(data = species_IUCN, 
                       aes(x = long, y = lat, group = group),
                       colour = "black",
                       fill = NA) +
          geom_point(data = species_dat,
                     aes(x = decimalLongitude, y = decimalLatitude),
                     colour = "blue"))
  dev.off()
}





