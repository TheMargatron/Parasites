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

GBIF_Data <- GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.01 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(basisOfRecord, "_SPECIMEN")) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE))) %>%
  filter(basisOfRecord != "MATERIAL_SAMPLE")

## Mapping ##

Host_Synonyms

GBIF_Plots <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter)
names(GBIF_Plots) <- Host_Synonyms$IUCNName

write.csv(Host_Synonyms, file = here::here("GBIF cleaning/GBIF_issues.csv"), row.names = TRUE)

## Going through them manually and finding species with sus data points

Sus_GBIF_Data <- c("Leopardus geoffroyi / 3 / Geoffroy's cat in the USA                 / probably outliers test",
  "Capreolus capreolus                  / 4 / Roe deer in the USA and Korea             / probably outliers test",
  "Cervus elaphus                       / 5 / Subspecies issues                         / IUCN polygons then maybe outliers",
  "Meles meles                          / 6 / Subspecies issues                         / IUCN polygons then maybe outliers",
  "Panthera pardus                      / 8 / Leopards in Europe                        / Tricky one",
  "Alces alces                          / 9 / Mooses in the UK and Europe               / Tricky one",
  "Canis lupus                          / 10/ Wolves everywhere                         / Tricky one",
  "Mustela erminea                      / 11/ Stoats in Southern Europe                 / Uncertain if sus or okay",
  "Mustela putorias                     / 13/ Polecats in the Azores, Canaries, Oceania / maybe outliers?",
  "Nyctereutes procyonoides             / 14/ Raccoon dogs in the Ireland and Europe    / Tricky one",
  "Puma concolor                        / 20/ Cougar in Europe                          / outliers",
  "Panthera leo                         / 22/ Lion in New Zealand                       / outliers? but south Asia",
  "Odocoileus hemionus                  / 25/ Mule deer in Florida                      / outliers? but maybe ok",
  "Odocoileus virginianus               / 26/ White tailed deer in Eurasia              / IUCN with buffer?",
  "Mustela vison                        / 27/ American mink everywhere                  / IUCN with buffer?",
  "Procyon lotor                        / 28/ Raccoons in Eurasia                       / IUCN with buffer?",
  "Bison bison                          / 29/ Bison in Europe and Africa                / outliers",
  "Mephitis mephitis                    / 32/ Striped skunk in Europe                   / outliers",
  "Rangifer tarandus                    / 34/ Reindeer in UK and Europe (Excl. North)   / Tricky one",
  "Ursus arctos                         / 35/ Brown bears in the UK                     / Tricky one",
  "Vulpes lagopus                       / 38/ Arctic foxes in EU, USA, Asia             / IUCN with buffer?",
  "Rupicapra rupicapra                  / 39/ Chamois in New Zealand                    / IUCN with buffer?",
  "Urocyon cinereoargenteus             / 40/ Gray fox in Europe                        / outliers",
  "Dama dama                            / 41/ Fallow deer all over the place            / IUCN with buffer?",
  "Genetta genetta                      / 43/ Genets in Europe but within IUCN          / Non-native IUCN?",
  "Felis silvestris                     / 45/ Wildcats in East Asia and southern UK     / outliers and ignore UK?",
  "Hyaena hyaena                        / 47/ Striped hyena in Namibia                  / IUCN with buffer?",
  "Mustela nivalis                      / 48/ Least weasel in Oceania                   / IUCN with buffer?",
  "Canis aureus                         / 49/ Golden jackal in Africa                   / Tricky one",
  "Cervus nippon                        / 51/ Sika deer all over the place              / Tricky one, IUCN?",
  "Martes melampus                      / 52/ Japanese martens in South Korea           / Tricky one, IUCN?",
  "Aepyceros melampus                   / 54/ Impala in USA and West Africa             / outliers or IUCN",
  "Crocuta crocuta                      / 58/ Spotted hyena in East Asia                / outliers or IUCN",
  "Taurotragus oryx                     / 61/ Eland in West Africa                      / outliers or IUCN",
  "Equus quagga                         / 62/ Plains zebra in West Africa               / outliers or ignore")



