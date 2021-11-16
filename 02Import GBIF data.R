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
library(taxize)
library(tidyverse)
library(beepr)
source(here::here("Functions.R"))

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
GMPD_Data_res <- read.csv(here::here("Data/Data back ups/GMPD_Data_res_01.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Data_cln <- read.csv(here::here("Data/Data back ups/GMPD_Data_cln_01.csv"), header = TRUE, stringsAsFactors = FALSE)

Legend_Text <- sort(unique(read.csv(here::here("Data/Data back ups/Native_DF_01.csv"), header = TRUE)$Status))

IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))
IUCN_Orders <- readRDS(here::here("Data/Data back ups/IUCN_Orders_01"))

GMPD_Data_res$FilterType <- "res"
GMPD_Data_cln$FilterType <- "cln"

GMPD_Data_full <- bind_rows(GMPD_Data_res, GMPD_Data_cln)
Hostlist <- sort(unique(GMPD_Data_full$HostCorrectedName))

# Getting taxon keys ##########################################################################################
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

# Actual download #############################################################################################
# warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")
# 
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
nrow(GBIF_Raw_Data) #2950766

GBIF_Raw_Plots_00 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Raw_Data, data_type = "base", range.polygon = IUCN_Native_Data)
names(GBIF_Raw_Plots_00) <- Host_Synonyms$IUCNName

pdf(file = here::here('GBIF Cleaning/GBIF_Raw_Plots_00.pdf'), width = 10, height = 7)
GBIF_Raw_Plots_00
dev.off()

rm(Download_Get)
# Cleaning data ###############################################################################################

GBIF_Data <- filter(GBIF_Raw_Data, countryCode != "" & countryCode != "XK" & countryCode != "ZZ")
GBIF_Data$countryCode <- countrycode(GBIF_Data$countryCode, origin = "iso2c", destination = "iso3c")
nrow(GBIF_Data) #2950037

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

nrow(GBIF_Data) #2757576
gc()

GBIF_bor_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "bor", range.polygon = IUCN_Native_Data)
names(GBIF_bor_Plots_01) <- Host_Synonyms$IUCNName

pdf(file = here::here('GBIF Cleaning/GBIF_bor_Plots_01.pdf'), width = 10, height = 7)
GBIF_bor_Plots_01
dev.off()

GBIF_base_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base", range.polygon = IUCN_Native_Data)
names(GBIF_base_Plots_01) <- Host_Synonyms$IUCNName

pdf(file = here::here('GBIF Cleaning/GBIF_base_Plots_01.pdf'), width = 10, height = 7)
GBIF_base_Plots_01
dev.off()

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
  hist(main = "Histogram of coordinate uncertainty, no max") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 10000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 10,000m") %>%
  print()

GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 5000) %>%
  pull(coordinateUncertaintyInMeters) %>%
  hist(main = "Histogram of coordinate uncertainty, max 5,000m") %>%
  print()


### Coordinate precision ####

# Raster resolution used later in analysis is 2.5 arcminutes (roughly 0.042 degrees)
# I could use similar reasoning as I did for coordinate uncertainty and apply 0.05 as a cut off
# Don't remember why I've used 0.01 later on ...

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

### Event date ####
 ## Not sure whether to filter by date
GBIF_Data %>%
  filter(!is.na(eventDate)) %>%
  pull(eventDate) %>%
  as.Date() %>%
  hist(breaks = "years",
       main = "Histogram of event dates") %>%
  print()

### Basis of record ####

# https://data-blog.gbif.org/post/living-specimen-to-preserved-specimen-understanding-basis-of-record/

# Using Chrysocyon brachyurus as an example for specimens in basis of record

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
  filter(coordinatePrecision <= 0.01 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(basisOfRecord, "_SPECIMEN")) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE))) %>%
  filter(basisOfRecord != "MATERIAL_SAMPLE") %>%
  filter(!str_detect(basisOfRecord, "UNKNOWN"))

nrow(GBIF_Data) #2237483

GBIF_Base_Plots_02 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base", range.polygon = IUCN_Native_Data)
names(GBIF_Base_Plots_02) <- Host_Synonyms$IUCNName

pdf(file = here::here('GBIF Cleaning/GBIF_Base_Plots_02.pdf'), width = 10, height = 7)
GBIF_Base_Plots_02
dev.off()

## Awkward species ####

### Aepyceros melampus ####
Aepyceros_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Aepyceros melampus", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Aepyceros melampus",]

pdf(file = here::here('GBIF Cleaning/Aepyceros melampus.pdf'), width = 10, height = 7)

GBIF_Base_Plots_02[["Aepyceros melampus"]] +
  geom_polypath(data = fortify(Aepyceros_remove_IUCN ), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) + 
  coord_map(xlim = c(0, 45), ylim = c(-40, 5))

GBIF_Base_Plots_02[["Aepyceros melampus"]] +
  geom_polypath(data = fortify(Aepyceros_remove_IUCN ), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) + 
  coord_map(xlim = c(10, 25), ylim = c(-30, -15))

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Aepyceros melampus")

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Aepyceros_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.5,
                       value = "flagged")) #%>%
filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Aepyceros melampus",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Aepyceros_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(10, 25), ylim = c(-30, -15)) +
  
  ggtitle("Aepyceros melampus")

dev.off()

GBIF_Data <- bind_rows(GBIF_Data)
rm(Aepyceros_remove_IUCN)

### Antilocapra americana ####
Antilocapra_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Antilocapra americana", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Antilocapra americana",]
summary(Antilocapra_remove_IUCN)[["bbox"]]

GBIF_Base_Plots_02[["Antilocapra americana"]] +
  geom_polypath(data = fortify(Antilocapra_remove_IUCN ), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA,
                size = 0.15) + 
  coord_map(xlim = c(-160, -155), ylim = c(20, 25))

### Canis lupus ####

GBIF_Base_Plots_02[["Canis lupus"]]

# just need to remove domestic dogs, dingos, and subsp rufus (counted as separate species by iucn red list)
Canis_lupus_dat <- GBIF_Data %>%
  filter(species == "Canis lupus") %>%
  filter(!str_detect(verbatimScientificName, "familiaris|dingo|rufus")) %>%
  filter(infraspecificEpithet != "familiaris") %>%
  filter(scientificName != "familiaris") %>%   
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
                           TRUE                                               ~ "lupus")) %>%
  mutate(subsp = as.factor(subsp))

Canis_lupus_dat %>%             # Double checking for missed subspecies
  filter(subsp == "lupus") %>% 
  dplyr::select(verbatimScientificName, infraspecificEpithet, scientificName) %>% 
  mutate(across(c(1:3), as_factor)) %>% 
  summary()

base_map +
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Canis lupus", ]), 
                aes(x = long, y = lat, group = group),
                colour = "white",
                fill = "white") +
  
  geom_point(data = Canis_lupus_dat,
             aes(x = decimalLongitude, y = decimalLatitude, colour = subsp, shape = subsp)) +
  
  scale_shape_manual(values = c(18,0,17,1,16,2,3,1,4,5,8,15)) +
  ggtitle("Canis lupus")
beep(2)
rm(Canis_lupus_dat)

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Canis lupus")
GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  filter(!str_detect(verbatimScientificName, "familiaris|dingo|rufus")) %>%
  filter(infraspecificEpithet != "familiaris") %>%
  filter(scientificName != "familiaris")

base_map +
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Canis lupus", ]), 
                aes(x = long, y = lat, group = group),
                colour = "forestgreen",
                fill = NA) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  ggtitle("Canis lupus")
beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
gc()

#Although it hasn't removed all dubious samples, I can now fun species_cleaner with the rest.

### Cervus elaphus ####
# Not sure where to put this but saving it for later when justifying polygon merge

GBIF_Base_Plots_02[["Cervus elaphus"]]

base_map +
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Cervus elaphus", ]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen3",
               fill = "palegreen3") +
  
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Cervus canadensis", ]), 
               aes(x = long, y = lat, group = group),
               colour = "forestgreen",
               fill = "forestgreen") +
  
  geom_point(data = filter(GBIF_Data, species == "Cervus elaphus"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  geom_point(data = filter(GMPD_Raw_Data , HostCorrectedName == "Cervus elaphus"),
             aes(x = Longitude, y = Latitude),
             colour = "navy") +
  
  ggtitle("Cervus elaphus")
beep(2)

### Genetta genetta ####
Genetta_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Genetta genetta", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Genetta genetta",]

GBIF_Base_Plots_02[["Genetta genetta"]] +
  geom_polypath(data = fortify(Genetta_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = NA)
# Samples in Europe are introduced
# Should make a better map for this

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Genetta genetta")
GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Genetta_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 3,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Genetta genetta",]), 
                aes(x = long, y = lat, group = group),
                colour = "forestgreen",
                fill = NA) +
  
  geom_polypath(data = fortify(Genetta_remove_IUCN),
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = NA) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  ggtitle("Genetta genetta")
beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
rm(Genetta_remove_IUCN)
gc()

# Still a few left in Europe but they are caught by outlier fun

### Giraffa camelopardalis ####
Giraffa_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Giraffa camelopardalis", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Giraffa camelopardalis",]

GBIF_Base_Plots_02[["Giraffa camelopardalis"]] +
  geom_polypath(data = fortify(Giraffa_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = NA) + 
  coord_map(xlim = c(25, 45), ylim = c(-5, 5))

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Giraffa camelopardalis")
GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Giraffa_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.2,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Giraffa camelopardalis",]), 
                aes(x = long, y = lat, group = group),
                colour = "forestgreen",
                fill = NA) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_polypath(data = fortify(Giraffa_remove_IUCN),
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = NA) +
  
  coord_map(xlim = c(25, 45), ylim = c(-5, 5)) +
  
  ggtitle("Giraffa camelopardalis")
beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
rm(Genetta_remove_IUCN)
gc()

# Still need to use outlier fun to get ones in West 

### Meles meles ####
GBIF_Base_Plots_02[["Meles meles"]]

# Meles meles is the only one recorded in GMPD_Raw_Data
GMPD_Raw_Data %>% filter(str_detect(HostCorrectedName, "Meles")) %>% pull(HostCorrectedName) %>% unique()

# Recorded separately in GBIF
taxize::get_gbifid_(c("Meles meles", "Meles anakuma", "Meles leucurus"), method = "backbone") %>% 
  bind_rows()

# M anakuma, M. Leucurus
base_map +
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Meles anakuma", ]), 
               aes(x = long, y = lat, group = group),
               colour = "navy",
               fill = NA) +
  
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Meles leucurus", ]), 
               aes(x = long, y = lat, group = group),
               colour = "firebrick",
               fill = NA) +
  
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Meles meles", ]), 
               aes(x = long, y = lat, group = group),
               colour = "forestgreen",
               fill = NA) +
  
  geom_point(data = filter(GBIF_Data, species == "Meles meles"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_point(data = filter(GMPD_Raw_Data , HostCorrectedName == "Meles meles"),
             aes(x = Longitude, y = Latitude),
             colour = "orange") +
  
  ggtitle("Meles sp.")
beep(2)

Meles_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Meles leucurus", ] + 
  IUCN_Orders[IUCN_Orders$binomial == "Meles anakuma", ]
Meles_remove_IUCN$binomial <- "Meles meles"

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Meles meles")
GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Meles_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Meles meles", ]), 
               aes(x = long, y = lat, group = group),
               colour = "forestgreen",
               fill = NA) +
  
  geom_polypath(data = fortify(Meles_remove_IUCN),
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = NA) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  ggtitle("Meles meles")
beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
rm(Meles_remove_IUCN)
gc()

nrow(GBIF_Data) #2077658

#Although it hasn't removed all dubious samples, I can now fun species_cleaner with the rest.

### Mustela erminea ####
Mustela_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Mustela erminea", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela erminea",]

# Firstly there are non-native samples recorded in New Zealand
GBIF_Base_Plots_02[["Mustela erminea"]] +
  geom_polypath(data = fortify(Mustela_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA)

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Mustela erminea")
GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Mustela_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 1,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela erminea",]), 
                aes(x = long, y = lat, group = group),
                colour = "white",
                fill = "white") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  ggtitle("Mustela erminea")
beep(2)


# Also Mustela erminea can be confused with M. nivalis. 
# The M. erminea samples outside its range generally fall within the M nivalis range

base_map +
  geom_point(data = filter(GBIF_Data[["TRUE"]], species == "Mustela erminea"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela erminea",]), 
               aes(x = long, y = lat, group = group),
               colour = "navy",
               fill = NA) +
  
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela nivalis",]), 
               aes(x = long, y = lat, group = group),
               colour = "firebrick",
               fill = NA) +
  
  ggtitle("Mustela erminea")


# Using the same buffer method as CoordinateCleaner for consistency
Mustela_remove_IUCN <- rgeos::gBuffer(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela erminea",], byid = TRUE, width = 1)
Mustela_remove_IUCN <- IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela nivalis",] - Mustela_remove_IUCN 
Mustela_remove_IUCN$binomial <- "Mustela erminea"

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Mustela_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.5,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela nivalis",]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen4",
               fill = "palegreen4") +
  
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela erminea",]), 
               aes(x = long, y = lat, group = group),
               colour = "white",
               fill = "white") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  ggtitle("Mustela erminea")
beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
rm(Mustela_remove_IUCN)
gc()

# Doesn't catch the ones in America but fun will

### Mustela nivalis ####

Mustela_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Mustela nivalis", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela nivalis",]
GBIF_Base_Plots_02[["Mustela nivalis"]] +
  geom_polypath(data = fortify(Mustela_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) + 
  coord_map(xlim = c(-40, 40), ylim = c(25, 50))

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Mustela nivalis")

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Mustela_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.1,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Mustela nivalis",]), 
                aes(x = long, y = lat, group = group),
                colour = "white",
                fill = "white") +
  
  geom_polypath(data = fortify(Mustela_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  coord_map(xlim = c(-40, 40), ylim = c(25, 50)) +
  
  ggtitle("Mustela nivalis")
beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
rm(Mustela_remove_IUCN)
gc()

### Ovis canadensis ####
Ovis_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Ovis canadensis", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Ovis canadensis",]
GBIF_Base_Plots_02[["Ovis canadensis"]] +
  geom_polypath(data = fortify(Ovis_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(-130, -100), ylim = c(20, 55))

GBIF_Base_Plots_02[["Ovis canadensis"]] +
  geom_polypath(data = fortify(Ovis_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(-110, -100), ylim = c(32, 40))

GBIF_Base_Plots_02[["Ovis canadensis"]] +
  geom_polypath(data = fortify(Ovis_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(-118, -108), ylim = c(25, 33))

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Ovis canadensis")

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Ovis_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.2,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Ovis canadensis",]), 
                aes(x = long, y = lat, group = group),
                colour = "white",
                fill = "white") +
  
  geom_polypath(data = fortify(Ovis_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  coord_map(xlim = c(-110, -100), ylim = c(32, 40)) +
  
  ggtitle("Ovis canadensis")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Ovis canadensis",]), 
                aes(x = long, y = lat, group = group),
                colour = "white",
                fill = "white") +
  
  geom_polypath(data = fortify(Ovis_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick",
                fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  coord_map(xlim = c(-118, -108), ylim = c(25, 33)) +
  
  ggtitle("Ovis canadensis")

beep(2)

GBIF_Data <- bind_rows(GBIF_Data)
rm(Ovis_remove_IUCN)
gc()

### Rangifer tarandus ####
Rangifer_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Rangifer tarandus", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Rangifer tarandus",]
GBIF_Base_Plots_02[["Rangifer tarandus"]] +
  geom_polypath(data = fortify(Rangifer_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA)

# Neither cc_outl nor cc_iucn capture reindeer well
# cc_outl doesn't capture any
# cc_iucn excludes ones in Scandinavia which I would include
# if the buffer is expanded enough to include these it then captures those in central Europe

# All countrycodes are valid so I can use them
# Can also capture introduced records in Iceland and South Georgia
summary(as.factor(filter(GBIF_Data, species == "Rangifer tarandus")$countryCode))

countries_remove <- str_c(c("AUT", "BEL", "CHE", 
                            "CZE", "DEU", "FRA", 
                            "GBR", "HRV", "HUN", 
                            "IRL", "MDA", "NLD",
                            "POL", "ROU", "UKR", 
                            "ISL", "SGS"), collapse = "|")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Rangifer tarandus",]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen3",
               fill = "palegreen3") +
  
  geom_point(data = filter(GBIF_Data, species == "Rangifer tarandus"),
             aes(x = decimalLongitude, y = decimalLatitude, colour = str_detect(countryCode, countries_remove))) +
  
  ggtitle("Rangifer tarandus")

GBIF_Data <- GBIF_Data %>%
  filter(!(species == "Rangifer tarandus" & str_detect(countryCode, countries_remove)))

rm(countries_remove, Rangifer_remove_IUCN)

# happily gets all of them

### Rupicapra rupicapra ####
Rupicapra_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Rupicapra rupicapra", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Rupicapra rupicapra",]

pdf(file = here::here('GBIF Cleaning/Rupicapra rupicapra.pdf'), width = 10, height = 7)

GBIF_Base_Plots_02[["Rupicapra rupicapra"]] +
 geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
               aes(x = long, y = lat, group = group), 
               colour = "firebrick", fill = NA) +
  coord_map(xlim = c(-10, 50), ylim = c(30, 55))

GBIF_Base_Plots_02[["Rupicapra rupicapra"]] +
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(0, 7), ylim = c(42, 47))

GBIF_Base_Plots_02[["Rupicapra rupicapra"]] +
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(10, 21), ylim = c(46, 52))

GBIF_Base_Plots_02[["Rupicapra rupicapra"]] +
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(17, 23), ylim = c(37, 42))

GBIF_Base_Plots_02[["Rupicapra rupicapra"]] +
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA) +
  coord_map(xlim = c(12, 22), ylim = c(41, 46))


GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Rupicapra rupicapra")

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Rupicapra_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.2,
                       value = "flagged")) #%>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Rupicapra rupicapra",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(0, 7), ylim = c(42, 47)) +
  
  ggtitle("Rupicapra rupicapra")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Rupicapra rupicapra",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(10, 21), ylim = c(46, 52)) +
  
  ggtitle("Rupicapra rupicapra")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Rupicapra rupicapra",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(17, 23), ylim = c(37, 42)) +
  
  ggtitle("Rupicapra rupicapra")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Rupicapra rupicapra",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Rupicapra_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(12, 22), ylim = c(41, 46)) +
  
  ggtitle("Rupicapra rupicapra")

dev.off()

GBIF_Data <- bind_rows(GBIF_Data)
rm(Rupicapra_remove_IUCN)

### Syncerus caffer ####
Syncerus_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Syncerus caffer", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Syncerus caffer",]

GBIF_Base_Plots_02[["Syncerus caffer"]] +
  geom_polypath(data = fortify(Syncerus_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA,
                size = 0.15) +
  coord_map(xlim = c(-20, 50), ylim = c(-40, 20))

GBIF_Base_Plots_02[["Syncerus caffer"]] +
  geom_polypath(data = fortify(Syncerus_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA,
                size = 0.15) +
  coord_map(xlim = c(10, 20), ylim = c(-25, -15))

GBIF_Base_Plots_02[["Syncerus caffer"]] +
  geom_polypath(data = fortify(Syncerus_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA,
                size = 0.15) +
  coord_map(xlim = c(25, 35), ylim = c(-25, -35))

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Syncerus caffer")

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Syncerus_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0.2,
                       value = "flagged")) #%>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Syncerus caffer",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Syncerus_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick",
                size = 0.15) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(10, 20), ylim = c(-25, -15)) +
  
  ggtitle("Syncerus caffer")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Syncerus caffer",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Syncerus_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick",
                size = 0.15) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  coord_map(xlim = c(25, 35), ylim = c(-25, -35)) +
  
  ggtitle("Syncerus caffer")

GBIF_Data <- bind_rows(GBIF_Data)
rm(Syncerus_remove_IUCN)

### Ursus arctos ####

# polygon arithmetic wasn't working for this one. 
# Will sort later as no samples in removed polygon anyway
GBIF_Base_Plots_02[["Ursus arctos"]] + coord_map(xlim = c(30, 40), ylim = c(30, 35))

### Vulpes velox ####
GBIF_Base_Plots_02[["Vulpes velox"]]

## plotted IUCN polygon + GMPD_Data +  GBIF_Data
# Not much correspondance between IUCN and GBIF. 

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes velox",]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen3",
               fill = "palegreen3") +
  
  geom_point(data = filter(GBIF_Raw_Data, species == "Vulpes velox"),
             aes(x = decimalLongitude, y = decimalLatitude, colour = basisOfRecord)) +
  
  geom_point(data = filter(GBIF_Data, species == "Vulpes velox"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  ggtitle("Vulpes velox")

## Plotted GBIF_Raw_Data, aes(colour = basisOfRecord)
# Seems like basically all of the records are filtered out by _SPECIMEN so I will take a closer look at them

## Read wiki for Vulpes macrotis after googling "Vulpes velox fossil" which hinted they might be same
## Plotted V macrotis polygon with data
# Turns out the ones at the bottom are probably Vulpes macrotis (Kit fox) because some people treat them as one sp.

base_map +
  geom_point(data = filter(GBIF_Raw_Data, species == "Vulpes velox"),
             aes(x = decimalLongitude, y = decimalLatitude, colour = basisOfRecord)) +
  
  geom_point(data = filter(GBIF_Data, species == "Vulpes velox"),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes velox",]), 
               aes(x = long, y = lat, group = group),
               colour = "navy",
               fill = NA) +
  
  geom_polypath(data = fortify(IUCN_Orders[IUCN_Orders$binomial == "Vulpes macrotis", ]), 
               aes(x = long, y = lat, group = group),
               colour = "darkred",
               fill = NA) +
  
  ggtitle("Vulpes velox")

## Will remove some filters for Vulpes velox to regain data and then filter based on polygon overlap

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Vulpes velox")
GBIF_Data[["TRUE"]] <- GBIF_Raw_Data %>%
  filter(species == "Vulpes velox") %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.01 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE))) %>%
  filter(basisOfRecord != "UNKNOWN" & basisOfRecord != "FOSSIL_SPECIMEN")

GBIF_Data[["TRUE"]]$countryCode <- countrycode(GBIF_Data[["TRUE"]]$countryCode, origin = "iso2c", destination = "iso3c")

GBIF_Data[["TRUE"]] <- clean_coordinates(x = GBIF_Data[["TRUE"]],
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

Vulpes_remove_IUCN <- rgeos::gBuffer(IUCN_Orders[IUCN_Orders$binomial == "Vulpes macrotis", ], byid = TRUE, width = 2)
Vulpes_remove_IUCN <- Vulpes_remove_IUCN - rgeos::gBuffer(IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes velox",], byid = TRUE, width = 1)
Vulpes_remove_IUCN$binomial <- "Vulpes velox"

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Vulpes_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 0,
                       value = "flagged")) %>%
  filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes velox",]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen3",
               fill = "palegreen3") +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  ggtitle("Vulpes velox")

GBIF_Data <- bind_rows(GBIF_Data)
nrow(GBIF_Data) #2061904
rm(Vulpes_remove_IUCN)
gc()

### Vulpes vulpes ####
Vulpes_remove_IUCN <- IUCN_Orders[IUCN_Orders$binomial == "Vulpes vulpes", ] - IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes vulpes",]

GBIF_Base_Plots_02[["Vulpes vulpes"]] +
  geom_polypath(data = fortify(Vulpes_remove_IUCN), 
                aes(x = long, y = lat, group = group), 
                colour = "firebrick", fill = NA,
                size = 0.15)

GBIF_Data <- split(GBIF_Data, GBIF_Data$species == "Vulpes vulpes")

GBIF_Data[["TRUE"]] <- GBIF_Data[["TRUE"]] %>%
  mutate(out = cc_iucn(x = GBIF_Data[["TRUE"]],
                       range = Vulpes_remove_IUCN,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       species = "species",
                       buffer = 5,
                       value = "flagged")) #%>%
filter(!out) %>%
  dplyr::select(-out)

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes vulpes",]), 
                aes(x = long, y = lat, group = group),
                colour = "white", fill = "white",
                size = 0.15) +
  
  geom_polypath(data = fortify(Vulpes_remove_IUCN), 
                aes(x = long, y = lat, group = group),
                colour = "firebrick", fill = "firebrick",
                size = 0.15) +
  
  geom_point(data = GBIF_Data[["TRUE"]],
             aes(x = decimalLongitude, y = decimalLatitude, colour = out)) +
  
  ggtitle("Vulpes vulpes")

GBIF_Data <- bind_rows(GBIF_Data)
rm(Vulpes_remove_IUCN)

## Outliers ####

GBIF_Issues <- read.csv(here::here("GBIF cleaning/GBIF_issues.csv"), header = TRUE, stringsAsFactors = FALSE) 

GBIF_Data_Test <- apply(GBIF_Issues, MARGIN = 1, species_cleaner, dat = GBIF_Data)
# GBIF_Data_Test <- bind_rows(GBIF_Data_Test) # ran out of memory

# back up for now
saveRDS(GBIF_Data_Test, file = here::here("Data/Data back ups/GBIF_Data_Test"))

GBIF_Data_Test <- readRDS(here::here("Data/Data back ups/GBIF_Data_Test"))
names(GBIF_Data_Test) <- Host_Synonyms$GBIFName

GBIF_Outlier_Plots <- apply(GBIF_Issues, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data_Test, data_type = "outlier", range.polygon = IUCN_Native_Data)
names(GBIF_Outlier_Plots) <- Host_Synonyms$IUCNName

# Saving all outlier plots
pdf(file = here::here('GBIF cleaning/outliers.pdf'))
GBIF_Outlier_Plots
dev.off()

# Saving only ones I'm unsure of
pdf(file = here::here('GBIF cleaning/unsure_outliers.pdf'))
GBIF_Outlier_Plots[which(GBIF_Issues$finished == "Unsure", )] 
dev.off()

GBIF_Data_Cleaned <- GBIF_Data_Test %>%
  bind_rows() %>%
  filter(is.na(outlier) | outlier) %>%
  dplyr::select(-outlier)

nrow(GBIF_Data_Cleaned) #2053888

write.csv(GBIF_Data_Cleaned, file = here::here("Data/Data back ups/GBIF_Data.csv"), row.names = FALSE)

## Temporary plotting to get outlier parameters ####

# iucn test
Species <- c("Nyctereutes procyonoides",
             "Odocoileus hemionus",
             "Odocoileus virginianus",
             "Otocyon megalotis",
             "Ovis canadensis",
             "Panthera leo",
             "Panthera onca",
             "Panthera pardus",
             "Pecari tajacu",
             "Pelea capreolus",
             "Philantomba monticola",
             "Procyon lotor",
             "Puma concolor",
             "Rangifer tarandus",
             "Raphicerus campestris",
             "Redunca arundinum",
             "Redunca fulvorufula",
             "Rupicapra rupicapra",
             "Spilogale gracilis",
             "Sylvicapra grimmia",
             "Syncerus caffer",
             "Taxidea taxus",
             "Tragelaphus angasii",
             "Tragelaphus oryx",
             "Tragelaphus scriptus",
             "Tragelaphus spekii",
             "Tragelaphus strepsiceros",
             "Urocyon cinereoargenteus",
             "Urocyon littoralis",
             "Ursus americanus",
             "Ursus arctos",
             "Ursus maritimus",
             "Vulpes lagopus",
             "Vulpes velox",
             "Vulpes vulpes"
)

species_dat <- GBIF_Data %>%
  mutate(species = case_when(species == "Mustela vison"    ~ "Neovison vison",
                             species == "Taurotragus oryx" ~ "Tragelaphus oryx",
                             species == "Pekania pennanti" ~ "Martes pennanti",
                             TRUE                          ~ species)) %>%
  filter(species == Species[1]) 

species_dat$out <-  cc_iucn(x = rename(species_dat, binomial = species),
                            range = IUCN_Native_Data[IUCN_Native_Data$binomial == Species[1],],
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            species = "binomial",
                            buffer = 1,
                            value = "flagged")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == Species[1],]), 
               aes(x = long, y = lat, group = group),
               colour = "palegreen3",
               fill = "palegreen3") +
  
  geom_point(data = filter(species_dat, out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "navy") +
  
  geom_point(data = filter(species_dat, !out),
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "orange") +
  
  #coord_map(xlim = c(120, 150), ylim = c(25, 45))+
  
  ggtitle(Species[1])
beep(2)


#outlier test
Species <- "Giraffa camelopardalis"

species_dat <- filter(GBIF_Data, species == Species)

species_dat$out <- cc_outl(x = species_dat,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           method = "quantile",
                           mltpl = 10,
                           value = "flagged")

base_map +
  geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Giraffa camelopardalis",]), 
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
  
  ggtitle(Species)
beep(2)

