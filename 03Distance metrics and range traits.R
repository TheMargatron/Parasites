# Calculating distance metrics (two methods), extracting range traits
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ####

library(here)               #
library(tidyverse)          #
library(sf)

source(here::here("Functions.R"))

GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data_01.csv"), header = TRUE, stringsAsFactors = FALSE)

GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Data_02.csv"), header = TRUE, stringsAsFactors = FALSE)

IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))

Projection_String <- sf::st_crs(IUCN_Native_Data)
sf::sf_use_s2(FALSE) # For "invalid spherical geometry" errors


# Prep ####
# simplifying data to essentials

# GMPD_Data <- GMPD_Data %>%
#   dplyr::select(-countrycode,
#                 -mapname,
#                 -rowname,
#                 -LocationName,
#                 -HostReportedSubspecies,
#                 -HostEnvironment,
#                 -PopulationType)

GBIF_Data <- GBIF_Data %>%
  dplyr::select(species, decimalLongitude, decimalLatitude) 

# Distance metrics and range traits ####
# IUCN
Distances_Data_res_all <- lapply(unique(GMPD_Data[GMPD_Data$RestrAll, "HostCorrectedName"]),
                                 range_distances_host,
                                 dat = GMPD_Data[GMPD_Data$RestrAll,], 
                                 range.object = IUCN_Native_Data, 
                                 method = "iucn")

Distances_Data_res_all_dat <- lapply(Distances_Data_res_all, function(x){x[[1]]}) %>% bind_rows() 
  
Distances_Data_res_all_range <- lapply(Distances_Data_res_all, function(x){x[[2]]}) %>% bind_rows()

# GBIF
Distances_Data_cln_all <- lapply(unique(GMPD_Data[GMPD_Data$CleanAll, "HostCorrectedName"]),
                                 range_distances_host,
                                 dat = GMPD_Data[GMPD_Data$CleanAll,],
                                 range.object = GBIF_Data, 
                                 method = "gbif")

GMPD_Distances_Data <- lapply(Distances_Data_cln_all, function(x){x[[1]]}) %>% 
  bind_rows() %>% 
  full_join(Distances_Data_res_all_dat, 
            by = names(.)[!grepl("Prop|Dist|Above", names(.), ignore.case = TRUE)],
            suffix = c("_gbif", "_iucn")) %>% 
  mutate(CleanAll = case_when(is.na(AboveMedn_gbif) ~ FALSE,
                              TRUE ~ CleanAll))

Range_Traits <- lapply(Distances_Data_cln_all, function(x){x[[2]]}) %>% 
  bind_rows() %>% 
  bind_rows(Distances_Data_res_all_range)

write.csv(GMPD_Distances_Data, file = here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), row.names = FALSE)
write.csv(Range_Traits, file = here::here("Data/Data back ups/Range_Traits_03.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)
