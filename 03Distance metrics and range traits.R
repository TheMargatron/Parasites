# Calculating distance metrics (two methods), extracting range traits
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ####

library(here)               #
library(tidyverse)          #

source(here::here("Functions.R"))

GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data_01.csv"), header = TRUE, stringsAsFactors = FALSE)

#GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Data.csv"), header = TRUE, stringsAsFactors = FALSE)

GBIF_Subgroups <- read.csv(here::here("Data/Data back ups/GBIF_Subgroups_02.csv"), header = TRUE, stringsAsFactors = FALSE)

IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))

# Prep ####
# simplifying data to essentials

GMPD_Data_Temp <- GMPD_Data %>%
  dplyr::select(HostCorrectedName, Group,
                HostOrder, HostFamily,
                
                ParasiteCorrectedName, ParType, 
                ParPhylum, ParClass,
                
                Citation, 
                LocationName, PopulationType,
                Longitude, Latitude, 
                 
                SamplingBasis, SamplingType, HostsSampled, 
                Prevalence, 
                HostSex, HostAge, 
                subgroup, 
                
                CleanAll, CleanSub,
                RestrAll, RestrSub)

GBIF_Data_Temp <- GBIF_Subgroups %>%
  dplyr::select(species, taxonRank, scientificName,
                subgroup, decimalLongitude, decimalLatitude) %>%
  dplyr::mutate(species = case_when(species == "Pekania pennanti" ~ "Martes pennanti", # TODO: remove after rerunning 02.1
                                    species == "Mustela vison" ~ "Neovison vison",
                                    species == "Taurotragus oryx" ~ "Tragelaphus oryx",
                                    TRUE ~ species))

# Distance metrics and range traits ####
Distances_Data_res_all <- range_distances(dat = GMPD_Data_Temp[GMPD_Data_Temp$RestrAll,], range.pol = IUCN_Native_Data, method = "iucn", subsp = FALSE)

Distances_Data_res_sub <- range_distances(dat = GMPD_Data_Temp[GMPD_Data_Temp$RestrSub,], range.pol = IUCN_Native_Data, method = "iucn", subsp = TRUE) 

Distances_Data_cln_all <- range_distances(dat = GMPD_Data_Temp[GMPD_Data_Temp$CleanAll,], range.dat = GBIF_Data_Temp, method = "gbif", subsp = FALSE)

Distances_Data_cln_sub <- range_distances(dat = GMPD_Data_Temp[GMPD_Data_Temp$CleanSub,], range.dat = GBIF_Data_Temp, method = "gbif", subsp = TRUE) 

# Merge into one df
GMPD_Distances_Data <- dplyr::bind_rows(Distances_Data_res_all$DistanceMetrics,
                                        Distances_Data_res_sub$DistanceMetrics,
                                        Distances_Data_cln_all$DistanceMetrics,
                                        Distances_Data_cln_sub$DistanceMetrics)

Range_Traits <- dplyr::bind_rows(Distances_Data_res_all$RangeTraits,
                                 Distances_Data_res_sub$RangeTraits,
                                 Distances_Data_cln_all$RangeTraits,
                                 Distances_Data_cln_sub$RangeTraits)


write.csv(GMPD_Distances_Data, file = here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), row.names = FALSE)
write.csv(Range_Traits, file = here::here("Data/Data back ups/Range_Traits_03.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)
