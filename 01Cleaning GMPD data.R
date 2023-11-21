# Reading and cleaning GMPD data by sample and location
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# takes:
## Raw GMPD data from https://doi.org/10.1002/ecy.1799
## Terrestrial mammal IUCN range polygons https://www.iucnredlist.org/resources/spatial-data-download

# makes: 
## GMPD_Data_res = cleaned GMPD data, restricted to points within native range polygons
## GMPD_Data_cln = cleaned GMPD_Data, cleaned to remove samples in non-native regions
## GMPD_Location_Data = associated location data for GMPD samples
## Native_DF = df indicating which polygons in a species range are native or non-native
## IUCN_Native_Data = SPDF of native IUCN ranges restricted according to Native_DF

# Libraries and data ##########################################################################################

library(beepr)
# library(CoordinateCleaner)  # cleaning geographic data  TODO: uses sp
# library(geosphere)          # calculating distances     TODO: uses sp
library(here)               #
library(maps)               # iso 3166 country codes and mapnames ## suggests sp but does not import
# library(maptools)             # TODO: depends sp
library(rnaturalearth)      # river data 
# library(rnaturalearthdata)  # countries data # TODO: depends sp
# library(raster)             # Use terra instead
# library(rgdal)              # read shapefiles
# library(rgeos)              # Just for gBuffer # TODO: retired
library(sf)
library(stringr)            #
library(tidyverse)          # beware of conflicts (mainly with raster)
library(tmap)
library(xlsx)

source(here::here("Functions.R"))

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), 
                          header = TRUE, 
                          stringsAsFactors = FALSE) 
nrow(GMPD_Raw_Data); length(unique(GMPD_Raw_Data$HostCorrectedName)) #Beginning with 24323 rows and 462 hosts

IUCN_Mammals <- sf::read_sf(dsn = here::here("Data/IUCN"), layer = "MAMMALS")
Projection_String <- sf::st_crs(IUCN_Mammals)

River_Data50 <- rnaturalearth::ne_load(scale = 50,
                                       type = "rivers_lake_centerlines",
                                       category = "physical",
                                       destdir = here::here("Data/Extras/ne_rivers"),
                                       returnclass = "sf")

# sf::sf_use_s2(FALSE) # For "invalid spherical geometry" errors
tmap_mode("view")

# Basic data cleaning #########################################################################################

## Removing/adjusting unuseable data ##########################################################################

# removing: 
## domestic species - may skew parasite measurements
## species with obfuscated location data (protected species) - inaccurate range position
## data with two reported species
## marine species - no associated abiotic data or likely poor relationship with abiotic data
## data that is unuseable in future analyses (missing variables)
## unique host-parasite pairs - uninformative in models

GMPD_Data <- GMPD_Raw_Data %>%
  dplyr::filter(HostCorrectedName != "Ovis aries" &
                  HostCorrectedName != "Bos frontalis" &
                  HostCorrectedName != "Bos grunniens" &
                  HostCorrectedName != "no binomial name" &
                  HostCorrectedName != "Dama dama" &
                  HostCorrectedName != "Diceros bicornis" &
                  HostCorrectedName != "Ceratotherium simum") %>%
  dplyr::filter(!str_detect(HostReportedName, " and ")) %>%
  dplyr::filter(!is.na(Prevalence)) %>%
  dplyr::filter(HostEnvironment != "marine") %>%
  dplyr::filter(NativeRange != "No" &
                  !is.na(NativeRange)) %>%
  dplyr::filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  dplyr::filter(!is.na(NumSamples) | !is.na(HostsSampled))
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 12060 rows of 200 hosts

GMPD_Data <- GMPD_Data %>%
  tidyr::separate_wider_delim(HostReportedName, 
                              delim = " ",
                              names = c("HostReportedGenus", 
                                        "HostReportedSpecies", 
                                        "HostReportedSubspecies"), 
                              cols_remove = FALSE,
                              too_few = "align_start",
                              too_many = "drop") 

# For recording information on subspecies decisions
# subsp_info <- GMPD_Data %>% 
#   select(HostReportedName, 
#          HostReportedGenus, 
#          HostReportedSpecies, 
#          HostReportedSubspecies, 
#          HostCorrectedName) %>% 
#   distinct()

# xlsx::write.xlsx(subsp_info, here::here("Data/Extras/Name Check.xlsx), row.names = FALSE, showNA = FALSE)

# Adding subspecies info and correcting HostCorrectedName according to IUCN
# Taxonomic justifications and citations in "Name Check.xlsx"
# Most recent update 05/07/2023
GMPD_Data <- GMPD_Data %>%
  dplyr::mutate(HostReportedSubspecies = case_when(HostReportedName  == "Alcelaphus cokii"           ~ "cokii", 
                                                   HostCorrectedName == "Alcelaphus lichtensteinii"  ~ "lichtensteinii", 
                                                   HostCorrectedName == "Axis axis"                  ~ "axis", 
                                                   
                                                   HostReportedName  == "Canis latrans Say"          ~ NA_character_, 
                                                   HostReportedName  == "Black-back Jackal"          ~ NA_character_, 
                                                   HostReportedName  == "Capra ibex ibex"            ~ NA_character_, 
                                                   HostReportedName  == "Capra i. ibex"              ~ NA_character_, 
                                                   HostReportedName  == "Cervus elaphus nelsoni"     ~ "canadensis", 
                                                   HostReportedName  == "Cervus elaphus hispanicus"  ~ "elaphus", 
                                                   HostReportedName  == "Cervus elaphus hippelaphus" ~ "elaphus", 
                                                   HostReportedName  == "Cervus nippon centralis"    ~ "nippon", 
                                                   
                                                   HostReportedName  == "Damaliscus korrigum"        ~ "korrigum", 
                                                   HostReportedName  == "Damaliscus dorcas dorcas"   ~ "pygargus", 
                                                   HostReportedName  == "Damaliscus pygargus dorcas" ~ "pygargus", 
                                                   HostCorrectedName == "Equus burchellii"           ~ "burchellii", 
                                                   str_detect(HostReportedName, "Felis libyca")      ~ "libyca", 
                                                   HostReportedName  == "Felis silvestris gordoni"   ~ "libyca", 
                                                   
                                                   HostReportedName  == "Giraffa reticulata"         ~ "reticulata", 
                                                   HostReportedName  == "Hyaena hyaena dubbah"       ~ NA_character_, 
                                                   HostReportedName  == "Kobus defassa"              ~ "defassa", 
                                                   HostReportedName  == "Lynx rufus floridanus"      ~ "rufus", 
                                                   
                                                   HostReportedName  == "Martes caurina"             ~ "caurina", 
                                                   HostReportedName  == "Meles meles anakuma"        ~ NA_character_, 
                                                   HostReportedName  == "Melogale moschata subauantiaca" ~ "subaurantiaca", 
                                                   HostReportedName  == "Mustela itatsi sho"         ~ NA_character_, 
                                                   HostReportedName  == "Neovison vison mink"        ~ NA_character_, 
                                                   
                                                   HostReportedName  == "Oryx gazella gazella"       ~ NA_character_, 
                                                   HostReportedName  == "Ourebia ourebi cottoni"     ~ NA_character_, 
                                                   HostReportedName  == "Ovibos moschatus moschatus" ~ NA_character_, 
                                                   HostReportedName  == "Ovibos moschatus wardi"     ~ NA_character_, 
                                                   HostReportedName  == "Ovis canadensis cremnobates" ~ "nelsoni", 
                                                   HostReportedName  == "Ovis canadensis mexicana"   ~ "nelsoni", 
                                                   
                                                   HostReportedName  == "Felis leo senegalensis"     ~ "leo", 
                                                   HostReportedName  == "Puma concolor coryi"        ~ "couguar", 
                                                   HostReportedName  == "Puma concolor stanleyana"   ~ "couguar", 
                                                   HostReportedName  == "Felis concolor coryi"       ~ "couguar", 
                                                   HostReportedName  == "Felis concolor vancouverensis" ~ "couguar", 
                                                   
                                                   HostReportedName  == "Spilogale gracilis amphiala" ~ "amphialus", 
                                                   HostReportedName  == "Urocyon cinereoargenteus texensis" ~ "scottii", 
                                                   HostReportedName  == "Ursus americanus pallas"    ~ NA_character_, 
                                                   HostReportedName  == "Ursus arctos marsicanus"    ~ "arctos", 
                                                   HostReportedName  == "Viverra civetta schwartzi"  ~ "schwarzi", 
                                                   HostReportedName  == "Vulpes fulva"               ~ "fulvus", 
                                                   HostReportedName  == "Vulpes vulpes schrencki"    ~ "schrenckii", 
                                                   
                                                   TRUE ~ HostReportedSubspecies)) %>%
  dplyr::mutate(HostCorrectedName = case_when(HostCorrectedName == "Alcelaphus lichtensteinii" ~ "Alcelaphus buselaphus", 
                                              HostCorrectedName == "Alces americanus"          ~ "Alces alces", 
                                              
                                              HostReportedName  == "Canis rufus"               ~ "Canis rufus", 
                                              HostReportedName  == "Cervus elaphus nannodes"   ~ "Cervus canadensis", 
                                              HostReportedName  == "Cervus elaphus nelsoni"    ~ "Cervus canadensis", 
                                              HostReportedName  == "Cervus elaphus roosevelti" ~ "Cervus canadensis", 
                                              HostReportedName  == "Cervus elaphus canadensis" ~ "Cervus canadensis", 
                                              
                                              HostCorrectedName == "Equus burchellii"          ~ "Equus quagga", 
                                              HostCorrectedName == "Felis manul"               ~ "Otocolobus manul", 
                                              str_detect(HostReportedName, "Felis libyca")     ~ "Felis libyca", 
                                              HostReportedName  == "Felis silvestris gordoni"  ~ "Felis libyca", 
                                              
                                              HostCorrectedName == "Lama glama"                ~ "Lama guanicoe", 
                                              HostCorrectedName == "Leopardus pajeros"         ~ "Leopardus colocolo", 
                                              HostReportedName  == "Meles meles anakuma"       ~ "Meles anakuma", 
                                              HostReportedName  == "Martes sibirica"           ~ "Msutela itatsi", 
                                              HostCorrectedName == "Neotragus moschatus"       ~ "Nesotragus moschatus", 
                                              
                                              HostReportedName  == "Putorius eversmanni"       ~ "Mustela eversmanni", 
                                              HostCorrectedName == "Puma yagouaroundi"         ~ "Herpailurus yagouaroundi", 
                                              HostCorrectedName == "Taurotragus oryx"          ~ "Tragelaphus oryx", 
                                              str_detect(HostReportedName, "Viverra civetta")  ~ "Civettictis civetta", 
                                              TRUE ~ HostCorrectedName)) %>%
  dplyr::filter(HostReportedName != "Ovis ammon musimon") 
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 12033 and 202

GMPD_Data <- GMPD_Data %>%
  dplyr::select(-ParasiteReportedName,
                -HostReportedName,
                -HostReportedGenus,
                -HostReportedSpecies,
                -HasBinomialName, 
                -NativeRange, 
                -Intensity, 
                -IntensityMeasure, 
                -SampleNotes) %>%   # no longer used
  dplyr::distinct()                                                                                                            #remove duplicated rows
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 11975 and 202

GMPD_Data <- GMPD_Data %>%
  dplyr::group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  dplyr::filter(n()>1) %>% 
  dplyr::ungroup()
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 10124 and 138

## adjust prevalence data that has been reported as a percentage

GMPD_Data <- GMPD_Data %>%
  dplyr::mutate(Prevalence = case_when(Prevalence >1 ~ Prevalence/100,
                                       TRUE          ~ Prevalence))

Hostlist <- unique(GMPD_Data$HostCorrectedName)
# GMPD_base_plots_00 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, plot_type = "base")
# names(GMPD_base_plots_00) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_base_plots_00.pdf'), width = 10, height = 7)
# GMPD_base_plots_00
# dev.off()

## Sample cleaning ############################################################################################

# pseudo-sampling and duplicates
## Throw samples which have nothing to use for sample size
## rows where SamplingBasis was "samples" AND HostsSampled was NA are excluded

## rows which are missing HostsSampled: get data from the other one

## Filling in HostAge and HostSex reduces duplicated data
### e.g. data inputted once with NA for HostSex and again identically except HostSex is reported accurately
### This doesn't get rid of rows where HostSex differs between rows but still contains a value

## Different sampling method on the same sample group

# TODO: reassess
# SamplingBasis and HostsSampled 
# Some samples had pseudoreplication when subdivided by HostAge, or HostSex, or both
# Before tackling those I assumed all with a single sample for each Host, Parasite and Location within a Citation was distinct
GMPD_Data_Test <- GMPD_Data %>% 
  dplyr::mutate(NewSamp = case_when(is.na(HostsSampled) & 
                                      SamplingBasis == "Animals" ~ NumSamples,
                                    TRUE ~ NA),
                NewPrev = case_when(is.na(HostsSampled) & 
                                      SamplingBasis == "Animals" ~ Prevalence,
                                    TRUE ~ NA)) %>% 
  dplyr::group_by(Citation, HostCorrectedName, ParasiteCorrectedName, LocationName) %>% 
  dplyr::mutate(NewSamp = case_when(n() == 1 & !is.na(HostsSampled) ~ HostsSampled,
                                    n() == 1 & is.na(HostsSampled) ~ NaN,
                                    TRUE ~ NewSamp),
                NewPrev = case_when(n() == 1 ~ Prevalence,
                                    TRUE ~ NewPrev))

# Some samples which lacked HostsSampled had valid sample sizes in NumSamples
# I could not distinguish which were valid just based on the available data so resorted to filtering through papers
# E.g. Some which recorded SamplingType as "DirectFecal" were sampled from individual hosts, while others were from faeces in the environment
# I only kept those which could be associated directly with the specified number of hosts
Citations_temp <- GMPD_Data_Test %>% filter(is.nan(NewSamp)) %>% pull(Citation) %>% unique()
GMPD_Data_Test %>% 
  filter(Citation %in% Citations_temp) %>% 
  View()

GMPD_Data_Test2 <- GMPD_Data_Test %>% 
  dplyr::group_by(Citation, HostCorrectedName, ParasiteCorrectedName, LocationName) %>% 
  mutate(NewSamp = case_when(Citation == "Alexander et al. 2010" ~ NumSamples,
                             
                             Citation == "Almberg et al. 2009" &
                               is.na(HostAge) ~ NumSamples,
                             Citation == "Almberg et al. 2009" &
                               NumSamples == max(NumSamples) ~ max(NumSamples),
                             Citation == "Almberg et al. 2009" ~ 99999,
                             
                             Citation == "Archer et al. 1986" ~ 99999,
                             
                             Citation == "Biek et al. 2006" ~ NumSamples, 
                             
                             Citation == "Calderini et al. 2009" ~ 129,
                             
                             Citation == "de Lisle et al. 2008" ~ NumSamples,
                             
                             Citation == "Evans 2002" ~ NumSamples,
                             
                             Citation == "Fuglei et al. 2008" ~ 99999,
                             
                             Citation == "Gompper et al. 2003" ~ NumSamples,
                             
                             Citation == "Gudmundsdottir and Skirnisson 2005" ~ 99999,
                             
                             Citation == "Hill et al. 1998" &
                               is.na(HostsSampled) ~ 99999,
                             
                             Citation == "Jenkins et al 2006" ~ 99999,
                             
                             Citation == "Magnarelli et al. 1995" ~ 99999,
                             
                             Citation == "Pedersen et al. 2008" &
                               is.na(HostsSampled) ~ 99999,
                             
                             Citation == "Popiolek et al. 2007" ~ 99999,
                             
                             Citation == "Rosalino et al. 2006" ~ 99999,
                             
                             Citation == "Szczesna and Popiolek 2007" ~ 99999,
                             
                             Citation == "Szczesna et al. 2008" ~ 99999,
                             
                             Citation == "Truyen et al. 1998" ~ NumSamples,
                             
                             Citation == "Tsukada et al. 2000" ~ 99999,
                             
                             Citation == "Whitlaw and Lankester 1994" ~ 99999, 
                             
                             all(is.na(HostAge)) & all(is.na(HostSex)) &
                               !is.na(HostsSampled) ~ HostsSampled,
                             all(is.na(HostAge)) & all(is.na(HostSex)) &
                               is.na(HostsSampled) ~ NaN,
                             
                             TRUE ~ NewSamp),
         NewPrev = case_when(Citation == "Almberg et al. 2009" &
                               is.na(HostAge) ~ Prevalence,
                             Citation == "Almberg et al. 2009" &
                               NumSamples == max(NumSamples) ~ Prevalence,
                             Citation == "Almberg et al. 2009" ~ NaN,
                             
                             Citation == "Archer et al. 1986" ~ NaN,
                             
                             Citation == "Biek et al. 2006" ~ Prevalence, 
                            
                             Citation == "Fuglei et al. 2008" ~ NaN, 
                             
                             Citation == "Gudmundsdottir and Skirnisson 2005" ~ NaN,
                             
                             Citation == "Hill et al. 1998" &
                               is.na(HostsSampled) ~ NaN,
                             
                             Citation == "Jenkins et al 2006" ~ NaN,
                             
                             Citation == "Magnarelli et al. 1995" ~ NaN, 
                             
                             Citation == "Pedersen et al. 2008" & 
                               is.na(HostsSampled) ~ NaN,
                             
                             Citation == "Popiolek et al. 2007" ~ NaN,
                             
                             Citation == "Rosalino et al. 2006" ~ NaN, 
                             
                             Citation == "Szczesna and Popiolek 2007" ~ NaN,
                             
                             Citation == "Szczesna et al. 2008" ~ NaN,
                             
                             Citation == "Tsukada et al. 2000" ~ NaN,
                             
                             Citation == "Whitlaw and Lankester 1994" ~ NaN,
                             
                             all(is.na(HostAge)) & all(is.na(HostSex)) ~ Prevalence,
                             
                             TRUE ~ NewPrev))

Citations_temp2 <- GMPD_Data_Test2 %>% filter(is.nan(NewSamp)) %>% pull(Citation) %>% unique()
GMPD_Data_Test2 %>% 
  filter(Citation %in% Citations_temp2) %>% 
  View()

# Another situation like previous sort through
# Then do a fix of host age/sex/both
# Then another sort through situation

# TODO: Host age and sex

# host age and sex
GMPD_Data_Test2 <- GMPD_Data_Test %>% 
  # select(Citation, SamplingBasis, Prevalence, HostsSampled, HostSex, HostAge, NumSamples, SamplingType) %>% 
  # filter(!is.na(NumSamples) & HostsSampled != NumSamples) %>% 
  group_by(Citation, HostCorrectedName, ParasiteCorrectedName, LocationName) %>% 
  mutate(NewPrev = case_when(#HostsSampled == NumSamples | is.na(NumSamples) ~ Prevalence,
                             
                             n() == 1 ~ Prevalence,
                             
                             all(is.na(HostAge)) & all(is.na(HostSex)) ~ Prevalence,
                             
                             all(!is.na(HostAge) | !is.na(HostSex)) &
                               (sum(NumSamples)/HostsSampled)%%1 == 0 ~ weighted.mean(Prevalence, NumSamples),
                             
                             TRUE ~ NA),
         NewSamp = case_when(#HostsSampled == NumSamples | is.na(NumSamples) ~ HostsSampled,
                             
                             n() == 1 & !is.na(HostsSampled) ~ HostsSampled,
                             
                             all(is.na(HostAge)) & all(is.na(HostSex)) ~ HostsSampled,
                             
                             all(!is.na(HostAge) | !is.na(HostSex)) &
                               (sum(NumSamples)/HostsSampled)%%1 == 0 ~ HostsSampled,
                             
                             TRUE ~ NA)) %>% 
  group_by(Citation, HostCorrectedName, ParasiteCorrectedName, LocationName, HostsSampled) %>% 
  mutate(NewPrev = case_when(!is.na(NewPrev) ~ NewPrev,
                             
                             all(!is.na(HostAge) | !is.na(HostSex)) &
                               (sum(NumSamples)/HostsSampled)%%1 == 0 ~ weighted.mean(Prevalence, NumSamples),
                             
                             TRUE ~ NewPrev),
         NewSamp = case_when(!is.na(NewSamp) ~ NewSamp,
                             
                             all(!is.na(HostAge) | !is.na(HostSex)) &
                               (sum(NumSamples)/HostsSampled)%%1 == 0 ~ HostsSampled,
                             TRUE ~ NewSamp)) #%>% 
  
  # mutate(NewAge = HostAge,
  #        NewSex = HostSex) %>% 
  # dplyr::group_by(pick(-NewAge, -HostAge)) %>% 
  # tidyr::fill(NewAge, .direction = "updown") %>% 
  # dplyr::ungroup() %>% 
  # 
  # dplyr::group_by(pick(-NewSex, -HostSex)) %>% 
  # tidyr::fill(NewSex, .direction  = "updown") %>% 
  # dplyr::ungroup() %>% 
  # 
  # dplyr::group_by(pick(-NewAge, -NewSex, -HostAge, -HostSex)) %>%
  # tidyr::fill(c(NewSex, NewAge), .direction = "updown") %>% 

  # group_split() %>% 

GMPD_Data_Test3 <- GMPD_Data_Test2 %>% 
  mutate(NumericSex = case_when(HostSex == "Female" ~ 1,
                                HostSex == "Male" ~ -1,
                                HostSex == "All" ~ 0,
                                is.na(HostSex) ~ NA),
         NumericAge = case_when(HostAge == "Adult" ~ 1,
                                HostAge == "Juvenile" ~ -1,
                                HostAge == "All" ~ 0,
                                is.na(HostAge) ~ NA)) %>% 
  group_by(Citation, HostCorrectedName, ParasiteCorrectedName, LocationName) %>% 
  mutate(NewNewSamp = case_when(sum(NumericSex) == 0 & HostsSampled == NumSamples ~ HostsSampled,
                                sum(NumericSex) == 0 ~ 99999,
                                sum(NumericAge) == 0 & HostsSampled == NumSamples ~ HostsSampled,
                                sum(NumericAge) == 0 ~ 99999),
         NewNewPrev = case_when(sum(NumericSex) == 0 & HostsSampled == NumSamples ~ Prevalence,
                                sum(NumericSex) == 0 ~ 99999,
                                sum(NumericAge) == 0 & HostsSampled == NumSamples ~ Prevalence,
                                sum(NumericAge) == 0 ~ 99999))

# host age and sex NA duplications
GMPD_Data <- GMPD_Data %>%
  dplyr::group_by(pick(-HostAge)) %>% 
  tidyr::fill(HostAge, .direction = "updown") %>% 
  
  dplyr::group_by(pick(-HostSex)) %>% 
  tidyr::fill(HostSex, .direction  = "updown") %>% 
  
  dplyr::group_by(pick(-HostAge, -HostSex)) %>%
  tidyr::fill(c(HostSex, HostAge), .direction = "updown") %>% 
  dplyr::distinct() %>%
  dplyr::ungroup()

# HostsSampled vs NumSamples
GMPD_Data <- GMPD_Data %>%
  dplyr::mutate(HostsSampled = case_when(is.na(HostsSampled) ~ NumSamples,
                                         TRUE                ~ HostsSampled)) %>%
  dplyr::mutate(NumSamples = case_when(is.na(NumSamples) ~ HostsSampled,
                                       TRUE              ~ NumSamples)) %>%
  dplyr::select(-NumSamples)
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 9561 and 138

# differing sample method on the same sample group

GMPD_Data <- GMPD_Data %>% 
  dplyr::group_by(pick(-SamplingType, -Prevalence)) %>%
  dplyr::mutate(sample_temp = case_when(n() == 1 ~ "fine",
                                        length(unique(SamplingType)) == 1 ~ "fine",
                                        length(unique(Prevalence))   == 1 & !duplicated(Prevalence) ~ "fine",
                                        length(unique(Prevalence))   == 1 ~ "not fine",
                                        n() == length(unique(SamplingType)) & Prevalence == max(Prevalence) ~ "fine",
                                        n() == length(unique(SamplingType)) ~ "not fine",
                                        length(unique(Prevalence)) > 1 & length(unique(SamplingType)) > 1 & Prevalence == max(Prevalence) ~ "fine",
                                        length(unique(Prevalence)) > 1 & length(unique(SamplingType)) > 1 ~ "not fine")) %>%
  dplyr::filter(sample_temp == "fine") %>%
  dplyr::select(-sample_temp) %>%
  dplyr::ungroup()
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 9514 138

# Cleaning by location ########################################################################################

## Fixing country names #######################################################################################

# for matching LocationName with iso3166 country names
## regex adjustments to iso3166 map names deal with incorrect matches
## not future-proof

Country_Match <- iso3166$mapname
Country_Match[which(Country_Match == "France")]   <- "France(?!s)"
Country_Match[which(Country_Match == "Niger")]    <- "Niger(?!i)"
Country_Match[which(Country_Match == "India")]    <- "India(?!n)"
Country_Match[which(Country_Match == "Mongolia")] <- "(?<!inner )Mongolia"
Country_Match[which(Country_Match == "Jersey")]   <- "(?<!New )Jersey"
Country_Match[which(Country_Match == "Mexico")]   <- "(?<!New )Mexico"
Country_Match[which(Country_Match == "USA")]      <- "(?<![:alpha:])USA"
Country_Match[which(Country_Match == "Kenya")]    <- "Kenya(?! border)"
Country_Match[which(Country_Match == "UK(?!r)")]  <- "(?<![:alpha:])UK(?![:alpha:])"
Country_Match[which(Country_Match == "Georgia")]  <- "(?<!.)Georgia(?!.)"
Country_Match[which(Country_Match == "Saba")]     <- "(?<![:alpha:])Saba"
Country_Match[which(Country_Match == "Jordan")]   <- "(?<![:alpha:])Jordan"
Country_Match[which(Country_Match == "Oman")]     <- "(?<![:alpha:])Oman"


Country_Match <- stringr::str_c(Country_Match, collapse = "|")
State_Match <- stringr::str_c(state.name, collapse = "|")

# All GMPD location descriptions needing matched to a country, 1972 unique descriptions
GMPD_Location_Data <- GMPD_Data %>%
  dplyr::select(LocationName) %>%
  dplyr::distinct()

# extracting any country names in descriptions
Country_Data_Temp <- dplyr::full_join(GMPD_Location_Data %>%
                                        tibble::rownames_to_column(),
                                      GMPD_Location_Data$LocationName %>%
                                        stringr::str_to_lower() %>%
                                        stringr::str_extract_all(stringr::str_to_lower(Country_Match), simplify = TRUE) %>%
                                        as.data.frame() %>%
                                        tibble::rownames_to_column(),
                                      by = "rowname")

# joining multiple countries into single strings
Country_Data_Temp <- Country_Data_Temp %>%
  dplyr::mutate(across(where(is.factor), as.character)) %>% 
  mutate_at(vars(-rowname, -LocationName), list(~ na_if(., ""))) %>% ## TODO: mutate_at -> mutate
  tidyr::unite("countries", names(Country_Data_Temp)[-c(1,2)], sep = "|", remove = TRUE, na.rm = TRUE)

# lacking country name in description or misspelled country name (while already having one correct country or state)
Country_Data_Temp <- Country_Data_Temp %>%
  dplyr::mutate(countries = case_when(LocationName == "Masai Mara National Park, Nairobi National Park, Ngorongoro crater,  Serengeti National Park and Namibia" ~ "namibia|kenya|tanzania",
                                      LocationName == "Czech Republic and Slovatkia"                   ~ "czech republic|slovakia",
                                      LocationName == "East and West Azarbaijan, Ardebil, Markazi, Isfahan, and Khorassan, IRAN" ~ "azerbaijan|iran",
                                      LocationName == "Macedonia, Thrace, Epirus, Peloponnesus, Thessaly, Sterea Hellas, Lesbos, and Lefka" ~ "greece",
                                      LocationName == "La Canada-Flintridge"                           ~ "usa",
                                      LocationName == "Lebanon, Pennsylvania"                          ~ "usa",
                                      LocationName == "San Marino, Los Angeles County, California"     ~ "usa",
                                      LocationName == "Slovak/Hungary border region (Dunajska Streda)" ~ "slovakia|hungary",
                                      LocationName == "Saint Martin-sous-Vigouroux, FRANCE"            ~ "france",
                                      LocationName == "Amama and Trinidad, Department of Moreno, Santiago del Estero" ~ "argentina",
                                      
                                      LocationName == "Glacier National Park, Montana and British Columbia" ~ "usa|canada",
                                      
                                      countries == "japan|japan"     ~ "japan",
                                      countries == "namibia|namibia" ~ "namibia",
                                      countries == "spain|spain"     ~ "spain", 
                                      countries == "bolivia|bolivia" ~ "bolivia",
                                      TRUE                           ~ countries))

# extracting US state names from descriptions
State_Data_Temp <- dplyr::full_join(GMPD_Location_Data %>%
                                      tibble::rownames_to_column(),
                                    GMPD_Location_Data$LocationName %>%
                                      stringr::str_to_lower() %>%
                                      stringr::str_extract_all(stringr::str_to_lower(State_Match),simplify = TRUE) %>%
                                      as.data.frame() %>%
                                      tibble::rownames_to_column(),
                                    by = "rowname")

# joining states into single string
State_Data_Temp <- State_Data_Temp %>%
  dplyr::mutate(across(where(is.factor), as.character)) %>% 
  mutate_at(vars(-rowname, -LocationName), list(~ na_if(., ""))) %>% ## TODO: mutate_at -> mutate
  tidyr::unite("states", names(State_Data_Temp)[-c(1,2)], sep = "|", remove = TRUE, na.rm = TRUE)

# merge country and state data, then correct missing names (grouped roughly by continent)
GMPD_Location_Data <- dplyr::full_join(State_Data_Temp, Country_Data_Temp, by = c("rowname", "LocationName")) %>%
  dplyr::mutate(countries = case_when(countries != "" ~ countries,
                                      states    != "" ~ "usa",
                                      stringr::str_detect(LocationName, regex("yellowstone|channel islands|eastern us| CA(?!.)|yosemite|nebrask|califonia|orange|purdue|marion|susitna|orgeon|tallahal", ignore_case = TRUE)) ~ "usa",
                                      stringr::str_detect(LocationName, regex("ontario|saskatchewan|quebec|yukon|nova scotia|alberta|prince edward|british columbia|northwest territories|baffin|vancou|newfoundland|brunswick|hudson bay", ignore_case = TRUE)) ~ "canada",
                                      stringr::str_detect(LocationName, regex("orkendalen",                       ignore_case = TRUE)) ~ "greenland",
                                      stringr::str_detect(LocationName, "Arctic")                                                      ~ "usa|canada",
                                      
                                      stringr::str_detect(LocationName, regex("parana|leones|marajo|jequitinhonha|pantanal", ignore_case = TRUE)) ~ "brazil",
                                      stringr::str_detect(LocationName, regex("mexican",                          ignore_case = TRUE)) ~ "mexico",
                                      stringr::str_detect(LocationName, regex("kaa",                              ignore_case = TRUE)) ~ "bolivia",
                                      
                                      stringr::str_detect(LocationName, regex("scotland|england|wales|united kingdom|great britain|shire|hebrides", ignore_case = TRUE)) ~ "uk",
                                      stringr::str_detect(LocationName, regex("brandenburg|berlin",               ignore_case = TRUE)) ~ "germany",
                                      stringr::str_detect(LocationName, regex("reykjavik",                        ignore_case = TRUE)) ~ "iceland",
                                      stringr::str_detect(LocationName, regex("noway|svalbard|barents",           ignore_case = TRUE)) ~ "norway",
                                      stringr::str_detect(LocationName, regex("bialowie|polish|Puszcza",          ignore_case = TRUE)) ~ "poland",
                                      stringr::str_detect(LocationName, regex("swiss",                            ignore_case = TRUE)) ~ "switzerland",
                                      stringr::str_detect(LocationName, regex("copenhagen",                       ignore_case = TRUE)) ~ "denmark",
                                      stringr::str_detect(LocationName, regex("zilina|slova",                     ignore_case = TRUE)) ~ "slovakia",
                                      stringr::str_detect(LocationName, regex("budakeszi",                        ignore_case = TRUE)) ~ "hungary",
                                      stringr::str_detect(LocationName, regex("moravia|mim",                      ignore_case = TRUE)) ~ "czech republic",
                                      stringr::str_detect(LocationName, regex("italia|sondrio|brembana|belviso",  ignore_case = TRUE)) ~ "italy",
                                      stringr::str_detect(LocationName, regex("meurthe|french|bauges|savoy",      ignore_case = TRUE)) ~ "france",
                                      stringr::str_detect(LocationName, regex("sorbe|zaragoza|malaga|catalonia|sierras|madrid|jaen|pallars|aller", ignore_case = TRUE)) ~ "spain", #assumed spain for jaen (as opposed to peru) because of species
                                      stringr::str_detect(LocationName, regex("dalmatia",                         ignore_case = TRUE)) ~ "croatia",
                                      stringr::str_detect(LocationName, regex("vojvodina",                        ignore_case = TRUE)) ~ "serbia",
                                      stringr::str_detect(LocationName, regex("danubian",                         ignore_case = TRUE)) ~ "romania",
                                      
                                      stringr::str_detect(LocationName, regex("alpine areas",                     ignore_case = TRUE)) ~ "italy|switzerland|france",
                                      stringr::str_detect(LocationName, regex("pyrenees",                         ignore_case = TRUE)) ~ "spain|france",
                                      
                                      stringr::str_detect(LocationName, regex("cameroun|ngaoun",                  ignore_case = TRUE)) ~ "cameroon",
                                      stringr::str_detect(LocationName, regex("kenia|masai|nairobi|jogi|bungoma", ignore_case = TRUE)) ~ "kenya",
                                      stringr::str_detect(LocationName, regex("kruger|natal|skukuza|queenstown|eastern shores|karroid|KNP|rooiwal|rietvlei|potchefstroom|benfontein|pieter|transvaa|sabi|hluhluwe|kuruman|mbiyamiti|ntomeni|west coast national park|transkei", ignore_case = TRUE)) ~ "south africa",
                                      stringr::str_detect(LocationName, regex("serengeti|ngorongoro|selous|temi|ruaha|kaisho", ignore_case = TRUE)) ~ "tanzania",
                                      stringr::str_detect(LocationName, regex("bale|sidamo|urso",                 ignore_case = TRUE)) ~ "ethiopia",
                                      stringr::str_detect(LocationName, regex("zimbawe|hippo|mana pools|buffalo range", ignore_case = TRUE)) ~ "zimbabwe",
                                      stringr::str_detect(LocationName, regex("zaire",                            ignore_case = TRUE)) ~ "democratic republic of the congo",
                                      stringr::str_detect(LocationName, regex("adiopodoume",                      ignore_case = TRUE)) ~ "ivory coast",
                                      stringr::str_detect(LocationName, regex("etosha",                           ignore_case = TRUE)) ~ "namibia",
                                      stringr::str_detect(LocationName, regex("ankole|koja|jie",                  ignore_case = TRUE)) ~ "uganda",
                                      stringr::str_detect(LocationName, regex("umsalala",                         ignore_case = TRUE)) ~ "sudan",
                                      stringr::str_detect(LocationName, regex("bandia|saboya",                    ignore_case = TRUE)) ~ "senegal",
                                      stringr::str_detect(LocationName, regex("batie",                            ignore_case = TRUE)) ~ "burkina faso",
                                      stringr::str_detect(LocationName, regex("ndoki|republic of the congo",      ignore_case = TRUE)) ~ "republic of congo",
                                      
                                      stringr::str_detect(LocationName, regex("kerguelen",                        ignore_case = TRUE)) ~ "french southern and antarctic lands",
                                      
                                      stringr::str_detect(LocationName, regex("new zeland|flagstaff",             ignore_case = TRUE)) ~ "new zealand",
                                      stringr::str_detect(LocationName, regex("tokyo|kkaido|sobo|shimane|akita",  ignore_case = TRUE)) ~ "japan",
                                      stringr::str_detect(LocationName, regex("rahasthan",                        ignore_case = TRUE)) ~ "india",
                                      stringr::str_detect(LocationName, regex("karak",                            ignore_case = TRUE)) ~ "jordan",
                                      TRUE ~ countries))

# split countries at "|" then pivot into single column 
# Number of rows will increase slightly because of location descriptions for which there are multiple countries. 
# They'll be removed again by coordinate cleaner
GMPD_Location_Data <- GMPD_Location_Data %>%
  tidyr::separate(countries, into = c("A","B","C"), sep = "\\|") %>% ## TODO: separate -> separate_woder_delim
  tidyr::pivot_longer(cols = c("A","B","C"), values_to = "mapname", values_drop_na = TRUE) %>%
  dplyr::select(-name, -states) %>%
  dplyr::mutate(mapname = case_when(mapname == "uk"      ~ "uk(?!r)",
                             mapname == "norway"  ~ "norway(?!:bouvet|:svalbard|:jan mayen)",
                             mapname == "finland" ~ "finland(?!:aland)",
                             mapname == "china"   ~ "china(?!:hong kong|:macao)",
                             TRUE                 ~ mapname)) # ignore warning

GMPD_Location_Data <- iso3166 %>%
  dplyr::select(a3, mapname) %>%
  dplyr::mutate(mapname = tolower(mapname)) %>%
  dplyr::right_join(GMPD_Location_Data, by = "mapname") %>%
  dplyr::rename(countrycode = a3)

GMPD_Data <- dplyr::full_join(GMPD_Location_Data, GMPD_Data, by = "LocationName", relationship = "many-to-many")

rm(State_Match, Country_Match)
rm(list = ls(pattern = "_Temp$"))

## CoordinateCleaner tests ####################################################################################

# https://peerj.com/articles/9916.pdf
GMPD_Data <- GMPD_Data %>%
  CoordinateCleaner::clean_coordinates(lon = "Longitude",
                                       lat = "Latitude",
                                       species = "HostCorrectedName",
                                       countries = "countrycode",
                                       tests = c("capitals","centroids","institutions", "countries"),
                                       value = "clean")
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 8832 and 132

# Native and non-native ####
## Removing non-native polygons ####

# reducing to relevant subsets of IUCN data to save time and space 
# TODO: Maybe move hostlist to here??
IUCN_Data <- IUCN_Mammals[IUCN_Mammals$binomial %in% Hostlist, ]

# Correcting legend for Oreamnos americanus

### Oreamnos americanus
O_americanus_temp <- IUCN_Data[IUCN_Data$binomial == "Oreamnos americanus",]

# According to iucn red list, O. americanus was introduced to Chicagof Island
"Festa-Bianchet, M. 2022. Oreamnos americanus (errata version published in 2022). The IUCN Red List of Threatened Species 2022: e.T42680A211860282. Accessed on 27 September 2023."
# It is recorded in their spatial data as Extant (resident)
# Correcting to Extant & Introduced (resident)
# The same applies to Kodiak island
tm_shape(O_americanus_temp) + tm_polygons("legend") 

# Chicagof and Kodiak are missing data in dist_comm
tm_shape(O_americanus_temp) + tm_polygons("dist_comm") 

O_americanus_temp <- O_americanus_temp %>% 
  dplyr::mutate(legend = case_when(is.na(dist_comm) ~ "Extant & Introduced (resident)",
                                   TRUE ~ legend))

# Corrected map
tm_shape(O_americanus_temp) + tm_polygons("legend")

IUCN_Data <- IUCN_Data %>% 
  dplyr::filter(binomial != "Oreamnos americanus") %>% 
  rbind(O_americanus_temp)

### Ovibos moschatus
# O_moschatus_temp <- IUCN_Data[IUCN_Data$binomial == "Ovibos moschatus",]

# IUCN red list states:
# This has been corrected, I need to update my data
# Can state in reference that this one has a more recent download date because 
# the polygons were updataed by IUCN in a way that changed the legend
# "Muskoxen were introduced and are well established in West Greenland"

# tm_shape(O_moschatus_temp) + tm_polygons("legend")
# tm_shape(O_moschatus_temp) + tm_polygons("SHAPE_Area") # aggregated across Greenland
# 
# clip_points_temp <- as.data.frame(matrix(c(-54.86, 71.54,
#                                            -49.98, 66.64,
#                                            -47.69, 61.36),
#                                          ncol = 2, byrow = TRUE))
# clip_points_temp <- st_as_sf(clip_points_temp,
#                              crs = Projection_String,
#                              coords = c(1,2))
# 
# O_moschatus_temp <- st_cast(O_moschatus_temp)
# O_moschatus_temp <- O_moschatus_temp[clip_points_temp,]
# O_moschatus_temp$legend <- "Extant & Introduced (resident)"
# 
# O_moschatus_temp <- rbind(st_difference(IUCN_Data[IUCN_Data$binomial == "Ovibos moschatus",], O_moschatus_temp$geometry),
#                           O_moschatus_temp)
# 
# tm_shape(O_moschatus_temp) + tm_polygons("legend") + tm_shape(clip_points_temp) + tm_dots()
# 
# IUCN_Data <- IUCN_Data[IUCN_Data$binomial != "Ovibos moschatus",]
# IUCN_Data <- rbind(IUCN_Data, O_moschatus_temp)

# Rupicapra rupicapra
# IUCN red list:
# "The subspecies cartusiana is endemic to France, where it is restricted to a 350 km2 area of the Chartreuse limestone massif, centred around Grenoble, at the western edge of the French Alps."

# R_rupicapra_temp <- IUCN_Data[IUCN_Data$binomial == "Rupicapra rupicapra",]
# tm_shape(R_rupicapra_temp) + tm_polygons("legend")
# 
# clip_points_temp <- as.data.frame(matrix(c(2.935324, 45.199909),
#                                          ncol = 2, byrow = TRUE))
# clip_points_temp <- st_as_sf(clip_points_temp, 
#                              crs = Projection_String,
#                              coords = c(1,2))
# 
# R_rupicapra_temp <- sf::st_cast(R_rupicapra_temp)
# R_rupicapra_temp <- R_rupicapra_temp[clip_points_temp,]
# R_rupicapra_temp$legend <- "Extant & Reintroduced (Extant)"
# 
# R_rupicapra_temp <- rbind(sf::st_difference(IUCN_Data[IUCN_Data$binomial == "Rupicapra rupicapra",], 
#                                             R_rupicapra_temp$geometry), 
#                           R_rupicapra_temp)
# 
# tm_shape(R_rupicapra_temp) + tm_polygons("legend")
# 
# IUCN_Data <- IUCN_Data[IUCN_Data$binomial != "Rupicapra rupicapra",]
# IUCN_Data <- rbind(IUCN_Data, R_rupicapra_temp)
# 
# rm(list = ls(pattern = "_temp$"))

# Plotting with full polygons before restricting

Legend_Text <- sort(unique(IUCN_Data$legend))
# GMPD_plots_01 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Data, plot_type = "iucn")
# names(GMPD_plots_01) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_01.pdf'), width = 10, height = 7)
# GMPD_plots_01
# dev.off()

# removing unwanted polygons
# Justifications for exceptions based on IUCN Geographic Range descriptions
# Cervus elaphus
"In Greece, the small isolated subpopulations are the result of reintroductions 
into areas where it previously occurred. The last native population of Greek Red 
Deer is supposed to have survived in the Sithonia peninsula (Chalkidiki, 
north-eastern Greece) where it became extinct in the 1980s (Masseti 2012 and 
references therein)."
# Lovari, S., Lorenzini, R., Masseti, M., Pereladova, O., Carden, R.F., Brook, S.M. & Mattioli, S. 2018. Cervus elaphus (errata version published in 2019). The IUCN Red List of Threatened Species 2018: e.T55997072A142404453. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T55997072A142404453.en. Accessed on 06 July 2023.

# Cervus nippon
"Specifically, it was originally found in China (formerly from Manchuria south 
to Guangxi, and Sichuan to Anhui), North and South Korea (including Cheju 
Island) (but now probably extinct in both countries), Japan, Russia (a few 
places in Primorsky in the Far East), Taiwan (extinct in 1969, but subsequently 
re-introduced), and Viet Nam (probably now extinct)."
# Harris, R.B. 2015. Cervus nippon. The IUCN Red List of Threatened Species 2015: e.T41788A22155877. https://dx.doi.org/10.2305/IUCN.UK.2015-2.RLTS.T41788A22155877.en. Accessed on 06 July 2023.

# Lynx canadensis
# Area with legend "Presence Uncertain & Origin Uncertain" is adjacent to "Extant (resident)" and "Extant & Vagrant (seasonality uncertain)"
# Latitudinal span of this area is also within the bounds of the "Extant (resident)" range

# Vashon, J. 2016. Lynx canadensis. The IUCN Red List of Threatened Species 2016: e.T12518A101138963. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T12518A101138963.en. Accessed on 06 July 2023.

Native_DF <- IUCN_Data %>% 
  dplyr::select(binomial, legend) %>% 
  sf::st_drop_geometry() %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(keep = case_when(str_detect(binomial, "Cervus") & legend == "Extant & Introduced (resident)"  	    ~ TRUE,
                          binomial == "Lynx canadensis" & legend == "Presence Uncertain & Origin Uncertain" ~ TRUE,
                          str_detect(legend, "Introduced") 								                                  ~ FALSE,
                          legend == "Extant & Origin Uncertain (resident)"						                      ~ FALSE,
                          TRUE												                                                      ~ TRUE)) 

IUCN_Native_Data <- Native_DF %>% 
  dplyr::filter(keep) %>% 
  dplyr::select(-keep) %>% 
  {dplyr::right_join(IUCN_Data, ., by = c("binomial", "legend"), relationship = "many-to-one")}


# GMPD_plots_native_base_02 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_base_02) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_base_02.pdf'), width = 10, height = 7)
# GMPD_plots_native_base_02
# dev.off()

## Adjusting subspecies ####
# https://github.com/r-spatial/sf/wiki/Migrating
source(here::here("01.1Subspecies polygons.R"))

## Restricting by native IUCN polygon by species ####
Hostlist <- unique(GMPD_Spatial@data$HostCorrectedName)
GMPD_Data_res_all <- lapply(Hostlist, function(host, dat = GMPD_Spatial, range.polygon = IUCN_Native_Data, buff = 0, subsp = FALSE){
  range.polygon <- range.polygon[range.polygon$binomial == host,]
  out <- pip_test(host, dat, range.polygon, buff, subsp)
})

GMPD_Data_res_all <- GMPD_Data_res_all[which(lapply(GMPD_Data_res_all, function(x)nrow(x) != 0) == TRUE)]
GMPD_Data_res_all <- do.call(raster::bind, GMPD_Data_res_all)
GMPD_Data_res_all <- as.data.frame(GMPD_Data_res_all)

GMPD_Data_res_all <- GMPD_Data_res_all %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data_res_all); length(unique(GMPD_Data_res_all$HostCorrectedName)) # 7116 and 106

# GMPD_plots_native_restricted_03 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_res_all, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_restricted_03) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_restricted_03.pdf'), width = 10, height = 7)
# GMPD_plots_native_restricted_03
# dev.off()

## Restricting by native IUCN polygon by subspecies group ####

buffer_temp <- unique(GMPD_Subgroups[,c("HostCorrectedName", "subgroup")])
buffer_temp$buff <- 0

GMPD_Data_res_sub <- lapply(Hostlist, function(host, dat = GMPD_Spatial, range.polygon = IUCN_Native_Data, buff = buffer_temp){
  range.polygon <- range.polygon[range.polygon$binomial == host,]
  buff <- buff[buff$HostCorrectedName == host,]
  out <- pip_test(host, dat, range.polygon, buff)
})

GMPD_Data_res_sub <- GMPD_Data_res_sub[which(lapply(GMPD_Data_res_sub, is.null) == FALSE)]
GMPD_Data_res_sub <- do.call(raster::bind, GMPD_Data_res_sub)
GMPD_Data_res_sub <- as.data.frame(GMPD_Data_res_sub)

GMPD_Data_res_sub <- GMPD_Data_res_sub %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data_res_sub); length(unique(GMPD_Data_res_sub$HostCorrectedName)) # 7113 and 105

## Cleaning by native IUCN polygon ####
# Method in 01.1 already cleans species and subspecies

GMPD_Data_cln <- GMPD_Subgroups %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data_cln); length(unique(GMPD_Data_cln$HostCorrectedName)) # 8130 and 113

# GMPD_plots_native_clean_03 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_cln, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_clean_03) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_clean_03.pdf'), width = 10, height = 7)
# GMPD_plots_native_clean_03
# dev.off()

# Restricting each by proximity ####################################################################################
## IUCN restricted subgroup ####
Res_Temp <- restrict_deci(GMPD_Data_res_sub, subsp = TRUE)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_res_sub <- GMPD_Data_res_sub[GMPD_Data_res_sub$HostCorrectedName %in% Res_Temp$HostCorrectedName &
                                         GMPD_Data_res_sub$subgroup %in% Res_Temp$subgroup,]
GMPD_Data_res_sub <- GMPD_Data_res_sub %>%
  group_by(HostCorrectedName, subgroup) %>%
  filter(n() > 1) %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  filter(n() > 1) %>%
  ungroup() 
nrow(GMPD_Data_res_sub); length(unique(GMPD_Data_res_sub$HostCorrectedName)) # 6958 and 90

## IUCN restricted species ####
Res_Temp <- restrict_deci(GMPD_Data_res_all)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_res_all <- GMPD_Data_res_all[GMPD_Data_res_all$HostCorrectedName %in% Res_Temp$HostCorrectedName,]
GMPD_Data_res_all <- GMPD_Data_res_all %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  filter(n() > 1) %>%
  ungroup() 
nrow(GMPD_Data_res_all); length(unique(GMPD_Data_res_all$HostCorrectedName)) # 7054 and 94

## IUCN cleaned subgroup ####
Res_Temp <- restrict_deci(GMPD_Data_cln, subsp = TRUE)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_cln_sub <- GMPD_Data_cln[GMPD_Data_cln$HostCorrectedName %in% Res_Temp$HostCorrectedName &
                                 GMPD_Data_cln$subgroup %in% Res_Temp$subgroup,]
GMPD_Data_cln_sub <- GMPD_Data_cln_sub %>%
  group_by(HostCorrectedName, subgroup) %>%
  filter(n() > 1) %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  filter(n() > 1) %>%
  ungroup()
nrow(GMPD_Data_cln_sub); length(unique(GMPD_Data_cln_sub$HostCorrectedName)) # 7995 and 107

## IUCN cleaned species ####
Res_Temp <- restrict_deci(GMPD_Data_cln)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_cln_all <- GMPD_Data_cln[GMPD_Data_cln$HostCorrectedName %in% Res_Temp$HostCorrectedName,]
GMPD_Data_cln_all <- GMPD_Data_cln_all %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  filter(n() > 1) %>%
  ungroup()
nrow(GMPD_Data_cln_all); length(unique(GMPD_Data_cln_all$HostCorrectedName)) # 8110 and 110

rm(Res_Temp)

# TODO: sort out afterthought
# merging gmpd into one (afterthought) ####
GMPD_Data_cln_all$CleanAll <- TRUE
GMPD_Data_cln_sub$CleanSub <- TRUE
GMPD_Data_res_all$RestrAll <- TRUE
GMPD_Data_res_sub$RestrSub <- TRUE

GMPD_Data <- merge(GMPD_Data_cln_all, GMPD_Data_cln_sub, 
                   by = intersect(names(GMPD_Data_cln_all), names(GMPD_Data_cln_sub)),
                   all = TRUE)

GMPD_Data <- merge(GMPD_Data, GMPD_Data_res_all, 
                   by = intersect(names(GMPD_Data), names(GMPD_Data_res_all)),
                   all = TRUE)

GMPD_Data <- merge(GMPD_Data, GMPD_Data_res_sub, 
                   by = intersect(names(GMPD_Data), names(GMPD_Data_res_sub)),
                   all = TRUE)

GMPD_Data <- GMPD_Data %>%
  mutate(CleanSub = case_when(is.na(CleanSub) ~  FALSE, TRUE ~ CleanSub)) %>%
  mutate(RestrAll = case_when(is.na(RestrAll) ~  FALSE, TRUE ~ RestrAll)) %>%
  mutate(RestrSub = case_when(is.na(RestrSub) ~  FALSE, TRUE ~ RestrSub))

# Write files #################################################################################################
# narrowing down IUCN_Mammals to a more manageable size
IUCN_Orders <- IUCN_Mammals[IUCN_Mammals$order_ %in% unique(IUCN_Native_Data$order_), ]
rm(IUCN_Mammals)

write.csv(GMPD_Data, file = here::here("Data/Data back ups/GMPD_Data_01.csv"), row.names = FALSE)

write.csv(GMPD_Location_Data, file = here::here("Data/Data back ups/GMPD_Location_Data_01.csv"), row.names = FALSE)
write.csv(Native_DF, file = here::here("Data/Data back ups/Native_DF_01.csv"), row.names = FALSE)

saveRDS(IUCN_Native_Data, file = here::here("Data/Data back ups/IUCN_Native_Data_01")) # too large to commit 
saveRDS(IUCN_Orders, file = here::here("Data/Data back ups/IUCN_Orders_01")) # too large to commit 

