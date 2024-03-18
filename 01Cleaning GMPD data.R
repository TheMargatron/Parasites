# Reading and cleaning GMPD data by sample and location
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# takes:
## Raw GMPD data from https://doi.org/10.1002/ecy.1799 (Downloaded 2018)

## Terrestrial mammal IUCN range polygons https://www.iucnredlist.org/resources/spatial-data-download
# 07/03/2024

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
# library(rnaturalearth)      # river data 
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

IUCN_Mammals <- sf::read_sf(dsn = here::here("Data/IUCN/MAMMALS"), layer = "MAMMALS")
Projection_String <- sf::st_crs(IUCN_Mammals)

# River_Data50 <- rnaturalearth::ne_load(scale = 50,
#                                        type = "rivers_lake_centerlines",
#                                        category = "physical",
#                                        destdir = here::here("Data/Extras/ne_rivers"),
#                                        returnclass = "sf")

sf::sf_use_s2(FALSE) # For "invalid spherical geometry" errors
# tmap_mode("view")

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
  dplyr::filter(!is.na(NumSamples) | !is.na(HostsSampled)) %>% 
  dplyr::filter(!(HostsSampled == 0 & NumSamples == 0))

nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 12054 rows of 200 hosts

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
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 12027 and 202

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
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 11969 and 202

GMPD_Data <- GMPD_Data %>%
  dplyr::group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  dplyr::filter(n()>1) %>% 
  dplyr::ungroup()
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 10119 and 138

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
## Remove samples where animals were selected based on body condition
## Remove samples where feces was collected from an unknown number of hosts (e.g. environmental samples)
## Correct or remove rows which repeat data
## Remove samples from semi-domestic or captive populations
## Remove lower Prevalence estimates when different sampling method on the same sample group

## rows which are missing HostsSampled: get data from NumSamples


# Starting with DirectFecal because they are often sampled from the environment
Citations <- GMPD_Data %>% filter(SamplingType == "DirectFecal") %>% pull(Citation)
# GMPD_Data %>% filter(Citation %in% Citations) %>% View()

GMPD_Data_Filter <- GMPD_Data %>% 
  dplyr::mutate(SampleSize = case_when(
    # Kept
    Citation %in% c("Aramini et al. 1998",
                    "Aschfalk et al 2008",
                    "Bolt et al. 1992",
                    "Criffield et al. 2009",
                    "Davidson et al. 2006b",
                    "Deem and Emmons 2005",
                    "Dubinsky et al. 1999",
                    "Fiorello et al. 2006",
                    "Hamnes et al 2006",
                    "Handeland et al. 2008",
                    "Heidt et al. 1988",
                    "Horak et al 2003",
                    "Manville 1978",
                    "Miller et al. 1998",
                    "Moks et al. 2006",
                    "Nevarez et al. 2005",
                    "Raoul et al. 2001",
                    "Snyder 1988b",
                    "Shimalov et al. 2000",
                    "Watson et al. 1981") ~ HostsSampled,
    
    Citation %in% c("Evans 2002",
                    "Kidder et al. 1989") ~ NumSamples,
    
    # Rejected
    Citation %in% c("Archer et al. 1986",
                    "Aukstikalniene et al. 2007",
                    "Clancey et al. 2010",
                    "Courtenay et al. 2006",
                    "Gompper et al. 2003",
                    "Gudmundsdottir and Skirnisson 2005",
                    "Jenkins et al 2006",
                    "Kloch et al. 2005",
                    "Lankester et al 2007",
                    "Martinek et al. 2001a",
                    "Martinello et al. 1997",
                    "Newman et al. 2001",
                    "Novobilsky et al 2007",
                    "Popiolek et al. 2007",
                    "Rodriguez and Carbonell 1998",
                    "Rosalino et al. 2006",
                    "Sakai et al. 1998",
                    "Stien et al 2002b",
                    "Szczesna and Popiolek 2007",
                    "Szczesna et al. 2008",
                    "Vicente et al. 2005") ~ 99999,
    
    # Addressed
    Citation == "Alassad et al 2008" &
      HostsSampled == 2096 ~ HostsSampled,
    Citation == "Alassad et al 2008" ~ 99999,
    
    Citation == "Anwar et al. 2000" &
      ParasiteCorrectedName == "Eimeria melis" ~ NumSamples,
    Citation == "Anwar et al. 2000" &
      NumSamples == 30 ~ 100,
    Citation == "Anwar et al. 2000" &
      NumSamples == 9 ~ 23,
    Citation == "Anwar et al. 2000" &
      NumSamples == 13 ~ 26,
    Citation == "Anwar et al. 2000" &
      NumSamples == 39 ~ 110,
    
    Citation == "Bjork et al. 2000" &
      Prevalence > 0.07 ~ HostsSampled,
    Citation == "Bjork et al. 2000" ~ 99999,  
    
    Citation == "Bourque et al. 2005" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Bourque et al. 2005" ~ 99999,
    
    Citation == "Compton et al. 2008" &
      HostsSampled == 332 & Prevalence == 0.087 ~ HostsSampled,
    Citation == "Compton et al. 2008" &
      HostsSampled == 128 & Prevalence == 0.078 ~ HostsSampled,
    Citation == "Compton et al. 2008" &
      HostsSampled == 278 & Prevalence == 0.058 ~ HostsSampled,
    Citation == "Compton et al. 2008" ~ 99999,
    
    Citation == "Chambers et al. 2002" &
      Prevalence > 0.5 ~ HostsSampled,
    Citation == "Chambers et al. 2002" ~ 99999,
    
    Citation == "Deem et al. 2005" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Deem et al. 2005" ~ 99999, 
    
    Citation == "Foster et al. 2004a" &
      (Prevalence == 0.836 | Prevalence == 0.098) ~ HostsSampled,
    Citation == "Foster et al. 2004a" ~ 99999, 
    
    Citation == "Hirvela-Koski et al. 2003" &
      !is.na(HostsSampled) ~ HostsSampled,
    Citation == "Hirvela-Koski et al. 2003" ~ 99999, 
    
    Citation == "Kimball et al. 2003" &
      HostsSampled != 6 ~ HostsSampled,
    Citation == "Kimball et al. 2003" ~ 99999, 
    
    Citation == "Lassnig et al. 1998" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Lassnig et al. 1998" ~ 99999, 
    
    Citation == "Lesmeister et al. 2008" &
      HostsSampled == 29 ~ HostsSampled,
    Citation == "Lesmeister et al. 2008" ~ 99999, 
    
    Citation == "Magi et al. 2009b" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Magi et al. 2009b" ~ 99999,
    
    Citation == "Martinez-Carrasco et al. 2007" &
      HostsSampled == 55 ~ HostsSampled,
    Citation == "Martinez-Carrasco et al. 2007" ~ 99999, 
    
    Citation == "Miterpakova et al. 2009" &
      HostsSampled > 1000 ~ HostsSampled,
    Citation == "Miterpakova et al. 2009" ~ 99999, 
    
    Citation == "Muoria et al 2005b" & 
      HostsSampled == 51 ~ HostsSampled,
    Citation == "Muoria et al 2005b" ~ 99999,
    
    Citation == "Nonaka et al. 1998" &
      Prevalence < 0.1 ~ NumSamples,
    Citation == "Nonaka et al. 1998" ~ 99999,
    
    Citation == "Page et al. 2005" &
      HostsSampled == 212 ~ HostsSampled,
    Citation == "Page et al. 2005" ~ 99999,
    
    Citation == "Page et al. 2009" &
      HostsSampled == 307 ~ HostsSampled,
    Citation == "Page et al. 2009" ~ 99999,
    
    Citation == "Pedersen et al. 2008" &
      !is.na(HostsSampled) ~ HostsSampled,
    Citation == "Pedersen et al. 2008" ~ 99999,
    
    Citation == "Snyder and Fitzgerald 1987" &
      HostSex == "All" & HostAge == "All" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Snyder and Fitzgerald 1987" ~ 99999,
    
    Citation == "Smith and Kok 2006" &
      HostsSampled == 5 ~ HostsSampled,
    Citation == "Smith and Kok 2006" ~ 99999,
    
    Citation == "Shimalov and Shimalov 2003b" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Shimalov and Shimalov 2003b" ~ 99999, 
    
    Citation == "Sexsmith et al. 2009" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Sexsmith et al. 2009" ~ 99999,
    
    Citation == "Valdmann et al. 2004" &
      !is.na(HostsSampled) ~ HostsSampled, 
    Citation == "Valdmann et al. 2004" ~ 99999,
    
    Citation == "Willingham et al. 1996" &
      HostSex == "All" & HostAge == "All" &
      Prevalence < 0.856 ~ HostsSampled,
    Citation == "Willingham et al. 1996" ~ 99999,
    
    Citation == "Wolfe et al. 2001a" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Wolfe et al. 2001a" ~ 99999,
    
    TRUE ~ NA
    
  ))

# Now doing citations with only one row because they are unlikely to have pseudosampling
Citations <- GMPD_Data_Filter %>% filter(is.na(SampleSize)) %>% 
  group_by(Citation) %>% summarise(n = n()) %>% 
  filter(n == 1) %>% pull(Citation)
# GMPD_Data %>% filter(Citation %in% Citations) %>% View()

GMPD_Data_Filter <- GMPD_Data_Filter %>% 
  mutate(SampleSize = case_when(
    !is.na(SampleSize) ~ SampleSize,
    
    # Kept
    Citation %in% c("Alexander et al. 1993a",
                    "Allsopp et al 1999",
                    "Appleyard et al. 1998",
                    "Bain et al 1993",
                    "Bakker et al. 2006",
                    "Barker et al 1995",
                    "Beard et al. 1999",
                    "Beechler et al. 2009",
                    "Beldomenico et al. 2005",
                    "Bengis et al 1996",
                    "Blitvich et al. 2009",
                    "Bollinger and Welch 1996",
                    "Boomker 1986b",
                    "Boomker and Vermaak 1986a",
                    "Boomker et al 1989a",
                    "Borchers et al. 2005",
                    "Borchers et al. 2008",
                    "Brillhart et al 1994a",
                    "Briscoe et al. 1993",
                    "Butler and Khan 1992",
                    "Bwangamoi et al. 1993",
                    "Carpi et al. 2008",
                    "Casulli et al. 2001",
                    "Ceballos et al. 2006",
                    "Chambers et al. 2010",
                    "Clifton-Hadley et al. 1995",
                    "Coleman et al. 1994",
                    "Costello et al. 2006",
                    "Crist et al 1999",
                    "Dalley et al. 2008",
                    "Damien et al. 2002",
                    "Drewe 2010",
                    "Dubay et al 1995",
                    "Dubey and Speer 1985",
                    "Dubey et al 1999b",
                    "Dubey et al. 2004c",
                    "Duh et al. 2005",
                    "Dyer and Huffman 1999",
                    "East et al. 2004",
                    "Ferte et al 2000",
                    "Fiorello et al. 2007",
                    "Fischer et al 1995",
                    "Fischthal and Martin 1977",
                    "Fishman et al. 2004",
                    "Flowers 1996",
                    "Folstad et al 1991",
                    "Foreyt et al. 1999",
                    "Franklin et al. 2008",
                    "Frolich 1995",
                    "Frolich et al 2005a",
                    "Gallivan et al 1998",
                    "Garcia-Bocanegra et al. 2010",
                    "Garcia-Sanchez et al 2007",
                    "Gaydos et al. 2007a",
                    "Goff et al 1993",
                    "Golezardy and Horak 2007b",
                    "Gomes et al. 2007",
                    "Gortazar et al. 1998b",
                    "Guerrero et al. 2010",
                    "Guzman et al. 2008",
                    "Hamir and Dubey 2001",
                    "Hamir et al. 1993b",
                    "Hamir et al. 1999b",
                    "Horak et al 1987a",
                    "Hove and Mukaratirwa 2005",
                    "Ingebrigtsen et al 1986",
                    "Isogai et al. 1994",
                    "Izdebska and Fryderyk 2000",
                    "James et al. 2002",
                    "Jankovska et al. 2010",
                    "Johnson 1975",
                    "Jolles et al 2008",
                    "Joly and Messier 2004",
                    "Joly and Messier 2004b",
                    "Jordan and Bird 1958",
                    "Kadulski et al 1996",
                    "Karbowiak et al. 2009",
                    "Keppner 1971",
                    "Kharchenko et al. 2008",
                    "Kimura et al 1995",
                    "Kirkpatrick et al. 1986",
                    "Klenavic et al. 2008",
                    "Knapp et al. 2008",
                    "Kocan et al 1986a",
                    "Kock et al 1995b",
                    "Kollars and Ladine 1999",
                    "Konig et al. 2006",
                    "Krametter et al 2004",
                    "Krecek et al 1990",
                    "Kruger et al 1988",
                    "Ladd-Wilson et al 2000",
                    "Lainson et al. 1990",
                    "Lecount 1981",
                    "Li et al 2005",
                    "Lindsay et al 1988",
                    "Liz et al 2002",
                    "Loken et al 1985",
                    "Lyashchenko et al 2008",
                    "Maas 1993",
                    "Macdonald et al. 1999",
                    "Majlathova et al. 2007",
                    "Maldonado and Kirkland 1986",
                    "Mannelli et al. 1993",
                    "Mares et al 1984",
                    "Mathews et al. 2006",
                    "McKown et al. 1991",
                    "McNeill and Rau 1987",
                    "McOrist et al. 1991",
                    "Miterpakova et al. 2006",
                    "Moks et al 2008",
                    "Morandi et al 2006",
                    "Morley and Hugh-Jones 1989",
                    "Murphy et al. 2010",
                    "Myburgh et al 1990",
                    "Nava et al. 2008",
                    "Nebbia et al 2000",
                    "Newell et al. 1997",
                    "Newman et al. 2002",
                    "Nielsen et al 2000",
                    "Nutter et al. 1998",
                    "O'Brien et al 2001",
                    "O'Brien et al 2008",
                    "O'Toole et al. 1994",
                    "Oates et al 1999",
                    "Ogunremi et al 2002",
                    "Oporto et al 2003b",
                    "Otranto et al. 2009",
                    "Peirce and Neal 1974",
                    "Pence et al. 2000",
                    "Petavy et al. 1991",
                    "Pfukenyi et al. 2009",
                    "Pietrokovsky et al. 1991",
                    "Pietrzak and Pung 1998",
                    "Polley 1986",
                    "Pozio et al. 1997",
                    "Pozio et al. 2004",
                    "Prestrud et al 1992",
                    "Prestrud et al. 2007",
                    "Qureshi et al 1989",
                    "Raizman et al. 2009",
                    "Renter et al 2001",
                    "Rhyan et al 2001",
                    "Ribas et al. 2009",
                    "Richardson and Barger 2005",
                    "Sacks et al. 2003",
                    "Sanmartin et al. 1992",
                    "Santin-duran et al 2000",
                    "Secord et al. 1980",
                    "Sherrard-Smith et al. 2009",
                    "Siepierski et al 1990",
                    "Simpson 2000",
                    "Slajchert et al. 1997",
                    "Snyder et al. 1989c",
                    "Southey et al. 2002",
                    "Speck et al 2008",
                    "Sreter et al. 2003a",
                    "Stagg et al 1987",
                    "Steen et al 2005",
                    "Tackmann et al. 2006",
                    "Takahashi et al 2001a",
                    "Tampieri et al 2008",
                    "Telford et al 1988",
                    "Torrence et al. 1992",
                    "Torres et al. 1997",
                    "Trap 1993",
                    "Uni et al. 2004",
                    "Ursprung et al. 2006",
                    "Van Campen et al 2001",
                    "Walls et al 1997",
                    "Wanha et al. 2005",
                    "Wanyangu et al 1987",
                    "Wild et al. 2006",
                    "Wisnivesky-Colli et al. 1992",
                    "Witmer et al. 2010",
                    "Yimam et al. 2001",
                    "Zimmerman et al. 2008") ~ HostsSampled,
    
    # Rejected
   Citation %in% c("Alexander and Appel 1994",
                   "Alexander et al. 1996",
                   "Bates 2003",
                   "Beck et al. 2008",
                   "Berrilli et al 2002",
                   "Bildfell et al 2007",
                   "Black et al. 1996",
                   "Blake et al. 2006",
                   "Boomker and Taylor 2004",
                   "Breagoli et al. 2006",
                   "Briones et al. 2000",
                   "Burtscher and Url 2007",
                   "Byrd et al. 1967",
                   "Cancrini et al. 2008",
                   "Colborne 1985",
                   "Cook et al 1997",
                   "De Meneghi et al 2002",
                   "de Mera et al. 2008",
                   "Delpietro et al. 1997",
                   "Delpietro et al. 2009",
                   "Despres et al 1995",
                   "Drewe et al. 2009b",
                   "Dubey and Lin 1994",
                   "Dubey et al. 1990",
                   "Dubey et al. 1991",
                   "Dubey et al. 1996",
                   "Ferreyra et al. 2009",
                   "Fitzgerald et al. 2008",
                   "Foronda et al. 2007",
                   "Gall et al 2000",
                   "Gerhold et al. 2005",
                   "Gerhold et al. 2007",
                   "Giacometti et al 2002",
                   "Glass et al. 1994",
                   "Haas et al. 1996",
                   "Hamir 2010",
                   "Hamir et al. 1995a",
                   "Hammer et al. 2004",
                   "Handeland 2002",
                   "Hernandez et al 1986",
                   "Hewicker et al. 1990",
                   "Hofmeyr et al. 2004",
                   "Honer et al. 2006",
                   "Ingram 1941",
                   "Inoshima et al 2002",
                   "Jessup et al 1990",
                   "Keet et al 1996a",
                   "Langley et al. 1994",
                   "Lehmkuhl et al 2001",
                   "Letshwenyo et al 2006",
                   "Literak et al. 2006",
                   "Lopez-Pena et al. 1994",
                   "Machida et al. 1993",
                   "Magnarelli et al. 1995",
                   "Marco et al 2007",
                   "Marco et al 2009d",
                   "Maritim et al 1992",
                   "Martella et al. 2002",
                   "Martin-Atance et al. 2005",
                   "McCollough and Pollard 1993",
                   "Mech et al. 1997",
                   "Megid et al. 2009",
                   "Millan et al. 2008",
                   "Miller and Harkema 1960",
                   "Miller et al. 2007",
                   "Moravkova et al 2008",
                   "Ndiaye et al. 2003",
                   "Neiffer et al. 2002",
                   "Nietfeld and Pollock 2002",
                   "Ninomiya and Ogata 2005",
                   "Noviana et al. 2004",
                   "Oleaga et al 2008",
                   "Oleaga et al 2008b",
                   "Paez et al. 2005",
                   "Pandey et al 1994",
                   "Pence et al. 1995",
                   "Penzhorn et al. 1998",
                   "Perez et al. 2001",
                   "Perz and Le Blancq 2001",
                   "Pietsch et al. 2002",
                   "Pitt and Jordan 1994",
                   "Prestrud et al. 2008",
                   "Priemer and Lux 1994",
                   "Ryser-Degiorgis et al. 2005",
                   "Saito and Little 1997",
                   "Sargeant et al 1999",
                   "Schmitt et al. 1987",
                   "Schurr et al 1988",
                   "Segovia et al. 2001a",
                   "Sillero-Zubiri et al. 1996",
                   "Silva et al. 2000",
                   "Simpson 1996",
                   "Simpson and Gavier-Widen 2000",
                   "Smith et al. 1995",
                   "Snyder 1989",
                   "Stefancikova 1994",
                   "Taylor et al. 2002",
                   "Tessaro and Forbes 1986",
                   "Tsukada et al. 2000",
                   "Wade et al. 1989",
                   "Wahl and Bain 1995",
                   "Weiler et al. 1995",
                   "Whitlaw and Lankester 1994",
                   "Williams et al. 2005",
                   "Woods et al 1997",
                   "Worley et al. 1990",
                   "Xiao et al. 2005") ~ 99999,
   
    # Addressed
    Citation == "Calderini et al. 2009" ~ 129, # check
    Citation == "Courtenay et al. 1994" ~ 25,
    
    TRUE ~ NA
  ))

Citations <- GMPD_Data_Filter %>% filter(is.na(SampleSize) & HostsSampled != NumSamples) %>% 
  pull(Citation)
# GMPD_Data %>% filter(Citation %in% Citations) %>% View()

GMPD_Data_Filter <- GMPD_Data_Filter %>% 
  mutate(SampleSize = case_when(
    !is.na(SampleSize) ~ SampleSize,
    
    # Kept
    Citation %in% c("Anderson et al 1987b",
                    "Anisimova 2002",
                    "Boomker et al 1995",
                    "Dubey et al. 1994",
                    "Harrison et al. 2003",
                    "Horak et al 2007",
                    "Kingscote et al Bohac 1986",
                    "Kocan et al 1986b",
                    "Smith et al. 2003",
                    "Van Den Bussche et al. 1987") ~ NumSamples,
    
    Citation %in% c("Anderson 1982",
                    "Fellis et al. 2003",
                    "Horak and MacIvor 1987",
                    "Johnson et al. 1994",
                    "Skotarczak et al 2008",
                    "Tall Timbers 1992") ~ HostsSampled,
    
    # Rejected
    Citation %in% c("Campbell et al 1989",
                    "Campbell et al 1994",       
                    "Cassirer et al 2001",       
                    "Chambers et al. 2009",      
                    "Degiorgis et al 2000b",     
                    "Ferroglio et al 1998",      
                    "Horak et al 1983a",         
                    "Lavin et al 1997",          
                    "Nilssen and Rolf 1994",     
                    "Okada et al 1984",          
                    "Penzhorn 1984",             
                    "Riley et al 1987",          
                    "Rosatte et al. 2006",       
                    "Schaffer et al. 1981",      
                    "Spraker et al 1984",        
                    "Takai et al 2004",          
                    "Taylor et al 1996b",        
                    "Telford and Forrester 1991a",
                    "Wobeser et al 1985") ~ 99999,
    
    # Addressed
    # Citation %in% c()
    Citation == "Addison et al 1988" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Addison et al 1988" ~ 99999,
    
    Citation == "Ball et al 2001" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Ball et al 2001" ~ 99999,
    
    Citation == "Banks and Ashley 2000" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Banks and Ashley 2000" ~ 99999, 
    
    Citation == "Barnard et al 1989" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Barnard et al 1989" ~ 99999, 
    
    Citation == "Beringer et al 2000" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Beringer et al 2000" ~ 99999,
    
    Citation == "Boomker et al 1983" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Boomker et al 1983" ~ 99999, 
    
    Citation == "Boomker et al 1991b" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Boomker et al 1991b" ~ 99999, 
    
    Citation == "Boomker et al 1991e" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Boomker et al 1991e" ~ 99999, 
    
    Citation == "Boromisa and Grimstad 1987" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Boromisa and Grimstad 1987" ~ 99999, 
    
    Citation == "Chaparro et al 1990" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Chaparro et al 1990" ~ 99999,
    
    Citation == "Claxton et al 1992" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Claxton et al 1992" ~ 99999,
    
    Citation == "Conner et al 2000" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Conner et al 2000" ~ 99999,
    
    Citation == "Cypher et al. 1998" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Cypher et al. 1998" ~ 99999,
    
    Citation == "de la Fuente et al 2001" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "de la Fuente et al 2001" ~ 99999,
    
    Citation == "De Villiers et al 1985" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "De Villiers et al 1985" ~ 99999,
    
    Citation == "De Vos et al 2001" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "De Vos et al 2001" ~ 99999, 
    
    Citation == "Dubey and Speer 1986" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Dubey and Speer 1986" ~ 99999,
    
    Citation == "Dunbar et al 1986" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Dunbar et al 1986" ~ 99999,
    
    Citation == "Ferrer et al 1998" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Ferrer et al 1998" ~ 99999,
    
    Citation == "Fletcher et al 1985" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Fletcher et al 1985" ~ 99999,
    
    Citation == "Forrester et al 1994" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Forrester et al 1994" ~ 99999,
    
    Citation == "Fournier et al 1986" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Fournier et al 1986" ~ 99999,
    
    Citation == "Gallivan et al 1989" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Gallivan et al 1989" ~ 99999,
    
    Citation == "Gill et al 1993" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Gill et al 1993" ~ 99999,
    
    Citation == "Gortazar et al. 1994" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Gortazar et al. 1994" ~ 99999,
    
    Citation == "Goyal et al 1992" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Goyal et al 1992" ~ 99999,
    
    Citation == "Halvorsen and Bye 1999" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Halvorsen and Bye 1999" ~ 99999,
    
    Citation == "Hein et al 1991" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Hein et al 1991" ~ 99999,
    
    Citation == "Hoberg et al 2002" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Hoberg et al 2002" ~ 99999,
    
    Citation == "Irvine et al 2000" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Irvine et al 2000" ~ 99999,
    
    Citation == "Isogai et al 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Isogai et al 1996" ~ 99999,
    
    Citation == "Jarvinen and Hedberg 1993" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Jarvinen and Hedberg 1993" ~ 99999,
    
    Citation == "Johnson et al 1986a" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Johnson et al 1986a" ~ 99999,
    
    Citation == "Johnson et al 1986b" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Johnson et al 1986b" ~ 99999,
    
    Citation == "Jorgensen and Vigh-Larsen 1986" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Jorgensen and Vigh-Larsen 1986" ~ 99999,
    
    Citation == "Keet et al 1997" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Keet et al 1997" ~ 99999,
    
    Citation == "Kitamura et al 1997" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Kitamura et al 1997" ~ 99999,
    
    Citation == "Kollars 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Kollars 1996" ~ 99999,
    
    Citation == "Kutz et al 2000" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Kutz et al 2000" ~ 99999,
    
    Citation == "Lankester and Luttich 1988" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Lankester and Luttich 1988" ~ 99999,
    
    Citation == "Leon-Vizcaino et al 1999" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Leon-Vizcaino et al 1999" ~ 99999,
    
    Citation == "Li et al 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Li et al 1996" ~ 99999,
    
    Citation == "Lindsay et al 1991" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Lindsay et al 1991" ~ 99999,
    
    Citation == "Lindsay et al 1999" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Lindsay et al 1999" ~ 99999,
    
    Citation == "Lux et al 1997" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Lux et al 1997" ~ 99999,
    
    Citation == "Mayer et al 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Mayer et al 1996" ~ 99999,
    
    Citation == "Mclean et al 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Mclean et al 1996" ~ 99999,
    
    Citation == "Murphy 1989" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Murphy 1989" ~ 99999,
    
    Citation == "New et al 1993" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "New et al 1993" ~ 99999,
    
    Citation == "O'Brien et al 2002" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "O'Brien et al 2002" ~ 99999,
    
    Citation == "Oates et al 2000" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Oates et al 2000" ~ 99999,
    
    Citation == "Pederson et al 1985" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Pederson et al 1985" ~ 99999,
    
    Citation == "Perez et al 1999" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Perez et al 1999" ~ 99999,
    
    Citation == "Perry et al 1985" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Perry et al 1985" ~ 99999,
    
    Citation == "Pletcher et al 1988" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Pletcher et al 1988" ~ 99999,
    
    Citation == "Qureshi et al 1994" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Qureshi et al 1994" ~ 99999,
    
    Citation == "Riley et al. 2007" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Riley et al. 2007" ~ 99999,
    
    Citation == "San Miguel et al 2001" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "San Miguel et al 2001" ~ 99999,
    
    Citation == "Santin-duran et al 2001" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Santin-duran et al 2001" ~ 99999,
    
    Citation == "Schmitt et al 1997" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Schmitt et al 1997" ~ 99999,
    
    Citation == "Shulaw et al 1986" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Shulaw et al 1986" ~ 99999,
    
    Citation == "Slomke et al 1995" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Slomke et al 1995" ~ 99999,
    
    Citation == "Stuve 1986" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Stuve 1986" ~ 99999,
    
    Citation == "Thomas 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Thomas 1996" ~ 99999,
    
    Citation == "Vanek et al 1996" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Vanek et al 1996" ~ 99999,
    
    Citation == "Vicente and Gortazar 2001" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Vicente and Gortazar 2001" ~ 99999,
    
    Citation == "Webster and Frandsen 1994" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Webster and Frandsen 1994" ~ 99999,
    
    Citation == "Westrom and Anderson 1992" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Westrom and Anderson 1992" ~ 99999,
    
    Citation == "Westrom et al 1985" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Westrom et al 1985" ~ 99999,
    
    Citation == "Williams et al 1993" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Williams et al 1993" ~ 99999,
    
    Citation == "Yokohata and Suzuki 1993" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Yokohata and Suzuki 1993" ~ 99999,
    
    Citation == "Zarnke and Erickson 1990" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Zarnke and Erickson 1990" ~ 99999,
    
    
    
    Citation == "Pybus 1990" &
      HostCorrectedName == "Alces alces" &
      is.na(HostsSampled) & is.na(NumSamples) ~ NumSamples,
    
    Citation == "Pybus 1990" &
      HostCorrectedName == "Cervus canadensis" &
      ParasiteCorrectedName == "Dictyocaulus viviparus" ~ 41,
    Citation == "Pybus 1990" &
      HostCorrectedName == "Cervus canadensis" &
      is.na(HostsSampled) & is.na(NumSamples) ~ NumSamples,
    
    Citation == "Pybus 1990" &
      HostCorrectedName == "Odocoileus hemionus" ~ HostsSampled,
    
    Citation == "Pybus 1990" &
      HostCorrectedName == "Odocoileus virginianus" &
      ParasiteCorrectedName == "Dictyocaulus viviparus" ~ 54,
    Citation == "Pybus 1990" &
      HostCorrectedName == "Odocoileus virginianus" &
      is.na(HostsSampled) & is.na(NumSamples)~ NumSamples,
    
    Citation == "Pybus 1990" ~ 99999,
    
    
    
    Citation == "Anderson et al 1990" &
      HostAge == "Adult" ~  NumSamples,
    Citation == "Anderson et al 1990" ~ 99999, 
    
    Citation == "Atkinson et al 1993" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Atkinson et al 1993" ~ 99999,
    
    Citation == "Bagrade et al. 2009" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Bagrade et al. 2009" ~ 99999, 
    
    Citation == "Bogaczyk et al 1993" &
      NumSamples == 679 ~ HostsSampled,
    Citation == "Bogaczyk et al 1993" ~ 99999, 
    
    Citation == "Boomker 1991b" &
      HostsSampled == NumSamples &
      NumSamples > 1 ~ HostsSampled,
    Citation == "Boomker 1991b" ~ 99999,
    
    Citation == "Boomker et al 1989b" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled, 
    Citation == "Boomker et al 1989b" ~ 99999, 
    
    Citation == "Boomker et al 1996" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled, 
    Citation == "Boomker et al 1996" ~ 99999,
    
    Citation == "Borecka et al. 2008" &
      HostsSampled == 214 ~ HostsSampled,
    Citation == "Borecka et al. 2008" ~ 99999, 
    
    Citation == "Bye 1987" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Bye 1987" ~ 99999,
    
    Citation == "Chomel et al. 1998" &
      HostCorrectedName == "Ursus americanus" &
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Chomel et al. 1998" &
      HostCorrectedName == "Ursus americanus" ~ 99999,
    Citation == "Chomel et al. 1998" &
      LocationName != "Alaska" ~ HostsSampled,
    Citation == "Chomel et al. 1998" ~ 99999, 
    
    Citation == "Courtenay et al. 2002" &
      SamplingType == "Serology" ~ HostsSampled,
    Citation == "Courtenay et al. 2002" ~ 99999, 
    
    Citation == "Delgiudice et al 1997" &
      HostsSampled > 2 ~ HostsSampled,
    Citation == "Delgiudice et al 1997" ~ 99999, 
    
    Citation == "Drozdz et al. 1998" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Drozdz et al. 1998" ~ 99999,
    
    Citation == "Dunbar et al 1990" &
      Prevalence > 0.35 ~ HostsSampled,
    Citation == "Dunbar et al 1990" ~ 99999, 
    
    Citation == "Fuller 1986" &
      HostsSampled == 37 ~ HostsSampled,
    Citation == "Fuller 1986" ~ 99999, 
    
    Citation == "Gabriel et al. 2009a" &
      Prevalence > 0.2 ~ HostsSampled,
    Citation == "Gabriel et al. 2009a" ~ 99999,  
    
    Citation == "Gonzalez-Candela et al. 2007" &
      LocationName != "Spain" ~ HostsSampled,
    Citation == "Gonzalez-Candela et al. 2007" ~ 99999, 
    
    Citation == "Grimstad et al 1986" &
      is.na(HostAge) & is.na(HostSex) &
      Prevalence > 0.5 ~ HostsSampled,
    Citation == "Grimstad et al 1986" ~ 99999,  
    
    Citation == "Hamir et al. 1999a" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Hamir et al. 1999a" ~ 99999,
    
    Citation == "Handeland and Gibbons 2001" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Handeland and Gibbons 2001" ~ 99999, 
    
    Citation == "Hill et al. 1992" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Hill et al. 1992" ~ 99999, 
    
    Citation == "Hill et al. 1998" &
      is.na(HostsSampled) ~ NumSamples,
    Citation == "Hill et al. 1998" ~ HostsSampled,
    
    Citation == "Holzman et al. 1992" &
      ParasiteCorrectedName == "Dirofilaria immitis" &
      NumSamples != 17 ~ 99999,
    Citation == "Holzman et al. 1992" &
      ParasiteCorrectedName == "Toxoplasma gondii" &
      HostSex != "All" ~ 99999,
    Citation == "Holzman et al. 1992" ~ NumSamples,
    
    Citation == "Horak et al 1983b" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Horak et al 1983b" ~ 99999, 
    
    Citation == "Horak et al 1988b" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Horak et al 1988b" ~ 99999, 
    
    Citation == "Jefferies et al. 1990" &
      HostsSampled == 56 ~ HostsSampled,
    Citation == "Jefferies et al. 1990" ~ 99999, 
    
    Citation == "Junge et al. 2007" &
      HostsSampled > 157 ~ HostsSampled,
    Citation == "Junge et al. 2007" ~ 99999, 
    
    Citation == "Keet et al 1996b" &
      HostsSampled == NumSamples &
      HostsSampled != 10 ~ HostsSampled,
    Citation == "Keet et al 1996b" ~ 99999, 
    
    Citation == "Kinjo 1987" &
      HostsSampled == NumSamples &
      Prevalence < 0.08 ~ HostsSampled, # Might revisit
    Citation == "Kinjo 1987" ~ 99999,
    
    Citation == "Kopecna et al 2006a" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Kopecna et al 2006a" ~ 99999,
    
    Citation == "L'Heureux et al 1996" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "L'Heureux et al 1996" ~ 99999, 
    
    Citation == "Lankester and Fong 1998" &
      HostCorrectedName == "Alces alces" ~ NumSamples,
    Citation == "Lankester and Fong 1998" ~ 99999, 
    
    Citation == "Martin-Atance et al. 2006" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Martin-Atance et al. 2006" &
      HostsSampled == 118 ~ HostsSampled,
    Citation == "Martin-Atance et al. 2006" ~ 99999, 
    
    Citation == "Mech et al. 2008" &
      HostsSampled == 518 ~ HostsSampled,
    Citation == "Mech et al. 2008" ~ 99999, 
    
    Citation == "Nettleton et al 1986" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled, 
    Citation == "Nettleton et al 1986" ~ 99999, 
    
    Citation == "Roffe et al 1999" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled, 
    Citation == "Roffe et al 1999" ~ 99999, 
    
    Citation == "Santin-Duran et al 2008" &
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Santin-Duran et al 2008" ~ 99999, 
    
    Citation == "Singer et al 2000" &
      SamplingType == "Serology" ~ HostsSampled,
    Citation == "Singer et al 2000" ~ 99999, 
    
    Citation == "Sugar 1997" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Sugar 1997" ~ 99999,
    
    Citation == "Thieking et al. 1992" &
      LocationName != "northern Minnesota and Wisconsin" ~ HostsSampled,
    Citation == "Thieking et al. 1992" ~ 99999, 
    
    Citation == "Valcarcel and Romero 2002" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Valcarcel and Romero 2002" ~ 99999, 
    
    Citation == "Webb et al 1987" &
      NumSamples == 170 ~ 180,
    Citation == "Webb et al 1987" &
      HostsSampled == NumSamples ~ HostsSampled,
    Citation == "Webb et al 1987" ~ 99999, 
    
    TRUE ~ NA
  ),
  Prevalence = case_when(
    Citation == "Kinjo 1987" &
      SampleSize != 99999 &
      ParasiteCorrectedName == "Leptospira interrogans" ~ 42/404, 
    Citation == "Kinjo 1987" &
      SampleSize != 99999 &
      ParasiteCorrectedName == "Toxoplasma gondii" ~ 41/765, 
    
    Citation == "Martin-Atance et al. 2006" &
      HostsSampled == 118 ~ 5/118,
    
    Citation == "Pybus 1990" &
      HostCorrectedName == "Cervus canadensis" &
      ParasiteCorrectedName == "Dictyocaulus viviparus" ~ 5/41,
    Citation == "Pybus 1990" &
      HostCorrectedName == "Odocoileus virginianus" &
      ParasiteCorrectedName == "Dictyocaulus viviparus" ~ 3/54,
    Citation == "Pybus 1990" &
      HostCorrectedName == "Odocoileus virginianus" &
      ParasiteCorrectedName == "Fascioloides magna" ~ 2/140,
    Citation == "Pybus 1990" &
      HostCorrectedName == "Odocoileus virginianus" &
      ParasiteCorrectedName == "Taenia hydatigena" ~ 7/147,
    
    TRUE ~ Prevalence
  ))

Citations <- GMPD_Data_Filter %>% filter(is.na(SampleSize) & HostsSampled < 5) %>% 
  pull(Citation)
# GMPD_Data %>% filter(Citation %in% Citations) %>% View()

GMPD_Data_Filter <- GMPD_Data_Filter %>% 
  mutate(SampleSize = case_when(
    !is.na(SampleSize) ~ SampleSize,

    # Kept
    Citation %in% c("Alexander et al. 1994",     
                    "Aranaz et al. 2004",        
                    "Artois and Remond 1994",    
                    "Bechara et al 2000",        
                    "Beck et al. 2009",          
                    "Boomker et al 1986",        
                    "Boomker et al 1989c",       
                    "Brillhart et al 1994b",     
                    "Brown et al. 1993",         
                    "Chomel et al. 2004",        
                    "Conder and Loveless 1978",  
                    "Conti 1984",                
                    "Crosbie et al 1997",        
                    "Dalimi et al. 2002",        
                    "Dubey et al. 2004a",        
                    "Filoni et al. 2006",        
                    "Gallivan and Surgeoner 1995",
                    "Golezardy and Horak 2006a", 
                    "Hamblin et al 1990",        
                    "Hirvela-Koski et al. 1985", 
                    "Hoberg et al 2002b",        
                    "Horak et al 1987b",         
                    "Horak et al. 1999",         
                    "Hubert et al. 1980",        
                    "Hurnikova et al. 2007",     
                    "Hwang et al. 2002",         
                    "Jeffery et al. 2004",       
                    "Kollars et al 2000",        
                    "Kollars et al. 1999",       
                    "Krivanec et al. 1976",      
                    "Laurenson et al. 1998",     
                    "Leiby et al. 1970",         
                    "Luaces et al. 2008",        
                    "Madic et al. 1993",         
                    "Malczewski et al. 2008",    
                    "Malek et al. 1961",         
                    "Marchiondo et al. 1986",    
                    "Mech and Goyal 1995",       
                    "Meshgi et al. 2009",        
                    "Millan et al. 2007b",       
                    "Millan et al. 2009c",       
                    "Miller et al. 2000a",       
                    "Miller et al. 2008b",       
                    "Miller et al. 2009a",       
                    "Morgan et al 2005",         
                    "Ohashi et al. 2001b",       
                    "Oliver et al 1999",         
                    "Osofsky et al. 1996",       
                    "Papadopoulos et al. 1997",  
                    "Patrick and Harrison 1995", 
                    "Penzhorn et al. 2002",      
                    "Pozio et al. 2001",         
                    "Quinn et al. 1976",         
                    "Rajsky et al 2002",         
                    "Rickard and Foreyt 1992",   
                    "Riemann et al. 1975",       
                    "Riemann et al. 1978",       
                    "Robbins and Deem 2002",     
                    "Rogers 1975",               
                    "Sacks 1998",                
                    "Sato et al. 1999b",         
                    "Senger et al. 1955",        
                    "Simpson et al. 2005",       
                    "Smith and Frenkel 1995",    
                    "Smith et al. 1992a",        
                    "Sobrino et al. 2007",       
                    "Stuen et al 2002",          
                    "Taylor et al 2005",         
                    "Thorne et al 1988",         
                    "Torres et al. 2001",        
                    "Tryland et al. 2005c",      
                    "Waldrup et al 1989a",       
                    "Wanyangu et al 1989",       
                    "Yabsley and Noblet 2002",   
                    "Yabsley et al 2002",        
                    "Zarnke and Ballard 1987",   
                    "Zieger et al 1998b") ~ HostsSampled,        
    
    # Rejected
    Citation %in% c("Barros et al. 1990",        
                    "Bartsch and Ward 1976",     
                    "Baszler et al 2000",        
                    "Boyce et al 1999",          
                    "Burcham et al. 2010",
                    "Camicas et al. 1972",
                    "Carnieli et al. 2009",      
                    "Carpenter et al. 1998",     
                    "Chae et al 1999",           
                    "Chen et al. 2008",          
                    "Cleaveland et al. 2005",    
                    "Daoust et al. 1996",        
                    "De Bosschere et al. 2007",  
                    "Degiorgis et al 2000a",     
                    "Degiorgis et al. 2001",     
                    "Dubay et al 2000",          
                    "Edmunds et al. 2008",       
                    "Favoretto et al. 2006",     
                    "Ferroglio et al 2000",      
                    "Fitzgerald et al 2000",     
                    "Foreyt et al 1994",         
                    "Foreyt et al. 1996",        
                    "Gasser et al 1999b",        
                    "Hamir et al. 1998",         
                    "Hirama et al. 2004",        
                    "Holman et al 1988",         
                    "Honour and Hickling 1993",  
                    "Huchzermeyer et al 2001",   
                    "Iori and Lanfranchi 1996",  
                    "Jenkins et al 2001",        
                    "Jessup et al 1993a",        
                    "Kamiya et al. 2003",        
                    "Kat et al. 1996",           
                    "Keet et al 1994",           
                    "Keet et al 2001",           
                    "Kelly and Sleeman 2003",    
                    "Kimber et al 2002",         
                    "Kingston et al 1985",       
                    "Kulonen and Boldina 1993",  
                    "Last et al. 1994",          
                    "Lemberger et al. 2005",     
                    "Little and Howerth 1999",   
                    "Lopez et al. 2009",         
                    "Lucientes and Castillo 1990",
                    "Machida et al. 1992",       
                    "MacIvor 1985",              
                    "MacIvor et al 1987",        
                    "Macko and Birova 1989",     
                    "Majoros and Sztojkov 1994", 
                    "Manna et al 1991",          
                    "Marco et al 2000",          
                    "Marucci et al. 2009",       
                    "Ngeranwa et al 1998",       
                    "Nilssen and Gjershaug 1988",
                    "Noon et al 2002",           
                    "Okaeme 1987",               
                    "Otranto et al. 2007a",      
                    "Priemer et al. 2002",       
                    "Rhyan et al 1994",          
                    "Ribas et al. 2004",         
                    "Rotstein et al. 1999b",     
                    "Ryser-Degiorgis et al 2009",
                    "Schultheiss et al 2007",
                    "Snyder et al. 1989a",       
                    "Sorensen et al. 2005",      
                    "Stuen et al 2006",          
                    "Sundar et al. 2008a",       
                    "Takahashi et al 2001b",     
                    "Theberge et al. 1994",      
                    "Thorne et al 1987",         
                    "Tryland et al. 2005",       
                    "Uzal et al. 2007",          
                    "Walker et al 1993",         
                    "Warsame and Steen 1989",    
                    "Whitby et al. 1997",        
                    "Winder et al. 2006") ~ 99999,       
    
    # Addressed
    Citation == "Aramini et al. 1999" & 
      LocationName %in% c("Courtenay, Vancouver Island",
                          "Ladysmith, Vancouver Island",
                          "Nanaimo, Vancouver Island",
                          "Port Alberni, Vancouver Island",
                          "Parksville, Vancouver Island, Canada",
                          "Sayward, Vancourver Island",
                          "Sooke, Vancourver Island",
                          "Vancouver Island, British Columbia") ~ HostsSampled,
    Citation == "Aramini et al. 1999" & 
      (HostSex == "All" | HostAge == "All") ~ HostsSampled,
    Citation == "Aramini et al. 1999" ~ 99999, 
    
    Citation == "Baeten et al 2007" &
      HostsSampled > 2 ~ HostsSampled,
    Citation == "Baeten et al 2007" ~ 99999, 
    
    Citation == "Ballard and Krausman 1997" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Ballard and Krausman 1997" ~ 99999, 
    
    Citation == "Borman et al. 2009" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Borman et al. 2009" ~ 99999, 
    
    Citation == "Bukva 1987" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Bukva 1987" ~ 99999, 
    
    Citation == "Carreno et al 2001" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Carreno et al 2001" ~ 99999, 
    
    Citation == "Clark et al 1985" &
      ParasiteCorrectedName == "Orbivirus Bluetongue virus" &
      LocationName =="Jacumba Mountains, California" &
      is.na(SamplingType) ~ 99999, 
    Citation == "Clark et al 1985" &
      ParasiteCorrectedName == "Respirovirus Human parainfluenza virus 3" &
      is.na(SamplingType) ~ 99999,
    Citation == "Clark et al 1985" ~ HostsSampled, 
    
    Citation == "Creel et al. 1997" &
      HostsSampled == 22 ~ HostsSampled,
    Citation == "Creel et al. 1997" ~ 99999,
    
    Citation == "Deem et al. 2002" &
      HostsSampled == 94 ~ HostsSampled, 
    Citation == "Deem et al. 2002" ~ 99999, 
    
    Citation == "Dietrich et al. 2005" &
      HostsSampled == 5 ~ HostsSampled,
    Citation == "Dietrich et al. 2005" ~ 99999,
    
    Citation == "Dipineto et al. 2007" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Dipineto et al. 2007" ~ 99999, 
    
    Citation == "Dubinsky et al. 2006" &
      (LocationName != "Poland" |
         LocationName != "Slovakia") ~ HostsSampled,
    Citation == "Dubinsky et al. 2006" ~ 99999,
    
    Citation == "Dunbar et al 1985" & 
      HostsSampled != 66 ~ HostsSampled,
    Citation == "Dunbar et al 1985" ~ 99999,
    
    Citation == "Durden and Horak 2004" & 
      HostsSampled != 322 ~ HostsSampled,
    Citation == "Durden and Horak 2004" ~ 99999,
    
    Citation == "El-Shehabi et al. 1999" &
      HostAge == "All" ~ HostsSampled,
    Citation == "El-Shehabi et al. 1999" ~ 99999,
    
    Citation == "Elnaiem et al. 2001" &
      HostsSampled == 14 ~ HostsSampled,
    Citation == "Elnaiem et al. 2001" ~ 99999,
    
    Citation == "Farid et al. 2010" &
      LocationName == "Yarmouth county, Nova Scotia, CANADA" ~ HostsSampled,
    Citation == "Farid et al. 2010" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Farid et al. 2010" ~ 99999,
    
    Citation == "Foster et al. 1998" &
      HostAge == "All" & is.na(HostSex) ~ HostsSampled,
    Citation == "Foster et al. 1998" ~ 99999,
    
    Citation == "Foster et al. 2003" &
      HostsSampled == 26 ~ HostsSampled,
    Citation == "Foster et al. 2003" ~ 99999,
    
    Citation == "Franklin et al. 2007a" &
      SamplingType == "Serology" ~ HostsSampled,
    Citation == "Franklin et al. 2007a" ~ 99999,
    
    Citation == "Franklin et al. 2007b" &
      HostCorrectedName == "Leopardus pardalis" &
      HostsSampled == 12 ~ HostsSampled,
    Citation == "Franklin et al. 2007b" &
      HostCorrectedName == "Lynx rufus" &
      HostsSampled == 43 ~ HostsSampled,
    Citation == "Franklin et al. 2007b" &
      HostCorrectedName == "Puma concolor" &
      HostsSampled == 31 ~ HostsSampled,
    Citation == "Franklin et al. 2007b" ~ 99999, 
    
    Citation == "Gascoyne et al. 1993a" &
      HostsSampled == 12 ~ HostsSampled,
    Citation == "Gascoyne et al. 1993a" ~ 99999,
    
    Citation == "Gavier-Widen et al. 2001" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Gavier-Widen et al. 2001" ~ 99999,
    
    Citation == "Grinder and Krausman 2001" &
      HostsSampled == 22 ~ HostsSampled,
    Citation == "Grinder and Krausman 2001" ~ 99999,
    
    Citation == "Guedegbe et al 1992" &
      is.na(NumSamples) ~ HostsSampled,
    Citation == "Guedegbe et al 1992" ~ 99999,
    
    Citation == "Harrison et al. 2004" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Harrison et al. 2004" ~ 99999,
    
    Citation == "Henn et al. 2009" &
      LocationName == "Israel" ~ HostsSampled,
    Citation == "Henn et al. 2009" ~ 99999,
    
    Citation == "Herrera et al 2008a" &
      HostsSampled == 123 ~ HostsSampled,
    Citation == "Herrera et al 2008a" ~ 99999,
    
    Citation == "Hofmeyer et al. 2000" &
      HostsSampled == 8 ~ HostsSampled,
    Citation == "Hofmeyer et al. 2000" ~ 99999,
    
    Citation == "Horak et al 1986a" &
      HostsSampled != 1 ~ HostsSampled,
    Citation == "Horak et al 1986a" ~ 99999,
    
    Citation == "Horak et al 2006" &
      LocationName != "South Africa" ~ HostsSampled,
    Citation == "Horak et al 2006" ~ 99999,
    
    Citation == "Hurtado et al 2004" &
      HostsSampled == 17 ~ HostsSampled,
    Citation == "Hurtado et al 2004" ~ 99999,
    
    Citation == "Jakubek et al. 2007" &
      LocationName != "Hungary" ~ HostsSampled,
    Citation == "Jakubek et al. 2007" ~ 99999,
    
    Citation == "Jenkins et al 2007" &
      HostsSampled == 56 ~ HostsSampled,
    Citation == "Jenkins et al 2007" ~ 99999,
    
    Citation == "Kikuchi et al. 2004" &
      HostsSampled != 96 ~ HostsSampled,
    Citation == "Kikuchi et al. 2004" ~ 99999,
    
    Citation == "Koizumi et al. 2009" &
      HostsSampled != 124 &
      SamplingType == "Serology" ~ HostsSampled,
    Citation == "Koizumi et al. 2009" ~ 99999,
    
    Citation == "Konig et al. 2005" &
      !str_detect(LocationName, "Southern Bavaria") ~ HostsSampled,
    Citation == "Konig et al. 2005" ~ 99999,
    
    Citation == "Kresta et al. 2009" &
      ParasiteCorrectedName %in% 
      c("Baylisascaris procyonis",
        "Gnathostoma procyonis",
        "Gyrosoma singularis",
        "Heterobilharzia americana") ~ HostsSampled,
    Citation == "Kresta et al. 2009" &
      LocationName != "Texas" ~ HostsSampled,
    Citation == "Kresta et al. 2009" ~ 99999,
    
    Citation == "Lahmar et al. 2009" &
      SamplingType == "DirectOther" ~ HostsSampled,
    Citation == "Lahmar et al. 2009" ~ 99999,
    
    Citation == "Leutenegger et al. 1999" &
      ParasiteCorrectedName %in% c("Alphacoronavirus Alphacoronavirus 1", 
                                   "Varicellovirus Felid alphaherpesvirus 1",
                                   "Vesivirus Feline calicivirus") &
      SamplingType == "Serology" ~ HostsSampled,
    Citation == "Leutenegger et al. 1999" &
      ParasiteCorrectedName %in% c("Gammaretrovirus Feline leukemia virus",
                                   "Lentivirus Feline immunodeficiency virus") &
      SamplingType == "PCR" ~ HostsSampled,
    Citation == "Leutenegger et al. 1999" ~ 99999,
    
    Citation == "Magi et al. 2009a" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Magi et al. 2009a" ~ 99999,
    
    Citation == "McIntosh et al 2007b" &
      HostsSampled == 29 ~ HostsSampled,
    Citation == "McIntosh et al 2007b" ~ 99999,
    
    Citation == "Mitchell et al. 2002" &
      (is.na(HostSex) | HostSex == "All") ~ HostsSampled,
    Citation == "Mitchell et al. 2002" ~ 99999,
    
    Citation == "Mock et al. 1991" &
      HostsSampled == 26 ~ HostsSampled, 
    Citation == "Mock et al. 1991" ~ 99999,
    
    Citation == "Morner et al. 2005" &
      HostsSampled == 7 ~ HostsSampled,
    Citation == "Morner et al. 2005" ~ 99999,
    
    Citation == "Munson et al. 2004" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Munson et al. 2004" ~ 99999,
    
    Citation == "Olson and Lindzey 2002" &
      HostsSampled == 16 ~ HostsSampled,
    Citation == "Olson and Lindzey 2002" ~ 99999,
    
    Citation == "Ostrowski et al. 2003" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Ostrowski et al. 2003" ~ 99999,
    
    Citation == "Pena et al. 2006" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Pena et al. 2006" ~ 99999,
    
    Citation == "Quadros et al. 2009" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Quadros et al. 2009" ~ 99999,
    
    Citation == "Riemann et al. 1975" &
      LocationName == "California (14 unspecified geographic areas)" ~ HostsSampled,
    Citation == "Riemann et al. 1975" ~ 99999,
    
    Citation == "Roelke et al. 2009" &
      HostsSampled != 64 ~ HostsSampled,
    Citation == "Roelke et al. 2009" ~ 99999,
    
    Citation == "Rossiter et al 2006" &
      HostsSampled != 344 ~ HostsSampled,
    Citation == "Rossiter et al 2006" ~ 99999,
    
    Citation == "Segovia et al. 2001b" &
      ParasiteCorrectedName %in% c("Ancylostoma caninum",
                                   "Angiostrongylus vasorum",
                                   "Capillaria plica",
                                   "Dipylidium caninum",
                                   "Mesocestoides litteratus",
                                   "Unicinaria stenocephala") &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Segovia et al. 2001b" &
      ParasiteCorrectedName %in% c("Dirofilaria immitis",
                                   "Taenia hydatigena",
                                   "Taenia multiceps",
                                   "Toxascaris leonina",
                                   "Toxocara canis",
                                   "Trichuris vulpis") &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Segovia et al. 2001b" ~ 99999,
    
    Citation == "Sobrino et al. 2006" &
      LocationName == "Northern Spain" ~ HostsSampled,
    Citation == "Sobrino et al. 2006" ~ 99999,
    
    Citation == "Sobrino et al. 2008a" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Sobrino et al. 2008a" ~ 99999,
    
    Citation == "Souza et al. 2009" &
      LocationName != "Eastern Tennessee" ~ HostsSampled,
    Citation == "Souza et al. 2009" ~ 99999,
    
    Citation == "Telford and Forrester 1991b" &
      Prevalence == (0.38|0.16) ~ 99999,
    Citation == "Telford and Forrester 1991b" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Telford and Forrester 1991b" ~ 99999,
    
    Citation == "Thomas and Hughes 1992" &
      LocationName != "Los Angeles County, California" ~ HostsSampled,
    Citation == "Thomas and Hughes 1992" ~ 99999,
    
    Citation == "Thompson et al. 2010" &
      HostsSampled == 58 ~ HostsSampled,
    Citation == "Thompson et al. 2010" ~ 99999,
    
    Citation == "Tocidlowski et al. 1997" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Tocidlowski et al. 1997" ~ 99999,
    
    Citation == "Upshall et al 1987" &
      HostsSampled != 40 ~ HostsSampled,
    Citation == "Upshall et al 1987" ~ 99999,
    
    Citation == "VerCauteren et al. 2008" &
      HostSex == "All" ~ HostsSampled,
    Citation == "VerCauteren et al. 2008" ~ 99999,
    
    Citation == "Ward et al 1997" &
      ParasiteCorrectedName %in% c("Bibersteinia trehalosi",
                                   "Pasteurella multocida") ~ HostsSampled,
    Citation == "Ward et al 1997" &
      ParasiteCorrectedName == "Mannheimia haemolytica" &
      Prevalence > 0.5 ~ HostsSampled,
    Citation == "Ward et al 1997" ~ 99999,
    
    Citation == "Woodroffe et al. 2005" & 
      HostSex == "All" ~ HostsSampled, 
    Citation == "Woodroffe et al. 2005" ~ 99999,
    
    Citation == "Zarnke et al. 2001b" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Zarnke et al. 2001b" ~ 99999,
    
    
    
    Citation == "Baeten et al 2007" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Baeten et al 2007" ~ 99999, 
    
    Citation == "Boomker et al 1984" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Boomker et al 1984" ~ 99999, 
    
    Citation == "Boomker et al 1991a" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Boomker et al 1991a" ~ 99999, 
    
    Citation == "Boomker et al 1991c" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Boomker et al 1991c" ~ 99999, 
    
    Citation == "Boomker et al 1991d" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Boomker et al 1991d" ~ 99999, 
    
    Citation == "Boomker et al 2000" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Boomker et al 2000" ~ 99999, 
    

    TRUE ~ NA
    ))

Citations <- GMPD_Data_Filter %>% filter(is.na(SampleSize) & !is.na(HostAge)) %>% 
  pull(Citation)
# GMPD_Data %>% filter(Citation %in% Citations) %>% View()

GMPD_Data_Filter <- GMPD_Data_Filter %>% 
  mutate(SampleSize = case_when(
    !is.na(SampleSize) ~ SampleSize,
    # Kept
    Citation %in% c("Adamska 2008",            
                    "Akerstedt et al. 2010",   
                    "Balestrieri et al. 2006", 
                    "Beldomenico et al. 2005a",
                    "Binninger et al. 1980",   
                    "Bridger et al. 2009",     
                    "Butler and Grundmann 1953",
                    "Ching et al. 2000",   
                    "Clifford et al. 2006",   
                    "Davidson et al. 1992",
                    "Dunbar and Giordano 2003",
                    "Durden and Richardson 2003",
                    "Foster et al. 2004b",
                    "Foster et al. 2006",
                    "Foster et al. 2007",
                    "Frolich et al 2006",
                    "Garcelon et al. 1992",
                    "Gicik et al. 2009",
                    "Goble and Cook 1942",
                    "Hersteinsson et al. 1993",
                    "Hwang et al. 2007",
                    "Jordan and Hayes 1959",
                    "Kita et al 2003",
                    "Kocan et al. 1999",
                    "Koppel et al 2007b",
                    "Labelle et al. 2001",
                    "Lamm et al. 1997",
                    "Lindsay et al 2005",
                    "Lindsay et al. 1996",
                    "Lindsay et al. 2001a",
                    "Machackova et al 2006",
                    "Magnarelli et al. 1991",
                    "Manangan et al 2007",
                    "Manfredi et al 2007b",
                    "Marsilio et al. 1997",
                    "McGee et al. 2006",
                    "McKinny et al 2006",
                    "McLaughlin et al. 1993",
                    "McNeil and Krogsdale 1952",
                    "Mech and Goyal 1993",
                    "Mech and Tracy 2001",
                    "Millan and Ferroglio 2001",
                    "Millan et al. 2004b",
                    "Miller et al. 2006b",
                    "Moll et al. 1995",
                    "Namekata et al. 2009",
                    "O'Toole et al. 1993",
                    "Pappas and Lunzman 1985",
                    "Paul-Murphy et al. 1994",
                    "Pence et al. 2003",
                    "Perez et al 2003",
                    "Pollock et al 2009",
                    "Pryor 1956",
                    "Radomski and Pence 1993",
                    "Reperant et al. 2007",
                    "Richardson and Gauthier 2003",
                    "Roelke et al. 1993",
                    "Roemer et al. 2000",
                    "Sacks and Caswell-Chen 2003",
                    "Simmons et al. 1980",
                    "Skirnisson et al. 1993",
                    "Smith 1943",
                    "Snyder et al. 1990",
                    "Spencer and Morkel 1993",
                    "Spencer et al. 1999",
                    "Stone and Pence 1978",
                    "Stuen et al 2002b",
                    "Stuht and Youatt 1972",
                    "Thalwitzer et al. 2010",
                    "Thornton et al. 1974",
                    "Vikoren et al. 2006",
                    "Waid and Pence 1988") ~ HostsSampled,
    
    Citation %in% c("Biek et al. 2006") ~ NumSamples,
    
    # Rejected
    Citation %in% c("Besser et al 2008",
                    "Fourie & Vrahimis 1989",
                    "Greenwood et al. 1997",
                    "Jaworski et al 1993",
                    "Marco et al 2009",
                    "Nyberg et al. 1992",
                    "Rehbinder et al 2004",
                    "Rickard et al 1993",
                    "Stuen et al 2002b") ~ 99999,
    
    
    # Addressed
    Citation == "Aguirre et al 1995" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Aguirre et al 1995" ~ 99999,
    
    Citation == "Anderson & Rowe 1998" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Anderson & Rowe 1998" ~ 99999,
    
    
    
    Citation == "Almberg et al. 2009" &
      (HostAge == "All" | is.na(HostAge)) ~ NumSamples,
    Citation == "Almberg et al. 2009" ~ 99999,
    
    Citation == "Almeria et al. 2002" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Almeria et al. 2002" ~ 99999,
    
    Citation == "Almeria et al. 2007" &
      HostCorrectedName == "Capreolus capreolus" ~ HostsSampled,
    Citation == "Almeria et al. 2007" &
      ParasiteCorrectedName == "Toxoplasma gondii" ~ HostsSampled,
    Citation == "Almeria et al. 2007" &
      LocationName != "Spain" ~ HostsSampled,
    Citation == "Almeria et al. 2007" ~ 99999,
    
    
    Citation == "Arens et al 2003" &
      SamplingType == "PCR" ~ HostsSampled, # Not choosing the highest prev this time because:
    Citation == "Arens et al 2003" ~ 99999, # "presumably caused by cross-reacting antibodies"
    
    Citation == "Arjo et al. 2003" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Arjo et al. 2003" ~ 99999,
    
    Citation == "Brochier et al. 2007" & 
      HostAge == "All" ~ HostsSampled,
    Citation == "Brochier et al. 2007" ~ 99999,
    
    Citation == "Fischer et al. 2005" & 
      HostAge == "All" ~ HostsSampled,
    Citation == "Fischer et al. 2005" ~ 99999,
    
    
    
    
    Citation == "Bagrade et al. 2008" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Bagrade et al. 2008" ~ 99999, 
    
    Citation == "Ballard et al. 2001" &
      Prevalence > 0.02 ~ HostsSampled,
    Citation == "Ballard et al. 2001" ~ 99999,
    
    Citation == "Barnard 1993" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Barnard 1993" ~ 99999,
    
    
    
    
    Citation == "Beldomenico et al. 2005b" &
      Prevalence > 0.3 ~ HostsSampled,
    Citation == "Beldomenico et al. 2005b" ~ 99999,
    
    Citation == "Biek et al. 2003" &
      HostsSampled == 52 ~ HostsSampled,
    Citation == "Biek et al. 2003" ~ 99999,
    
    Citation == "Bornstein et al. 2006" & 
      HostsSampled == 88 ~ HostsSampled,
    Citation == "Bornstein et al. 2006" ~ 99999,
    
    Citation == "Brossard et al. 2007" &
      is.na(HostSex) & is.na(HostAge) ~ HostsSampled,
    Citation == "Brossard et al. 2007" ~ 99999,
    
    
    
    Citation == "Cattet et al. 2004" & 
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Cattet et al. 2004" ~ 99999,
    
    Citation == "Chang et al. 1999" & 
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Chang et al. 1999" ~ 99999,
    
    Citation == "Chang et al. 2000" & 
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Chang et al. 2000" ~ 99999,
    
    Citation == "Ferroglio et al. 2009" & 
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Ferroglio et al. 2009" ~ 99999,
    
    Citation == "Foreyt and Lagerquist 1993" & 
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Foreyt and Lagerquist 1993" ~ 99999,
    
    Citation == "Foreyt et al. 2009" & 
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Foreyt et al. 2009" ~ 99999,
    
    
  
    Citation == "Chomel et al 1994" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Chomel et al 1994" ~ 99999,
    
    Citation == "Chomel et al. 1995" &
      HostCorrectedName == "Ursus americanus" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Chomel et al. 1995" &
      HostCorrectedName == "Ursus arctos" &
      HostSex == "All" & HostAge %in% c("All", "Adult") ~ HostsSampled,
    Citation == "Chomel et al. 1995" ~ 99999,
    
    Citation == "Citterio et al. 2003" &
      (HostAge == "All" | is.na(HostAge)) ~ HostsSampled,
    Citation == "Citterio et al. 2003" ~ 99999,
    
    Citation == "Clifton-Hadley et al. 1993" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Clifton-Hadley et al. 1993" ~ 99999, 
    
    Citation == "Crawshaw et al. 2008" &
      Prevalence != 0.559 & Prevalence > 0.1 ~ HostsSampled,
    Citation == "Crawshaw et al. 2008" ~ 99999,
    
    Citation == "Creel et al. 1995" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Creel et al. 1995" ~ 99999,
    
    Citation == "Davidson et al 1985" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Davidson et al 1985" ~ 99999,
    
    Citation == "Davidson et al. 2006a" &
      LocationName != "Norway" ~ HostsSampled,
    Citation == "Davidson et al. 2006a" ~ 99999,
    
    Citation == "Delahay et al 2007" &
      HostCorrectedName != "Vulpes vulpes" ~ HostsSampled,
    Citation == "Delahay et al 2007" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Delahay et al 2007" ~ 99999,
    
    Citation == "Di Cerbo et al. 2008b" &
      (is.na(HostAge) | HostAge == "All") ~ HostsSampled,
    Citation == "Di Cerbo et al. 2008b" ~ 99999,
    
    Citation == "Driciru et al. 2006" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Driciru et al. 2006" ~ 99999,
    
    Citation == "du Preez and Moeng 2004" &
      HostAge == "All" ~ HostsSampled,
    Citation == "du Preez and Moeng 2004" ~ 99999,
    
    Citation == "Dubey and Kistner 1985" &
      is.na(HostAge) ~ HostsSampled,
    Citation == "Dubey and Kistner 1985" ~ 99999,
    
    Citation == "Dubey et al. 1995" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Dubey et al. 1995" ~ 99999,
    
    Citation == "Dunbar et al 1999" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Dunbar et al 1999" ~ 99999,
    
    Citation == "Foley et al. 1999" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Foley et al. 1999" ~ 99999,
    
    
    Citation == "Fromont et al. 2000" & 
      ParasiteCorrectedName == "Gammaretrovirus Feline leukemia virus" &
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Fromont et al. 2000" & 
      ParasiteCorrectedName == "Lentivirus Feline immunodeficiency virus" &
      Prevalence < 0.2 &
      HostAge == "All" & HostSex == "All" ~ HostsSampled,
    Citation == "Fromont et al. 2000" ~ 99999,
    
    Citation == "Gallivan et al 1996" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Gallivan et al 1996" ~ 99999,
    
    Citation == "Ganley-Leal et al. 2007" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Ganley-Leal et al. 2007" ~ 99999,
    
    Citation == "Gehrt et al. 2010" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Gehrt et al. 2010" &
      HostSex == "All" & HostAge == "All" ~ NumSamples,
    Citation == "Gehrt et al. 2010" ~ 99999,
    
    Citation == "Gese et al. 1991" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Gese et al. 1991" ~ 99999,
    
    Citation == "Gese et al. 1997" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Gese et al. 1997" ~ 99999,
    
    Citation == "Gese et al. 2004" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Gese et al. 2004" ~ 99999,
    
    Citation == "Gonzalez-Acuna et al. 2007" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Gonzalez-Acuna et al. 2007" ~ 99999,
    
    Citation == "Guislain et al. 2008" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Guislain et al. 2008" ~ 99999,
    
    Citation == "Hagiwara et al. 2009" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Hagiwara et al. 2009" ~ 99999,
    
    Citation == "Hanosset et al. 2008" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Hanosset et al. 2008" ~ 99999,
    
    Citation == "Jensen et al. 2010" &
      LocationName != "Svalbard, Norway" ~ HostsSampled,
    Citation == "Jensen et al. 2010" ~ 99999,
    
    Citation == "Jessup et al 1993b" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Jessup et al 1993b" ~ 99999,
    
    Citation == "Jordan et al 2003" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Jordan et al 2003" ~ 99999,
    
    Citation == "Kapel and Nansen 1996" &
      LocationName != "Greenland" ~ HostsSampled,
    Citation == "Kapel and Nansen 1996" ~ 99999,
    
    Citation == "Katz et al. 2007" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Katz et al. 2007" ~ 99999,
    
    Citation == "Kelley and Horner 2008" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Kelley and Horner 2008" ~ 99999,
    
    Citation == "King et al. 2004" &
      HostAge == "All" ~ HostsSampled,
    Citation == "King et al. 2004" ~ 99999,
    
    Citation == "Kiraly and Egri 2007" &
      LocationName != "Hungary" ~ HostsSampled,
    Citation == "Kiraly and Egri 2007" ~ 99999,
    
    Citation == "Kollars et al 1997a" &
      (is.na(HostAge) | HostAge == "All") ~ HostsSampled,
    Citation == "Kollars et al 1997a" ~ 99999,
    
    Citation == "Krumm et al 2005" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Krumm et al 2005" ~ 99999,
    
    Citation == "Kutz et al 2001" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Kutz et al 2001" ~ 99999,
    
    Citation == "Laakkonen et al. 1998" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Laakkonen et al. 1998" ~ 99999,
    
    Citation == "Lepojev et al 1999" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Lepojev et al 1999" ~ 99999,
    
    Citation == "Losson et al. 2003" &
      (is.na(HostAge) & is.na(HostSex)) |
      (HostAge == "All" & HostSex == "All") ~ HostsSampled,
    Citation == "Losson et al. 2003" ~ 99999,
    
    Citation == "Malczewski et al. 1995" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Malczewski et al. 1995" ~ 99999,
    
    Citation == "Manas et al. 2005" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Manas et al. 2005" ~ 99999,
    
    Citation == "Marco et al 2008a" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Marco et al 2008a" ~ 99999,
    
    Citation == "Marco et al 2008b" &
      HostsSampled == 145 ~ HostsSampled,
    Citation == "Marco et al 2008b" ~ 99999,
    
    Citation == "Matsuura et al 2007b" &
      LocationName == "Hyogo, Japan" ~ HostsSampled,
    Citation == "Matsuura et al 2007b" ~ 99999,
    
    Citation == "McFadden et al. 2005" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "McFadden et al. 2005" ~ 99999,
    
    Citation == "Michel 1993" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Michel 1993" ~ 99999,
    
    Citation == "Millan and Rodriguez 2009" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Millan and Rodriguez 2009" ~ 99999,
    
    Citation == "Mitchell et al. 1999" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Mitchell et al. 1999" ~ 99999,
    
    Citation == "Molia et al. 2004" &
      LocationName == "Masai Mara National Park, Nairobi National Park, Ngorongoro crater,  Serengeti National Park and Namibia" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Molia et al. 2004" &
      HostCorrectedName == "Panthera leo" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Molia et al. 2004"  ~ 99999,
    
    Citation == "Monello and Gompper 2007" &
      HostSex == "All" & HostAge == "All" &
      Prevalence > 0.9 ~ HostsSampled,
    Citation == "Monello and Gompper 2007" ~ 99999,
    
    Citation == "Monello and Gompper 2009" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Monello and Gompper 2009" ~ 99999,
    
    Citation == "Mulvey et al 1991" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Mulvey et al 1991" ~ 99999,
    
    Citation == "Negovetich et al 2006" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Negovetich et al 2006" ~ 99999,
    
    Citation == "Nelson et al. 2003" &
      LocationName != "Illinois, USA" ~ HostsSampled,
    Citation == "Nelson et al. 2003" ~ 99999,
    
    Citation == "Newman et al. 2004" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Newman et al. 2004" ~ 99999,
        
    Citation == "Nyamsuren et al 2006" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Nyamsuren et al 2006" ~ 99999,
    
    Citation == "Oksanen et al. 2009" &
      (HostSex == "All" | is.na(HostSex)) & HostAge == "All" ~ HostsSampled,
    Citation == "Oksanen et al. 2009" ~ 99999,
    
    Citation == "Richards et al. 1995" &
      HostSex == "Male" & HostsSampled == 229 ~ 292,
    Citation == "Richards et al. 1995" ~ HostsSampled,
    
    Citation == "Richardson and Demarais 1992" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Richardson and Demarais 1992" ~ 99999,
    
    Citation == "Roddie et al. 2008" &
      HostSex == "All" & HostAge == "All" &
      Prevalence > 0.6 ~ HostsSampled,
    Citation == "Roddie et al. 2008" ~ 99999,
    
    Citation == "Rodwell et al 2001b" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Rodwell et al 2001b" ~ 99999,
    
    Citation == "Roelke-Parker et al. 1996" &
      HostsSampled == (77|34) ~ HostsSampled,
    Citation == "Roelke-Parker et al. 1996" ~ 99999,
    
    Citation == "Rotstein et al. 2000" &
      HostsSampled > 9 ~ HostsSampled,
    Citation == "Rotstein et al. 2000" ~ 99999,
    
    Citation == "Saeed and Kapel 2006" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Saeed and Kapel 2006" ~ 99999,
    
    Citation == "Saeed et al. 2006" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Saeed et al. 2006" ~ 99999,
    
    Citation == "Santin-Duran et al 2004" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Santin-Duran et al 2004" ~ 99999,
    
    Citation == "Santin-Duran et al 2004" &
      ParasiteCorrectedName == "Ostertagia drozdzi" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Santin-Duran et al 2004" &
      ParasiteCorrectedName %in% c("Ostertagia leptospicularis",
                                   "Spiculopteragia asymmetrica",
                                   "Trichostrongylus axei") &
      HostSex != "All" ~ HostsSampled,
    Citation == "Santin-Duran et al 2004" ~ 99999,
    
    Citation == "Sidorovich and Anisimova 1997" &
      LocationName != "Belarus" ~ HostsSampled,
    Citation == "Sidorovich and Anisimova 1997" ~ 99999,
    
    Citation == "Simpson et al. 2009" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Simpson et al. 2009" ~ 99999,
    
    Citation == "Stubblefield et al 1987" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Stubblefield et al 1987" ~ 99999,
    
    Citation == "Stuve 1987" &
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Stuve 1987" ~ 99999,
    
    Citation == "Takumi et al. 2008" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Takumi et al. 2008" ~ 99999,
    
    Citation == "Thiede et al 2002" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Thiede et al 2002" ~ 99999,
    
    Citation == "Thomas et al. 2008" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Thomas et al. 2008" ~ 99999,
    
    Citation == "Tryland et al 2004" &
      HostAge != "Juvenile" ~ HostsSampled,
    Citation == "Tryland et al 2004" ~ 99999,
    
    Citation == "Valcarcel et al. 2004" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Valcarcel et al. 2004" ~ 99999,
    
    Citation == "van der Giessen et al. 1999" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "van der Giessen et al. 1999" ~ 99999,
    
    Citation == "Vervaeke et al. 2003" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Vervaeke et al. 2003" ~ 99999,
    
    Citation == "Vikoren et al 2004" &
      LocationName != "Norway" ~ HostsSampled,
    Citation == "Vikoren et al 2004" &
      HostCorrectedName == "Rangifer tarandus" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Vikoren et al 2004" ~ 99999,
    
    Citation == "Wacker et al. 1999" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Wacker et al. 1999" ~ 99999,
    
    Citation == "Waid et al 1985" & 
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Waid et al 1985" ~ 99999,
    
    Citation == "Waldrup et al 1992" & 
      is.na(HostAge) & is.na(HostSex) ~ HostsSampled,
    Citation == "Waldrup et al 1992" ~ 99999,
    
    Citation == "Wixsom et al. 1991" &
      HostsSampled == (293|85) ~ HostsSampled,
    Citation == "Wixsom et al. 1991" ~ 99999,
    
    Citation == "Yamamoto et al. 1998" &
      LocationName != "California state, USA" ~ HostsSampled,
    Citation == "Yamamoto et al. 1998" ~ 99999,
    
    Citation == "Yeitz et al. 2009" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Yeitz et al. 2009" ~ 99999,
    
    Citation == "Yimam et al. 2002" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Yimam et al. 2002" ~ 99999,
    
    Citation == "Zanella et al 2008" &
      (HostSex == "All" & HostAge == "All") |
      (is.na(HostSex) & is.na(HostAge)) ~ HostsSampled,
    Citation == "Zanella et al 2008" ~ 99999,
    
    Citation == "Zanella et al. 2008b" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Zanella et al. 2008b" ~ 99999,
    
    Citation == "Zarnke and Evans 1989" &
      LocationName != "Alaska" ~ HostsSampled,
    Citation == "Zarnke and Evans 1989" ~ 99999,
    
    Citation == "Zarnke et al. 1995" &
      LocationName != "Alaska" ~ HostsSampled,
    Citation == "Zarnke et al. 1995" ~ 99999,
    
    Citation == "Zarnke et al. 2004b" &
      HostAge == "All" ~ HostsSampled,
    Citation == "Zarnke et al. 2004b" ~ 99999,
    
    TRUE ~ NA
  ))

Citations <- GMPD_Data_Filter %>% filter(is.na(SampleSize)) %>% 
  pull(Citation)
# GMPD_Data %>% filter(Citation %in% Citations) %>% View()

GMPD_Data_Filter <- GMPD_Data_Filter %>% 
  mutate(SampleSize = case_when(
    !is.na(SampleSize) ~ SampleSize,
    # Kept
    Citation %in% c("Alexander et al. 1995",
                    "Casulli et al. 2005",
                    "Criado-Fornelio et al. 2000",
                    "Delahay et al. 1998",
                    "Hackett and Walters 1980",
                    "Ikeda et al. 1999",
                    "Khan et al. 1991",
                    "King 1991",
                    "Modric and Huber 1993",
                    "Moro et al. 1998",
                    "Pung et al. 1996",
                    "Ramsauer et al. 2007",
                    "Salkeld et al. 2007",
                    "Shimalov and Shimalov 2002",
                    "Weber 1991",
                    "Willi et al. 2007") ~ HostsSampled,
    
    Citation %in% c("Alexander et al. 2010",
                    "de Lisle et al. 2008",
                    "Truyen et al. 1998") ~ NumSamples,
    
    # Rejected
    Citation %in% c("Bender and Hall 1996",
                    "Fuglei et al. 2008",
                    "Shimalov and Shimalov 2002a",
                    "Stieger et al. 2002") ~ 99999,
    
    # Addressed
    Citation == "Alados et al 1996" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Alados et al 1996" ~ 99999,
    
    Citation == "Brown et al. 2010" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Brown et al. 2010" ~ 99999,

    Citation == "Carlson and Nielson 1985" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Carlson and Nielson 1985" ~ 99999,   
    
    Citation == "Gallivan et al 1995" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Gallivan et al 1995" ~ 99999,
    
    Citation == "Hamilton et al. 2005" &
      LocationName != "United Kingdom" ~ HostsSampled,
    Citation == "Hamilton et al. 2005" ~ 99999,
    
    Citation == "Kelly et al 2008" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Kelly et al 2008" ~ 99999, 

    Citation == "Kelly et al. 2010" &
      LocationName != "Ireland" ~ HostsSampled,
    Citation == "Kelly et al. 2010" ~ 99999,
    
    Citation == "Magnaval et al. 2004" &
      !is.na(HostsSampled) ~ HostsSampled,
    Citation == "Magnaval et al. 2004" ~ 99999,
    
    Citation == "Matthee et al 1997" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Matthee et al 1997" ~ 99999,
    
    Citation == "Matthee et al 1998" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Matthee et al 1998" ~ 99999,
    
    Citation == "Morfeld et al 2001" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Morfeld et al 2001" ~ 99999,
    
    Citation == "Radwan et al. 2009" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Radwan et al. 2009" ~ 99999, 
    
    Citation == "Reinecke et al 1988" &
      is.na(HostSex) ~ HostsSampled,
    Citation == "Reinecke et al 1988" ~ 99999,
    
    Citation == "Ryser-Degiorgis et al. 2006" &
      LocationName != "North & Central Sweden" ~ HostsSampled,
    Citation == "Ryser-Degiorgis et al. 2006" ~ 99999,
    
    Citation == "Sangster et al. 2007" &
      Prevalence > 0.06 ~ HostsSampled,
    Citation == "Sangster et al. 2007" ~ 99999,
    
    Citation == "Stefanidesova et al 2008" &
      LocationName != "Slovakia" ~ HostsSampled,
    Citation == "Stefanidesova et al 2008" ~ 99999,
    
    Citation == "Verbisck-Bucker et al 2008" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Verbisck-Bucker et al 2008" ~ 99999,
    
    Citation == "Vervaeke et al. 2005" &
      HostSex == "All" & HostAge == "All" ~ HostsSampled,
    Citation == "Vervaeke et al. 2005" ~ 99999,
    
    Citation == "Zabiega 1996" &
      HostSex == "All" ~ HostsSampled,
    Citation == "Zabiega 1996" ~ 99999,
    
    TRUE ~ NA
  ))

GMPD_Data <- GMPD_Data_Filter %>% 
  filter(SampleSize != 99999)
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 6530 128

rm(GMPD_Data_Filter)

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

GMPD_Data <- GMPD_Data %>%
  dplyr::group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  dplyr::filter(n()>1) %>% 
  dplyr::ungroup()
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 5714 and 108

# Native and non-native ####
## Removing non-native polygons ####

# reducing to relevant subsets of IUCN data to save time and space 
Hostlist <- sort(unique(GMPD_Data$HostCorrectedName))
sf::sf_use_s2(FALSE) # make valid only works fully if geometry is planar
IUCN_Data <- IUCN_Mammals[IUCN_Mammals$sci_name %in% Hostlist, ] %>% 
  sf::st_make_valid()

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

## Species corrections and removals
# Lynx lynx
# some Lynx canadensis samples appear to be incorrectly recorded as Lynx lynx
GMPD_Data <- sf::st_as_sf(GMPD_Data,
                          coords = c("Longitude", "Latitude"),
                          crs = Projection_String)
tm_shape(GMPD_Data[GMPD_Data$HostCorrectedName == "Lynx lynx",]) +
  tm_dots("countrycode")

GMPD_Data <- GMPD_Data %>% 
  mutate(HostCorrectedName = case_when(HostCorrectedName == "Lynx lynx" &
                                         countrycode %in% c("CAN", "USA") ~ "Lynx canadensis",
                                       TRUE ~ HostCorrectedName))

# Urocyon littoralis
"Six distinct subspecies are recognized, one on each of the islands where they occur:

San Miguel Island Fox (Urocyon littoralis littoralis (Baird, 1858)), San Miguel Island,
Santa Rosa Island Fox (U. l. santarosae Grinnell & Linsdale, 1930), Santa Rosa Island,
Santa Cruz Island Fox (U. l. santacruzae Merriam, 1903), Santa Cruz Island,
Santa Catalina Island Fox (U. l. catalinae Merriam, 1903), Santa Caralina Island,
San Nicolas Island Fox (U. l. dickeyi Grinnell & Linsdale, 1930), San Nicolas Island, and
San Clemente Island Fox (U. l. clementae Merriam, 1903), San Clemente Island."
# Coonan, T., Ralls, K., Hudgens, B., Cypher, B. & Boser, C. 2013. Urocyon littoralis. The IUCN Red List of Threatened Species 2013: e.T22781A13985603. https://dx.doi.org/10.2305/IUCN.UK.2013-2.RLTS.T22781A13985603.en. Accessed on 07 March 2024.
# Population structure makes them not relevant to our analysis, particularly range position component


Native_DF <- IUCN_Data %>% 
  dplyr::select(sci_name, legend) %>% 
  sf::st_drop_geometry() %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(keep = case_when(str_detect(sci_name, "Cervus") & legend == "Extant & Introduced (resident)"  	    ~ TRUE,
                                 sci_name == "Lynx canadensis" & legend == "Presence Uncertain & Origin Uncertain"  ~ TRUE,
                                 str_detect(legend, "Introduced") 								                                  ~ FALSE,
                                 legend == "Extant & Origin Uncertain (resident)"						                        ~ FALSE,
                                 TRUE												                                                        ~ TRUE)) 

IUCN_Data <- Native_DF %>% 
  {dplyr::right_join(IUCN_Data, ., by = c("sci_name", "legend"), relationship = "many-to-one")}

IUCN_Native_Data <- Native_DF %>% 
  dplyr::filter(keep) %>% 
  dplyr::select(-keep) %>% 
  {dplyr::right_join(IUCN_Data, ., by = c("sci_name", "legend"), relationship = "many-to-one")}


# GMPD_plots_native_base_02 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_base_02) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_base_02.pdf'), width = 10, height = 7)
# GMPD_plots_native_base_02
# dev.off()

## Restricting by native IUCN polygon ####

GMPD_Data_res_all <- lapply(Hostlist, function(host, point.data = GMPD_Data, 
                                               range.polygon = IUCN_Native_Data, 
                                               buff = 0){
  range.polygon <- range.polygon[range.polygon$sci_name == host, ]
  point.data <- point.data[point.data$HostCorrectedName == host, ]
  out <- pip_test(point.data, range.polygon, buff)
})

GMPD_Data_res_all <- GMPD_Data_res_all[lapply(GMPD_Data_res_all, function(x)nrow(x) != 0) == TRUE]
GMPD_Data_res_all <- do.call(rbind, GMPD_Data_res_all)

nrow(GMPD_Data_res_all); length(unique(GMPD_Data_res_all$HostCorrectedName)) # 4746 and 100

# GMPD_plots_native_restricted_03 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_res_all, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_restricted_03) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_restricted_03.pdf'), width = 10, height = 7)
# GMPD_plots_native_restricted_03
# dev.off()

## Cleaning by native IUCN polygon ####

# lapply(Hostlist, plot_native)

# Canis aureus
"Recent studies based on mtDNA and morphology have shown that 'Golden Jackals' in Africa are larger in size than those from 
Eurasia and are actually more closely related to the Grey Wolf Canis lupus. African animals hence represent a previously 
overlooked distinct species, the African Wolf, Canis lupaster (see Rueness et al. 2011, Gaubert et al. 2012, Koepfli et al. 
2015, Viranta et al. 2017). However, the putative presence of Golden Jackal in the Sinai Peninsula of Egypt remains unclear 
(see Gaubert et al. 2012, Viranta et al. 2017)." 
# Hoffmann, M., Arnold, J., Duckworth, J.W., Jhala, Y., Kamler, J.F. & Krofel, M. 2018. Canis aureus (errata version published in 2020). The IUCN Red List of Threatened Species 2018: e.T118264161A163507876. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T118264161A163507876.en. Accessed on 07 March 2024.

# Ovis ammon 
# Sample locations don't correspond at all to Ovis ammon range

# Phacochoerus aethiopicus
# Some suspicious GMPD samples
# odd ones are probably introduced for hunting

# Rangifer tarandus
# Excluding outlying points in North Norway because of overlap with semi-domestic herds (see Sami wiki page)
# Although caribou in North America are hunted, they remain undomesticated


# come back to after gbif:
c("Antilocapra americana",
  "Bison bison",
  "Bison bonasus",
  "Equus quagga,",
  "Hippotragus niger",
  "Lynx pardinus", 
  "Lynx rufus")

Native_DF <- Native_DF %>% 
  mutate(gmpd_buffer = case_when(sci_name %in% c("Genetta genetta",
                                                 "Mustela erminea",
                                                 "Nyctereutes procyonoides") & !keep ~ 0,
                                 sci_name %in% c("Aepyceros melampus",
                                                 "Ovibos moschatus",
                                                 "Rupicapra rupicapra")      & !keep ~ 0.5,
                                 sci_name %in% c("Vulpes vulpes")            & !keep ~ 2, 
                                 sci_name %in% c("Canis aureus",
                                                 "Neovison vison",
                                                 "Ovis ammon")               & keep  ~ 0,
                                 sci_name %in% c("Rangifer tarandus",
                                                 "Cervus elaphus")           & keep  ~ 1.5,
                                 sci_name %in% c("Phacochoerus aethiopicus",
                                                 "Procyon lotor")            & keep  ~ 6.5,
                                 TRUE ~ NA_real_))


GMPD_Data_cln_all <- lapply(Hostlist, function(host, point.data = GMPD_Data, 
                                               range.polygon = IUCN_Data, 
                                               buff = 0){
  Native_DF <- Native_DF[Native_DF$sci_name == host,]
  
  # Start with ones that need no buffer
  if(all(is.na(Native_DF$gmpd_buffer))) {
    out <- point.data[point.data$HostCorrectedName == host, ]
    
  # Then do ones that have no buffer for kept but do for removed
  } else if(all(is.na(Native_DF[Native_DF$keep, "gmpd_buffer"]))){
    
    # filter range polygon to ones I'm removing
    range.polygon <- range.polygon[range.polygon$sci_name == host &
                                     !range.polygon$keep, ]
    point.data <- point.data[point.data$HostCorrectedName == host, ]
    
    # extract buffer
    buff <- unique(Native_DF[!Native_DF$keep, "gmpd_buffer"])[[1]]
    
    # buffer the polygon
    range.polygon <- sf::st_buffer(range.polygon, buff)
    
    # get the intersections
    kept.indices <- sf::st_intersects(point.data, range.polygon)
    
    # Keep only the points that have no intersections
    out <- point.data[lengths(kept.indices) == 0, ]
    
  } else {
    # get the full range polygon and the point data for species
    range.polygon <- range.polygon[range.polygon$sci_name == host, ]
    point.data <- point.data[point.data$HostCorrectedName == host, ]
    
    # get the buffer
    buff <- unique(Native_DF[Native_DF$keep, "gmpd_buffer"])[[1]]
    
    # TODO replace pip_test
    # use pip_test to get the points
    out <- pip_test(point.data, range.polygon, buff)
  }
  
  return(out)
})

GMPD_Data_cln_all <- GMPD_Data_cln_all[lapply(GMPD_Data_cln_all, function(x)nrow(x) != 0) == TRUE]
GMPD_Data_cln_all <- do.call(rbind, GMPD_Data_cln_all)

GMPD_Data_cln_all <- GMPD_Data_cln_all %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data_cln_all); length(unique(GMPD_Data_cln_all$HostCorrectedName)) # 5367 and 106

# GMPD_plots_native_clean_03 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_cln, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_clean_03) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_clean_03.pdf'), width = 10, height = 7)
# GMPD_plots_native_clean_03
# dev.off()


# Adjusting subspecies ####
## Restricting by native IUCN polygon by subspecies group ####

# https://github.com/r-spatial/sf/wiki/Migrating
# TODO: redo
# source(here::here("01.1Subspecies polygons.R"))

# TODO: redo

# buffer_temp <- unique(GMPD_Subgroups[,c("HostCorrectedName", "subgroup")])
# buffer_temp$buff <- 0
# 
# GMPD_Data_res_sub <- lapply(Hostlist, function(host, dat = GMPD_Spatial, range.polygon = IUCN_Native_Data, buff = buffer_temp){
#   range.polygon <- range.polygon[range.polygon$sci_name == host,]
#   buff <- buff[buff$HostCorrectedName == host,]
#   out <- pip_test(host, dat, range.polygon, buff)
# })
# 
# GMPD_Data_res_sub <- GMPD_Data_res_sub[which(lapply(GMPD_Data_res_sub, is.null) == FALSE)]
# GMPD_Data_res_sub <- do.call(raster::bind, GMPD_Data_res_sub)
# GMPD_Data_res_sub <- as.data.frame(GMPD_Data_res_sub)
# 
# GMPD_Data_res_sub <- GMPD_Data_res_sub %>%
#   group_by(ParasiteCorrectedName, HostCorrectedName) %>%
#   filter(n() > 1) %>% 
#   ungroup()
# nrow(GMPD_Data_res_sub); length(unique(GMPD_Data_res_sub$HostCorrectedName)) # 7113 and 105


# Restricting by proximity ####################################################################################
## IUCN restricted subgroup ####
# Res_Temp <- restrict_decimal(GMPD_Data_res_sub, subsp = TRUE)
# Res_Temp <- Res_Temp[Res_Temp$enough, ]
# GMPD_Data_res_sub <- GMPD_Data_res_sub[GMPD_Data_res_sub$HostCorrectedName %in% Res_Temp$HostCorrectedName &
#                                          GMPD_Data_res_sub$subgroup %in% Res_Temp$subgroup,]
# GMPD_Data_res_sub <- GMPD_Data_res_sub %>%
#   group_by(HostCorrectedName, subgroup) %>%
#   filter(n() > 1) %>%
#   group_by(HostCorrectedName, ParasiteCorrectedName) %>%
#   filter(n() > 1) %>%
#   ungroup() 
# nrow(GMPD_Data_res_sub); length(unique(GMPD_Data_res_sub$HostCorrectedName)) # 6958 and 90

## IUCN restricted species ####
GMPD_Data_res_all <- cbind(sf::st_drop_geometry(GMPD_Data_res_all), 
                           data.frame(sf::st_coordinates(GMPD_Data_res_all))) %>% 
  dplyr::rename(Longitude = X, Latitude = Y)
Res_Temp <- restrict_decimal(GMPD_Data_res_all)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_res_all <- GMPD_Data_res_all[GMPD_Data_res_all$HostCorrectedName %in% Res_Temp$HostCorrectedName,]
GMPD_Data_res_all <- GMPD_Data_res_all %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  filter(n() > 1) %>%
  ungroup() 
nrow(GMPD_Data_res_all); length(unique(GMPD_Data_res_all$HostCorrectedName)) # 4538 and 84

rm(Res_Temp)

## IUCN cleaned subgroup ####
# Res_Temp <- restrict_decimal(GMPD_Data_cln, subsp = TRUE)
# Res_Temp <- Res_Temp[Res_Temp$enough, ]
# GMPD_Data_cln_sub <- GMPD_Data_cln[GMPD_Data_cln$HostCorrectedName %in% Res_Temp$HostCorrectedName &
#                                  GMPD_Data_cln$subgroup %in% Res_Temp$subgroup,]
# GMPD_Data_cln_sub <- GMPD_Data_cln_sub %>%
#   group_by(HostCorrectedName, subgroup) %>%
#   filter(n() > 1) %>%
#   group_by(HostCorrectedName, ParasiteCorrectedName) %>%
#   filter(n() > 1) %>%
#   ungroup()
# nrow(GMPD_Data_cln_sub); length(unique(GMPD_Data_cln_sub$HostCorrectedName)) # 7995 and 107

## IUCN cleaned species ####
GMPD_Data_cln_all <- cbind(sf::st_drop_geometry(GMPD_Data_cln_all), 
                           data.frame(sf::st_coordinates(GMPD_Data_cln_all))) %>% 
  dplyr::rename(Longitude = X, Latitude = Y)
Res_Temp <- restrict_decimal(GMPD_Data_cln_all)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_cln_all <- GMPD_Data_cln_all[GMPD_Data_cln_all$HostCorrectedName %in% Res_Temp$HostCorrectedName,]
GMPD_Data_cln_all <- GMPD_Data_cln_all %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  filter(n() > 1) %>%
  ungroup()
nrow(GMPD_Data_cln_all); length(unique(GMPD_Data_cln_all$HostCorrectedName)) # 5276 and 96

rm(Res_Temp)

# TODO: sort out afterthought
# merging gmpd into one (afterthought) ####
GMPD_Data_cln_all$CleanAll <- TRUE
# GMPD_Data_cln_sub$CleanSub <- TRUE
GMPD_Data_res_all$RestrAll <- TRUE
# GMPD_Data_res_sub$RestrSub <- TRUE

GMPD_Data <- merge(GMPD_Data_cln_all, GMPD_Data_res_all, 
                   by = intersect(names(GMPD_Data_cln_all), names(GMPD_Data_res_all)),
                   all = TRUE)

# GMPD_Data <- merge(GMPD_Data, GMPD_Data_cln_sub, 
#                    by = intersect(names(GMPD_Data), names(GMPD_Data_cln_sub)),
#                    all = TRUE)
# 
# GMPD_Data <- merge(GMPD_Data, GMPD_Data_res_sub, 
#                    by = intersect(names(GMPD_Data), names(GMPD_Data_res_sub)),
#                    all = TRUE)

GMPD_Data <- GMPD_Data %>%
  # mutate(CleanSub = case_when(is.na(CleanSub) ~  FALSE, TRUE ~ CleanSub)) %>%
  # mutate(RestrSub = case_when(is.na(RestrSub) ~  FALSE, TRUE ~ RestrSub)) %>% 
  mutate(RestrAll = case_when(is.na(RestrAll) ~  FALSE, TRUE ~ RestrAll)) 

# Write files #################################################################################################
# narrowing down IUCN_Mammals to a more manageable size
# IUCN_Orders <- IUCN_Mammals[IUCN_Mammals$order_ %in% unique(IUCN_Native_Data$order_), ]
# rm(IUCN_Mammals)
# 
write.csv(GMPD_Data, file = here::here("Data/Data back ups/GMPD_Data_01.csv"), row.names = FALSE)

write.csv(GMPD_Location_Data, file = here::here("Data/Data back ups/GMPD_Location_Data_01.csv"), row.names = FALSE)
write.csv(Native_DF, file = here::here("Data/Data back ups/Native_DF_01.csv"), row.names = FALSE)

saveRDS(IUCN_Native_Data, file = here::here("Data/Data back ups/IUCN_Native_Data_01")) # too large to commit
saveRDS(IUCN_Data, file = here::here("Data/Data back ups/IUCN_Data_01")) # too large to commit
# saveRDS(IUCN_Orders, file = here::here("Data/Data back ups/IUCN_Orders_01")) # too large to commit 

