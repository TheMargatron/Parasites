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
library(CoordinateCleaner)  # cleaning geographic data
library(geosphere)          # calculating distances
library(here)               #
library(maps)               # iso 3166 country codes and mapnames
library(maptools)
library(rnaturalearth)      # river data
library(rnaturalearthdata)  # countries data
library(raster)             # 
library(rgdal)              # read shapefiles
library(rgeos)              # Just for gBuffer
library(stringr)            #
library(tidyverse)          # beware of conflicts (mainly with raster)
library(tmap)

source(here::here("Functions.R"))

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
nrow(GMPD_Raw_Data); length(unique(GMPD_Raw_Data$HostCorrectedName)) #Beginning with 24323 rows and 462 hosts

IUCN_Mammals <- readOGR(here::here("Data/IUCN"), "MAMMALS") # plenty time to make a cup of tea

River_Data50 <- ne_load(scale = 50,
                        type = "rivers_lake_centerlines",
                        category = "physical",
                        destdir = here::here("Data/Extras/ne_rivers"))

sf::sf_use_s2(FALSE) # For "invalid spherical geometry" errors
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
  filter(HostCorrectedName != "Ovis aries" &
           HostCorrectedName != "Bos frontalis" &
           HostCorrectedName != "no binomial name" &
           HostCorrectedName != "Dama dama" &
           HostCorrectedName != "Diceros bicornis" &
           HostCorrectedName != "Ceratotherium simum") %>%
  filter(!str_detect(HostReportedName, " and ")) %>%
  filter(!is.na(Prevalence)) %>%
  filter(HostEnvironment != "marine") %>%
  filter(NativeRange != "No" &
           !is.na(NativeRange)) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  filter(!is.na(NumSamples) | !is.na(HostsSampled))
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 12061 rows of 201 hosts

GMPD_Data <- GMPD_Data %>%
  separate(HostReportedName, c("HostReportedGenus", "HostReportedSpecies", "HostReportedSubspecies"), remove = FALSE) 
GMPD_Data[1045, "HostReportedName"]

subsp_info <- unique(GMPD_Data[which(GMPD_Data$HostReportedName != GMPD_Data$HostCorrectedName), c("HostReportedName", 
                                                                                          "HostReportedGenus", 
                                                                                          "HostReportedSpecies",
                                                                                          "HostReportedSubspecies",
                                                                                          "HostCorrectedName")])

# Adding subspecies info and correcting HostCorrectedName according to IUCN
# Taxonomic justifications in Subspecies Info excel
GMPD_Data <- GMPD_Data %>%
  mutate(HostReportedSubspecies = case_when(HostReportedName == "Alcelaphus cokii"           ~ "cokii", 
                                            HostCorrectedName == "Alcelaphus lichtensteinii" ~ "lichtensteinii", 
                                            HostReportedName == "Canis latrans Say"          ~ NA_character_, 
                                            HostReportedName == "Black-back Jackal"          ~ NA_character_, 
                                            HostReportedName == "Capra ibex ibex"            ~ NA_character_, 
                                            HostReportedName == "Capra i. ibex"              ~ NA_character_, 
                                            HostReportedName == "Cervus elaphus nelsoni"     ~ "canadensis", 
                                            HostReportedName == "Cervus nippon centralis"    ~ "nippon", 
                                            HostReportedName == "Damaliscus korrigum"        ~ "korrigum", 
                                            HostReportedName == "Damaliscus dorcas dorcas"   ~ NA_character_, 
                                            HostReportedName == "Damaliscus pygargus dorcas" ~ NA_character_, 
                                            HostCorrectedName == "Equus burchellii"          ~ "burchellii", 
                                            str_detect(HostReportedName, "Felis libyca")     ~ "libyca", 
                                            HostReportedName == "Felis silvestris gordoni"   ~ "libyca", 
                                            HostReportedName == "Giraffa reticulata"         ~ "reticulata", 
                                            HostReportedName == "Hyaena hyaena dubbah"       ~ NA_character_, 
                                            HostReportedName == "Kobus defassa"              ~ "defassa", 
                                            HostReportedName == "Lynx rufus floridanus"      ~ "rufus", 
                                            HostReportedName == "Martes caurina"             ~ "caurina", 
                                            HostReportedName == "Meles meles anakuma"        ~ NA_character_, 
                                            HostReportedName == "Melogale moschata subauantiaca" ~ "subaurantiaca", 
                                            HostReportedName == "Mustela itatsi sho"         ~ NA_character_, 
                                            HostReportedName == "Neovison vison mink"        ~ NA_character_, 
                                            HostReportedName == "Oryx gazella gazella"       ~ NA_character_, 
                                            HostReportedName == "Ourebia ourebi cottoni"     ~ NA_character_, 
                                            HostReportedName == "Ovis canadensis cremnobates" ~ "nelsoni", 
                                            HostReportedName == "Ovis canadensis mexicana"   ~ "nelsoni", 
                                            HostReportedName == "Felis leo senegalensis"     ~ "leo", 
                                            HostReportedName == "Panthera pardus saxicolor"  ~ "tulliana", 
                                            HostReportedName == "Puma concolor coryi"        ~ "couguar", 
                                            HostReportedName == "Puma concolor stanleyana"   ~ "couguar", 
                                            HostReportedName == "Felis concolor coryi"       ~ "couguar", 
                                            HostReportedName == "Felis concolor vancouverensis" ~ "couguar", 
                                            HostReportedName == "Spilogale gracilis amphiala" ~ "amphialus", 
                                            HostReportedName == "Urocyon cinereoargenteus texensis" ~ "scottii", 
                                            HostReportedName == "Ursus americanus pallas"    ~ NA_character_, 
                                            HostReportedName == "Ursus arctos marsicanus"    ~ "arctos", 
                                            HostReportedName == "Viverra civetta schwartzi"  ~ "schwarzi", 
                                            HostReportedName == "Vulpes vulpes schrencki"    ~ "schrenckii",
                                            TRUE ~ HostReportedSubspecies)) %>%
  mutate(HostCorrectedName = case_when(HostCorrectedName == "Alcelaphus lichtensteinii" ~ "Alcelaphus buselaphus", 
                                       HostCorrectedName == "Alces americanus"          ~ "Alces alces", 
                                       HostReportedName == "Canis rufus"                ~ "Canis rufus", 
                                       HostReportedName == "Cervus elaphus nannodes"    ~ "Cervus canadensis", 
                                       HostReportedName == "Cervus elaphus nelsoni"     ~ "Cervus canadensis", 
                                       HostReportedName == "Cervus elaphus roosevelti"  ~ "Cervus canadensis", 
                                       HostReportedName == "Cervus elaphus canadensis"  ~ "Cervus canadensis", 
                                       HostCorrectedName == "Equus burchellii"          ~ "Equus quagga", 
                                       HostCorrectedName == "Felis manul"               ~ "Otocolobus manul", 
                                       HostCorrectedName == "Lama glama"                ~ HostReportedName, 
                                       HostCorrectedName == "Leopardus pajeros"         ~ "Leopardus colocolo", 
                                       HostReportedName == "Meles meles anakuma"        ~ "Meles anakuma", 
                                       HostCorrectedName == "Neotragus moschatus"       ~ "Nesotragus moschatus",     # IUCN name differs
                                       HostReportedName == "Putorius eversmanni"        ~ "Mustela eversmanni", 
                                       HostCorrectedName == "Puma yagouaroundi"         ~ "Herpailurus yagouaroundi", 
                                       HostCorrectedName == "Taurotragus oryx"          ~ "Tragelaphus oryx",         # IUCN name differs
                                       str_detect(HostReportedName, "Viverra civetta")  ~ "Civettictis civetta", 
                                       TRUE ~ HostCorrectedName)) %>%
  filter(HostReportedName != "Ovis ammon musimon") 
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 12034 and 202

GMPD_Data <- GMPD_Data %>%
  dplyr::select(-ParasiteReportedName,
                -HostReportedName,
                -HostReportedGenus,
                -HostReportedSpecies,
                -HasBinomialName, 
                -NativeRange, 
                -Intensity, 
                -IntensityMeasure, 
                -SampleNotes) %>%   #not used
  distinct()                                                                                                            #remove duplicated rows
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 11976 and 202

GMPD_Data <- GMPD_Data %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n()>1) %>% 
  ungroup()
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 10125 and 137

## adjust prevalence data that has been reported as a percentage

GMPD_Data <- GMPD_Data %>%
  mutate(Prevalence = case_when(Prevalence >1 ~ Prevalence/100,
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
## rows where SamplingBasis was "samples" AND HostsSampled was NA are excluded

## rows which are missing either HostsSampled or NumSamples: get data from the other one

## Filling in HostAge and HostSex reduces duplicated data
### e.g. data inputted once with NA for HostSex and again identically except HostSex is reported accurately
### This doesn't get rid of rows where HostSex differs between rows but still contains a value

## Different sampling method on the same sample group

# SamplingBasis and HostsSampled 
GMPD_Data <- GMPD_Data %>%
  filter(!is.na(HostsSampled) | SamplingBasis != "Samples" | is.na(SamplingBasis))
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 9971 and 137

# host age and sex NA duplications
GMPD_Data <- GMPD_Data %>%
  group_by_at(vars(-HostAge)) %>%
  fill(HostAge, .direction = "updown") %>% 
  
  group_by_at(vars(-HostSex)) %>%
  fill(HostSex, .direction  = "updown") %>% 
  
  group_by_at(vars(-HostAge, -HostSex)) %>%
  fill(c(HostSex, HostAge), .direction = "updown") %>% 
  distinct() %>%
  ungroup()

# HostsSampled vs NumSamples
GMPD_Data <- GMPD_Data %>%
  mutate(HostsSampled = case_when(is.na(HostsSampled) ~ NumSamples,
                                  TRUE                ~ HostsSampled)) %>%
  mutate(NumSamples = case_when(is.na(NumSamples) ~ HostsSampled,
                                TRUE              ~ NumSamples)) %>%
  dplyr::select(-NumSamples)
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 9562 and 137

# differing sample method on the same sample group

GMPD_Data <- GMPD_Data %>% 
  group_by_at(vars(-SamplingType, -Prevalence)) %>%
  mutate(sample_temp = case_when(n() == 1 ~ "fine",
                                 length(unique(SamplingType)) == 1 ~ "fine",
                                 length(unique(Prevalence)) == 1 & !duplicated(Prevalence) ~ "fine",
                                 length(unique(Prevalence)) == 1 ~ "not fine",
                                 n() == length(unique(SamplingType)) & Prevalence == max(Prevalence) ~ "fine",
                                 n() == length(unique(SamplingType)) ~ "not fine",
                                 length(unique(Prevalence)) > 1 & length(unique(SamplingType)) > 1 & Prevalence == max(Prevalence) ~ "fine",
                                 length(unique(Prevalence)) > 1 & length(unique(SamplingType)) > 1 ~ "not fine")) %>%
  filter(sample_temp == "fine") %>%
  dplyr::select(-sample_temp) %>%
  ungroup()
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 9515 137

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


Country_Match <- str_c(Country_Match, collapse = "|")
State_Match <- str_c(state.name, collapse = "|")

# All GMPD location descriptions needing matched to a country, 1972 unique descriptions
GMPD_Location_Data <- GMPD_Data %>%
  dplyr::select(LocationName) %>%
  distinct()

# extracting any country names in descriptions
Country_Data_Temp <- full_join(GMPD_Location_Data %>%
                         rownames_to_column(),
                        GMPD_Location_Data$LocationName %>%
                         str_to_lower() %>%
                         str_extract_all(str_to_lower(Country_Match), simplify=TRUE) %>%
                         as.data.frame() %>%
                         rownames_to_column(),
                       by = "rowname")

# joining multiple countries into single strings
Country_Data_Temp <- Country_Data_Temp %>%
  mutate_if(is.factor, as.character) %>%
  mutate_at(vars(-rowname, -LocationName), list(~ na_if(., ""))) %>%
  unite("countries", names(Country_Data_Temp)[-c(1,2)], sep = "|", remove = TRUE, na.rm = TRUE)

# lacking country name in description or misspelled country name (while already having one correct country or state)
Country_Data_Temp <- Country_Data_Temp %>%
  mutate(countries = case_when(LocationName == "Masai Mara National Park, Nairobi National Park, Ngorongoro crater,  Serengeti National Park and Namibia" ~ "namibia|kenya|tanzania",
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
State_Data_Temp <- full_join(GMPD_Location_Data %>%
                        rownames_to_column(),
                      GMPD_Location_Data$LocationName %>%
                        str_to_lower() %>%
                        str_extract_all(str_to_lower(State_Match),simplify = TRUE) %>%
                        as.data.frame() %>%
                        rownames_to_column(),
                      by = "rowname")

# joining states into single string
State_Data_Temp <- State_Data_Temp %>%
  mutate_if(is.factor, as.character) %>%
  mutate_at(vars(-rowname, -LocationName), list(~ na_if(., ""))) %>%
  unite("states", names(State_Data_Temp)[-c(1,2)], sep = "|", remove = TRUE, na.rm = TRUE)

# merge country and state data, then correct missing names (grouped roughly by continent)
GMPD_Location_Data <- full_join(State_Data_Temp, Country_Data_Temp, by = c("rowname", "LocationName")) %>%
  mutate(countries = case_when(countries != "" ~ countries,
                               states    != "" ~ "usa",
                               str_detect(LocationName, regex("yellowstone|channel islands|eastern us| CA(?!.)|yosemite|nebrask|califonia|orange|purdue|marion|susitna|orgeon|tallahal", ignore_case = TRUE)) ~ "usa",
                               str_detect(LocationName, regex("ontario|saskatchewan|quebec|yukon|nova scotia|alberta|prince edward|british columbia|northwest territories|baffin|vancou|newfoundland|brunswick|hudson bay", ignore_case = TRUE)) ~ "canada",
                               str_detect(LocationName, regex("orkendalen",                       ignore_case = TRUE)) ~ "greenland",
                               str_detect(LocationName, "Arctic")                                                      ~ "usa|canada",
                               
                               str_detect(LocationName, regex("parana|leones|marajo|jequitinhonha|pantanal", ignore_case = TRUE)) ~ "brazil",
                               str_detect(LocationName, regex("mexican",                          ignore_case = TRUE)) ~ "mexico",
                               str_detect(LocationName, regex("kaa",                              ignore_case = TRUE)) ~ "bolivia",
                               
                               str_detect(LocationName, regex("scotland|england|wales|united kingdom|great britain|shire|hebrides", ignore_case = TRUE)) ~ "uk",
                               str_detect(LocationName, regex("brandenburg|berlin",               ignore_case = TRUE)) ~ "germany",
                               str_detect(LocationName, regex("reykjavik",                        ignore_case = TRUE)) ~ "iceland",
                               str_detect(LocationName, regex("noway|svalbard|barents",           ignore_case = TRUE)) ~ "norway",
                               str_detect(LocationName, regex("bialowie|polish|Puszcza",          ignore_case = TRUE)) ~ "poland",
                               str_detect(LocationName, regex("swiss",                            ignore_case = TRUE)) ~ "switzerland",
                               str_detect(LocationName, regex("copenhagen",                       ignore_case = TRUE)) ~ "denmark",
                               str_detect(LocationName, regex("zilina|slova",                     ignore_case = TRUE)) ~ "slovakia",
                               str_detect(LocationName, regex("budakeszi",                        ignore_case = TRUE)) ~ "hungary",
                               str_detect(LocationName, regex("moravia|mim",                      ignore_case = TRUE)) ~ "czech republic",
                               str_detect(LocationName, regex("italia|sondrio|brembana|belviso",  ignore_case = TRUE)) ~ "italy",
                               str_detect(LocationName, regex("meurthe|french|bauges|savoy",      ignore_case = TRUE)) ~ "france",
                               str_detect(LocationName, regex("sorbe|zaragoza|malaga|catalonia|sierras|madrid|jaen|pallars|aller", ignore_case = TRUE)) ~ "spain", #assumed spain for jaen (as opposed to peru) because of species
                               str_detect(LocationName, regex("dalmatia",                         ignore_case = TRUE)) ~ "croatia",
                               str_detect(LocationName, regex("vojvodina",                        ignore_case = TRUE)) ~ "serbia",
                               str_detect(LocationName, regex("danubian",                         ignore_case = TRUE)) ~ "romania",
                              
                               str_detect(LocationName, regex("alpine areas",                     ignore_case = TRUE)) ~ "italy|switzerland|france",
                               str_detect(LocationName, regex("pyrenees",                         ignore_case = TRUE)) ~ "spain|france",
                               
                               str_detect(LocationName, regex("cameroun|ngaoun",                  ignore_case = TRUE)) ~ "cameroon",
                               str_detect(LocationName, regex("kenia|masai|nairobi|jogi|bungoma", ignore_case = TRUE)) ~ "kenya",
                               str_detect(LocationName, regex("kruger|natal|skukuza|queenstown|eastern shores|karroid|KNP|rooiwal|rietvlei|potchefstroom|benfontein|pieter|transvaa|sabi|hluhluwe|kuruman|mbiyamiti|ntomeni|west coast national park|transkei", ignore_case = TRUE)) ~ "south africa",
                               str_detect(LocationName, regex("serengeti|ngorongoro|selous|temi|ruaha|kaisho", ignore_case = TRUE)) ~ "tanzania",
                               str_detect(LocationName, regex("bale|sidamo|urso",                 ignore_case = TRUE)) ~ "ethiopia",
                               str_detect(LocationName, regex("zimbawe|hippo|mana pools|buffalo range", ignore_case = TRUE)) ~ "zimbabwe",
                               str_detect(LocationName, regex("zaire",                            ignore_case = TRUE)) ~ "democratic republic of the congo",
                               str_detect(LocationName, regex("adiopodoume",                      ignore_case = TRUE)) ~ "ivory coast",
                               str_detect(LocationName, regex("etosha",                           ignore_case = TRUE)) ~ "namibia",
                               str_detect(LocationName, regex("ankole|koja|jie",                  ignore_case = TRUE)) ~ "uganda",
                               str_detect(LocationName, regex("umsalala",                         ignore_case = TRUE)) ~ "sudan",
                               str_detect(LocationName, regex("bandia|saboya",                    ignore_case = TRUE)) ~ "senegal",
                               str_detect(LocationName, regex("batie",                            ignore_case = TRUE)) ~ "burkina faso",
                               str_detect(LocationName, regex("ndoki|republic of the congo",      ignore_case = TRUE)) ~ "republic of congo",
                               
                               str_detect(LocationName, regex("kerguelen",                        ignore_case = TRUE)) ~ "french southern and antarctic lands",
                               
                               str_detect(LocationName, regex("new zeland|flagstaff",             ignore_case = TRUE)) ~ "new zealand",
                               str_detect(LocationName, regex("tokyo|kkaido|sobo|shimane|akita",  ignore_case = TRUE)) ~ "japan",
                               str_detect(LocationName, regex("rahasthan",                        ignore_case = TRUE)) ~ "india",
                               str_detect(LocationName, regex("karak",                            ignore_case = TRUE)) ~ "jordan",
                               TRUE ~ countries))

# split countries at "|" then pivot into single column 
# Number of rows will increase slightly because of location descriptions for which there are multiple countries. 
# They'll be removed again by coordinate cleaner
GMPD_Location_Data <- GMPD_Location_Data %>%
  separate(countries, into = c("A","B","C"), sep = "\\|") %>%
  pivot_longer(cols = c("A","B","C"), values_to = "mapname", values_drop_na = TRUE) %>%
  dplyr::select(-name, -states) %>%
  mutate(mapname = case_when(mapname == "uk"      ~ "uk(?!r)",
                             mapname == "norway"  ~ "norway(?!:bouvet|:svalbard|:jan mayen)",
                             mapname == "finland" ~ "finland(?!:aland)",
                             mapname == "china"   ~ "china(?!:hong kong|:macao)",
                             TRUE                 ~ mapname)) # ignore warning

GMPD_Location_Data <- iso3166 %>%
  dplyr::select(a3, mapname) %>%
  mutate(mapname = tolower(mapname)) %>%
  right_join(GMPD_Location_Data, by = "mapname") %>%
  rename(countrycode = a3)

GMPD_Data <- full_join(GMPD_Location_Data, GMPD_Data, by = "LocationName")

rm(State_Match, Country_Match)
rm(list = ls(pattern = "_Temp$"))

## CoordinateCleaner tests ####################################################################################

GMPD_Data <- GMPD_Data %>%
  clean_coordinates(lon = "Longitude",
                    lat = "Latitude",
                    species = "HostCorrectedName",
                    countries = "countrycode",
                    tests = c("capitals","centroids","institutions", "countries"),
                    value = "clean")
nrow(GMPD_Data); length(unique(GMPD_Data$HostCorrectedName)) # 8840 and 132

# Native and non-native ####
## Removing non-native polygons ####

# reducing to relevant subsets of IUCN data to save time and space 

IUCN_Data <- IUCN_Mammals[IUCN_Mammals$binomial %in% Hostlist, ]

# Correcting legend for two species

# Oreamnos americanus
O_americanus_temp <- IUCN_Data[IUCN_Data$binomial == "Oreamnos americanus",]
tm_shape(O_americanus_temp) + tm_polygons("legend")
tm_shape(O_americanus_temp) + tm_polygons("dist_comm") # aggregated with Chichagof Island

O_americanus_temp <- raster::disaggregate(O_americanus_temp)
tm_shape(O_americanus_temp[1,]) + tm_polygons()

O_americanus_temp@data[1, "legend"] <- "Extant & Introduced (resident)"
O_americanus_temp <- raster::aggregate(O_americanus_temp, by = names(O_americanus_temp))
tm_shape(O_americanus_temp) + tm_polygons("legend")

IUCN_Data <- IUCN_Data[IUCN_Data$binomial != "Oreamnos americanus",]
IUCN_Data <- raster::bind(IUCN_Data, O_americanus_temp)

# Ovibos moschatus
O_moschatus_temp <- IUCN_Data[IUCN_Data$binomial == "Ovibos moschatus",]
tm_shape(O_moschatus_temp) + tm_polygons("legend")
tm_shape(O_moschatus_temp) + tm_polygons("SHAPE_Area") # aggregated across Greenland 

O_moschatus_temp <- raster::disaggregate(O_moschatus_temp)
tm_shape(O_moschatus_temp[c(53, 55, 58),]) + tm_polygons()

O_moschatus_temp@data[c(53, 55, 58), "legend"] <- "Extant & Introduced (resident)"
O_moschatus_temp <- raster::aggregate(O_moschatus_temp, by = names(O_moschatus_temp))
tm_shape(O_moschatus_temp) + tm_polygons("legend")

IUCN_Data <- IUCN_Data[IUCN_Data$binomial != "Ovibos moschatus",]
IUCN_Data <- raster::bind(IUCN_Data, O_moschatus_temp)

rm(list = ls(pattern = "_temp$"))

# Plotting with full polygons before restricting

Legend_Text <- sort(unique(IUCN_Data$legend))
# GMPD_plots_01 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Data, plot_type = "iucn")
# names(GMPD_plots_01) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_01.pdf'), width = 10, height = 7)
# GMPD_plots_01
# dev.off()

# removing unwanted polygons

Native_DF <- lapply(Hostlist, function(host) {
  host.levels <- unique(IUCN_Data[IUCN_Data$binomial == host, ]$legend)
  data.frame(HostCorrectedName = rep(host, length(host.levels)), 
             Status = host.levels)
}) %>%
  bind_rows() %>%
  mutate(keep = case_when(str_detect(HostCorrectedName, "Cervus") & Status == "Extant & Introduced (resident)"  	    ~ TRUE,
                          HostCorrectedName == "Lynx canadensis" & Status == "Presence Uncertain & Origin Uncertain" 	~ TRUE,
                          str_detect(Status, "Introduced") 								                                            ~ FALSE,
                          Status == "Extant & Origin Uncertain (resident)"						                                ~ FALSE,
                          TRUE												                                                                ~ TRUE))

IUCN_Native_Data <- raster::bind(lapply(Hostlist, drop_introduced, range.polygon = IUCN_Data, native.df = Native_DF))

# GMPD_plots_native_base_02 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Native_Data, plot_type = "iucn")
# names(GMPD_plots_native_base_02) <- sort(Hostlist)
# 
# pdf(file = here::here('Data/GMPD/GMPD_plots_native_base_02.pdf'), width = 10, height = 7)
# GMPD_plots_native_base_02
# dev.off()

## Adjusting subspecies ####

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
nrow(GMPD_Data_res_sub); length(unique(GMPD_Data_res_sub$HostCorrectedName)) # 7084 and 104

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
nrow(GMPD_Data_res_sub); length(unique(GMPD_Data_res_sub$HostCorrectedName))

## IUCN restricted species ####
# new method
Res_Temp <- restrict_deci(GMPD_Data_res_all)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_res_all <- GMPD_Data_res_all[GMPD_Data_res_all$HostCorrectedName %in% Res_Temp$HostCorrectedName,]
nrow(GMPD_Data_res_all); length(unique(GMPD_Data_res_all$HostCorrectedName))

## IUCN cleaned subgroup ####
Res_Temp <- restrict_deci(GMPD_Data_cln, subsp = TRUE)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_cln_sub <- GMPD_Data_cln[GMPD_Data_cln$HostCorrectedName %in% Res_Temp$HostCorrectedName &
                                 GMPD_Data_cln$subgroup %in% Res_Temp$subgroup,]
nrow(GMPD_Data_cln_sub); length(unique(GMPD_Data_cln_sub$HostCorrectedName))

## IUCN cleaned species ####
Res_Temp <- restrict_deci(GMPD_Data_cln)
Res_Temp <- Res_Temp[Res_Temp$enough, ]
GMPD_Data_cln_all <- GMPD_Data_cln[GMPD_Data_cln$HostCorrectedName %in% Res_Temp$HostCorrectedName,]
nrow(GMPD_Data_cln_all); length(unique(GMPD_Data_cln_all$HostCorrectedName))

rm(Res_Temp)

# Misc ########################################################################################################

# narrowing down IUCN_Mammals to a more manageable size
IUCN_Orders <- IUCN_Mammals[IUCN_Mammals$order_ %in% unique(IUCN_Native_Data$order_), ]
rm(IUCN_Mammals)

# Write files #################################################################################################

write.csv(GMPD_Data_res_all, file = here::here("Data/Data back ups/GMPD_Data_res_all_01.csv"), row.names = FALSE)
write.csv(GMPD_Data_res_sub, file = here::here("Data/Data back ups/GMPD_Data_res_sub_01.csv"), row.names = FALSE)

write.csv(GMPD_Data_cln_all, file = here::here("Data/Data back ups/GMPD_Data_cln_all_01.csv"), row.names = FALSE)
write.csv(GMPD_Data_cln_sub, file = here::here("Data/Data back ups/GMPD_Data_cln_sub_01.csv"), row.names = FALSE)

write.csv(GMPD_Location_Data, file = here::here("Data/Data back ups/GMPD_Location_Data_01.csv"), row.names = FALSE)
write.csv(Native_DF, file = here::here("Data/Data back ups/Native_DF_01.csv"), row.names = FALSE)

saveRDS(IUCN_Native_Data, file = here::here("Data/Data back ups/IUCN_Native_Data_01")) # too large to commit 
saveRDS(IUCN_Orders, file = here::here("Data/Data back ups/IUCN_Orders_01")) # too large to commit 

