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

library(CoordinateCleaner)  # cleaning geographic data
library(geosphere)          # calculating distances
library(here)               #
library(maps)               # iso 3166 country codes and mapnames
library(maptools)
library(rnaturalearth)      # river data
library(raster)             # 
library(rgdal)              # read shapefiles
library(rgeos)              # Just for gBuffer
library(stringr)            #
library(tidyverse)          # beware of conflicts (mainly with raster)
library(tmap)

source(here::here("Functions.R"))

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
nrow(GMPD_Raw_Data); length(unique(GMPD_Raw_Data$HostCorrectedName)) #Beginning with 24323 rows and 462 hosts

IUCN_Mammals <- readOGR(here::here("Data/IUCN"), "MAMMALS") #this takes a while

River_Data50 <- ne_load(scale = 50,
                        type = "rivers_lake_centerlines",
                        category = "physical",
                        destdir = here::here("Data/Extras/ne_rivers"))

sf::sf_use_s2(FALSE) #Not sure about keeping this here. May move. For "invalid spherical geometry" errors
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
GMPD_base_plots_00 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, plot_type = "base")
names(GMPD_base_plots_00) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_base_plots_00.pdf'), width = 10, height = 7)
GMPD_base_plots_00
dev.off()

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
GMPD_plots_01 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Data, plot_type = "iucn")
names(GMPD_plots_01) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_plots_01.pdf'), width = 10, height = 7)
GMPD_plots_01
dev.off()

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

GMPD_plots_native_base_02 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data, range.polygon = IUCN_Native_Data, plot_type = "iucn")
names(GMPD_plots_native_base_02) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_plots_native_base_02.pdf'), width = 10, height = 7)
GMPD_plots_native_base_02
dev.off()

## Removing subspecies ####

## Restricting by native IUCN polygon ####

GMPD_Data_res <- GMPD_Data %>%                        
  rename(binomial = HostCorrectedName) %>%
  cc_iucn(IUCN_Native_Data,
          lon = "Longitude",
          lat = "Latitude",
          species = "binomial") %>%
  rename(HostCorrectedName = binomial)
nrow(GMPD_Data_res); length(unique(GMPD_Data_res$HostCorrectedName)) # 7458 and 126

GMPD_Data_res <- GMPD_Data_res %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data_res); length(unique(GMPD_Data_res$HostCorrectedName)) # 7287 and 120

GMPD_plots_native_restricted_03 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_res, range.polygon = IUCN_Native_Data, plot_type = "iucn")
names(GMPD_plots_native_restricted_03) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_plots_native_restricted_03.pdf'), width = 10, height = 7)
GMPD_plots_native_restricted_03
dev.off()

## Cleaning by native IUCN polygon ####
Native_Clean <- Native_DF %>%
  filter(!keep) %>%
  filter(HostCorrectedName %in% unique(GMPD_Data$HostCorrectedName)) %>%
  mutate(buffer = case_when(str_detect(HostCorrectedName, "Capra ibex|Martes melampus") ~ 0,
                            HostCorrectedName == "Vulpes vulpes"                        ~ 2,
                            TRUE                                                        ~ 1))

GMPD_clean_plots <- lapply(sort(unique(Native_Clean$HostCorrectedName)), iucn_test)
names(GMPD_clean_plots) <- sort(unique(Native_Clean$HostCorrectedName))

GMPD_Data_cln <- iucn_cleaning(dat = GMPD_Data, native.df = Native_Clean)
GMPD_Data_cln <- GMPD_Data_cln %>%
  filter(is.na(out) | !out) %>%
  dplyr::select(-out)

GMPD_plots_native_clean_03 <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_cln, range.polygon = IUCN_Native_Data, plot_type = "iucn")
names(GMPD_plots_native_clean_03) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_plots_native_clean_03.pdf'), width = 10, height = 7)
GMPD_plots_native_clean_03
dev.off()

# Restricting each by proximity ####################################################################################
## IUCN restricted ####
# Creating nested data frame

Host_Par_Loc_Nest_res <- GMPD_Data_res %>%
  dplyr::select(HostCorrectedName, ParasiteCorrectedName, Longitude, Latitude) %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  nest(Location = c(Longitude, Latitude))
nrow(Host_Par_Loc_Nest_res) # 1221

# Restricting to those that occupy at least two 60/60 res grid squares

Rastr_60 <- raster(resolution = (60/60))

Host_Par_Loc_Nest_res <- Host_Par_Loc_Nest_res %>%
  mutate(Across60 = restrict(Location, rastr = Rastr_60)) %>%
  filter(Across60) %>%
  dplyr::select(-Across60)
nrow(Host_Par_Loc_Nest_res) # 875

GMPD_Data_res <- merge(GMPD_Data_res, Host_Par_Loc_Nest_res[c(1, 2)], by = c("HostCorrectedName", "ParasiteCorrectedName"), 
                   sort = FALSE, all.x = FALSE)
nrow(GMPD_Data_res) # 6135

# Plotting again

GMPD_plots_clean_04_restricted <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_res, range.polygon = IUCN_Native_Data, plot_type = "iucn")
names(GMPD_plots_clean_04_restricted) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_plots_clean_04_restricted.pdf'), width = 10, height = 7)
GMPD_plots_clean_04_restricted
dev.off()

## IUCN cleaned ####
# Creating nested data frame

Host_Par_Loc_Nest_cln <- GMPD_Data_cln %>%
  dplyr::select(HostCorrectedName, ParasiteCorrectedName, Longitude, Latitude) %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  nest(Location = c(Longitude, Latitude))
nrow(Host_Par_Loc_Nest_cln) # 1582

# Restricting to those that occupy at least two 60/60 res grid squares

Host_Par_Loc_Nest_cln <- Host_Par_Loc_Nest_cln %>%
  mutate(Across60 = restrict(Location, rastr = Rastr_60)) %>%
  filter(Across60) %>%
  dplyr::select(-Across60)
nrow(Host_Par_Loc_Nest_cln) # 1079

GMPD_Data_cln <- merge(GMPD_Data_cln, Host_Par_Loc_Nest_cln[c(1, 2)], by = c("HostCorrectedName", "ParasiteCorrectedName"), 
                   sort = FALSE, all.x = FALSE)
nrow(GMPD_Data_cln) # 7337

# Plotting again

GMPD_plots_clean_04_clean <- lapply(sort(Hostlist), FUN = gmpd_plotter, dat = GMPD_Data_cln, range.polygon = IUCN_Native_Data, plot_type = "iucn")
names(GMPD_plots_clean_04_clean) <- sort(Hostlist)

pdf(file = here::here('Data/GMPD/GMPD_plots_clean_04_clean.pdf'), width = 10, height = 7)
GMPD_plots_clean_04_clean
dev.off()

# Misc ########################################################################################################

# narrowing down IUCN_Mammals to a more manageable size
IUCN_Orders <- IUCN_Mammals[IUCN_Mammals$order_ %in% unique(IUCN_Native_Data$order_), ]
rm(IUCN_Mammals)

# Write files #################################################################################################

write.csv(GMPD_Data_res, file = here::here("Data/Data back ups/GMPD_Data_res_01.csv"), row.names = FALSE)
write.csv(GMPD_Data_cln, file = here::here("Data/Data back ups/GMPD_Data_cln_01.csv"), row.names = FALSE)
write.csv(GMPD_Location_Data, file = here::here("Data/Data back ups/GMPD_Location_Data_01.csv"), row.names = FALSE)
write.csv(Native_DF, file = here::here("Data/Data back ups/Native_DF_01.csv"), row.names = FALSE)

saveRDS(IUCN_Native_Data, file = here::here("Data/Data back ups/IUCN_Native_Data_01")) # too large to commit 
saveRDS(IUCN_Orders, file = here::here("Data/Data back ups/IUCN_Orders_01")) # too large to commit 

# Figuring things out ####

# going by subgroup rather than subspecies because of geographic groupings
IUCN_Native_Data@data$subgroup <- IUCN_Native_Data@data$subspecies

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(str_detect(subgroup, " ssp. ") ~ str_split(subgroup, "ssp. ", simplify = TRUE)[,2],
                                TRUE ~ subgroup))


## Acinonyx jubatus ####
# rationale: There are multiple subspecies that are quite geographically separated, 
# but GMPD data only occupies the range of A. j. jubatus
# therefore want to remove all polygons representing other subspecies. 
# I've removed polygons based on wiki distribution map
# https://en.wikipedia.org/wiki/Cheetah#/media/File:Acinonyx_jubatus_subspecies_range_IUCN_2015.png

curr.species <- "Acinonyx jubatus"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
tm_shape(sp.range.polygon) + tm_polygons("subgroup")

# upper
coords = matrix(c(33.0, 0,
                  -10, 0,
                  -10, 40,
                  60.0, 40,
                  60.0, 3.7,
                  45.0, 3.7,
                  35.5, 5,
                  33.0, 0), 
                ncol = 2, byrow = TRUE)

PolyRm <- Polygon(coords)
PolyRm <- SpatialPolygons(list(Polygons(list(PolyRm), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
PolyRm <- gUnion(PolyRm, sp.range.polygon[13, ])

sp.range.polygon_jubatus <- sp.range.polygon - PolyRm
sp.range.polygon_jubatus@data$subgroup <- "jubatus"
sp.range.polygon_otherssp <- raster::intersect(sp.range.polygon, PolyRm)
sp.range.polygon_try <- raster::bind(sp.range.polygon_jubatus, sp.range.polygon_otherssp)

tm_shape(sp.range.polygon_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sp.range.polygon_try)

rm(list = ls(pattern = "^sp.range.polygon"))

## Aepyceros melampus ####
# rationale: Two subspecies, only the common impala (subsp. melampus) is well represented in GMPD based on wiki and iucn maps
# Also black-faced impala seems geographically distinct 
# https://en.wikipedia.org/wiki/Impala#/media/File:Aepyceros_melampus.svg
# IUCN: "In Namibia, the Black-faced Impala is naturally confined to the Kaokoland in the north-west, and neighbouring south-western Angola"

curr.species <- "Aepyceros melampus"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
tm_shape(sp.range.polygon) + tm_polygons("subgroup")  

## Alcelaphus buselaphus ####
# rationale: 8 subgroup, GMPD data appears to represent major and cokii, and they are geographically distinct
# https://en.wikipedia.org/wiki/Hartebeest#/media/File:Alcelaphus_recent.png
# There is only one sample location for each subspecies so I'll have to drop this one :(

curr.species <- "Alcelaphus buselaphus"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]

## Alces alces ####
# rationale: There are multiple subspecies, GMPD represents shirasi, gigas, andersoni, and americana in North America
# and alces in Europe/Western Russia, but not buturlini, cameloides, or pfizenmayeri (east of Yenisei river)
# Splitting European polygon around Yenisei should do it. 
# Yenisei data source: https://doi.org/10.1016/j.dib.2018.09.016

curr.species <- "Alces alces"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]

# Used combination of Yenisei and Angara because it matched where holes were in polygon, split generously, and crossed the full polygon width
YeniAnga <- River_Data50[(River_Data50$name == "Yenisey" | River_Data50$name == "Angara"),]

YA_Temp <- disaggregate(YeniAnga)
YA_Temp$ID <- LETTERS[1:nrow(YA_Temp)]
tm_shape(YA_Temp) + tm_lines("ID", lwd = 2)
drop_ID <- "A|G|I|K"
tm_shape(YA_Temp[str_detect(YA_Temp$ID, drop_ID, negate = TRUE),]) + tm_lines("ID", lwd = 2)

YA_Temp <- YA_Temp[-grep(drop_ID, YA_Temp$ID),] # to avoid self intersections
YA_coords <- unlist(coordinates(YA_Temp), recursive = FALSE)
names(YA_coords) <- YA_Temp$ID
YA_coords[["C"]] <- YA_coords[["C"]][nrow(YA_coords[["C"]]):1,]
YA_coords[["L"]] <- YA_coords[["L"]][nrow(YA_coords[["L"]]):1,]
YA_coords <- YA_coords[c("D","B","E","C","F","H","J","L")]
YA_coords <- do.call(rbind, YA_coords)

YA_coords <- rbind(YA_coords,
                   matrix(c(YA_coords[nrow(YA_coords), 1], (sp.range.polygon@bbox[2,2] + 1),
                            (sp.range.polygon@bbox[1,2] + 1), (sp.range.polygon@bbox[2,2] + 1),
                            (sp.range.polygon@bbox[1,2] + 1), (sp.range.polygon@bbox[2,1] - 1),
                            YA_coords[1,1], (sp.range.polygon@bbox[2,1] - 1),
                            YA_coords[1,1], YA_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

PolyRm <- Polygon(YA_coords)
PolyRm <- SpatialPolygons(list(Polygons(list(PolyRm), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

tm_shape(PolyRm) + tm_polygons()

sp.range.polygon_try <- sp.range.polygon - PolyRm

alces <- Polygon(matrix(c(0, 75,
                          110, 75,
                          110, 30,
                          0, 30,
                          0, 75), 
                        ncol = 2, byrow = TRUE))
alces <- SpatialPolygons(list(Polygons(list(alces), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sp.range.polygon_alces <- raster::intersect(sp.range.polygon_try, alces)
sp.range.polygon_alces@data$subgroup <- "alces"
sp.range.polygon_americanus <- sp.range.polygon_try - alces
sp.range.polygon_americanus@data$subgroup <- "americana"

sp.range.polygon_try <- raster::bind(sp.range.polygon_alces, sp.range.polygon_americanus)
tm_shape(sp.range.polygon_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sp.range.polygon_try)

rm(list = ls(pattern = "^sp.range.polygon"))

## Antidorcas marsupilis ####
# rationale: There are three recognised subspecies (wiki) but their individual distributions are not well defined
# They are also not described by the IUCN
# I only seem to have A. m. marsupialis in GMPD data whose range lies south of Orange river (wiki) 
# In the Eastern part of this range it appears to be restricted by the Vaal rather than the Orange as wiki describes it as extending up to Kimberley which is North of Orange

curr.species <- "Antidorcas marsupialis"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
OrangeVaal <- River_Data50[(River_Data50$name == "Orange"| River_Data50$name == "Vaal"),]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

tm_shape(sp.range.polygon) + tm_polygons() + tm_shape(OrangeVaal) + tm_lines() +
  tm_shape(sp.gmpd.points) + tm_dots()

OV_Temp <- disaggregate(OrangeVaal)
OV_Temp$ID <- LETTERS[1:nrow(OV_Temp)]
tm_shape(OV_Temp) + tm_lines("ID", lwd = 2)
drop_ID <- "G|A|E|B"
tm_shape(OV_Temp[str_detect(OV_Temp$ID, drop_ID, negate = TRUE),]) + tm_lines("ID", lwd = 2)
OV_Temp <- OV_Temp[-grep(drop_ID, OV_Temp$ID),] # to avoid self intersections

OV_edit <- coordinates(OV_Temp[OV_Temp$ID == "F",])[[1]][[1]]
split_point <- nearestPointOnLine(OV_edit, 
                                  tail(coordinates(OV_Temp[OV_Temp$ID == "H",])[[1]][[1]], 1L))

split_row <- which(apply(OV_edit, 1, function(x) all(x == split_point)))
OV_edit <- OV_edit[split_row:nrow(OV_edit),]
tm_shape(OV_Temp) +tm_lines() + 
  tm_shape(SpatialLines(list(Lines(Line(OV_edit), ID = "6"))))  + tm_lines( col = "red")

OV_coords <- unlist(coordinates(OV_Temp), recursive = FALSE)
names(OV_coords) <- OV_Temp$ID
OV_coords[["F"]] <- OV_edit

# Turning them round
OV_coords <- OV_coords[c("J","D","I","C","H","F")]
OV_coords <- do.call(rbind, OV_coords)

OV_coords <- rbind(OV_coords,
                   matrix(c((sp.range.polygon@bbox[1,1] - 1), OV_coords[nrow(OV_coords), 2],
                            (sp.range.polygon@bbox[1,1] - 1), (sp.range.polygon@bbox[2,2] + 1),
                            (sp.range.polygon@bbox[1,2] + 1), (sp.range.polygon@bbox[2,2] + 1),
                            (sp.range.polygon@bbox[1,2] + 1), OV_coords[1,2], 
                            OV_coords[1,1], OV_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

PolyRm <- Polygon(OV_coords)
PolyRm <- SpatialPolygons(list(Polygons(list(PolyRm), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

tm_shape(PolyRm) + tm_polygons()

sp.range.polygon_marsupialis <- sp.range.polygon - PolyRm
sp.range.polygon_marsupialis@data$subgroup <- "marsupialis"
sp.range.polygon_otherssp <- raster::intersect(sp.range.polygon, PolyRm)
sp.range.polygon_try <-  raster::bind(sp.range.polygon_marsupialis, sp.range.polygon_otherssp)
tm_shape(sp.range.polygon_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sp.range.polygon_try)

rm(list = ls(pattern = "^sp.range.polygon"))

## Antilocapra americana ####
# rationale: three subspecies, but not all represented by GMPD
# "Mitochondrial DNA analyses since the early 1990s support the idea of clines within a wide-ranging species rather than separate subspecies (O’Gara and Yoakum 2004)."

curr.species <- "Antilocapra americana"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

tm_shape(sp.range.polygon) + tm_polygons("legend") +  tm_shape(sp.gmpd.points) + tm_dots()

sp.range.polygon_try <- raster::disaggregate(sp.range.polygon)
tm_shape(sp.range.polygon_try[17,]) + tm_polygons()
sp.range.polygon_try@data$subgroup <- "americana x sonoriensis"
sp.range.polygon_try@data[17,"subgroup"] <- "peninsularis"
tm_shape(sp.range.polygon_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sp.range.polygon_try)

rm(list = ls(pattern = "^sp.range.polygon"))

## Axis axis ####
# rationale: No description of subspecies in either IUCN or main wiki page, 
# but there's a page for sri lankan subspecies: https://en.wikipedia.org/wiki/Sri_Lankan_axis_deer
# and it's geographically isolated

curr.species <- "Axis axis"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
tm_shape(sp.range.polygon) + tm_polygons("island")

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & island == "Sri Lanka" ~ "ceylonensis",
                              binomial == curr.species ~ "axis",
                              TRUE ~ subgroup))

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons("subgroup")

## Bison bison ####
# rationale: "There are two recognized subspecies in North America: Bison bison bison and B. b. athabascae."

curr.species <- "Bison bison"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sp.range.polygon) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

# not sure how to split these for now because gbif data does not have clear delineation
# and subspecies data seems to be the other way round in gbif?
# probs just won't split them as they are not too geographically distinct anyway

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species ~ "bison x athabascae",
                              TRUE ~ subgroup))

## Bison bonasus ####
# rationale: "Three subspecies of European bison existed in the recent past, but only one, 
# the nominate subspecies (B. b. bonasus), survives today" (wiki)
# two subspecies are widely recognized as the Lowland Bison (Bison bonasus bonasus) and 
# the Caucasian Bison (Bison bonasus caucasicus) (Kowalczyk and Plumb 2020)" (IUCN)

# "European bison herds, scattered across Central and Eastern Europe, represent two genetic lines"
# the lowland line (Poland, Belarus, and Lithuania) and the lowland-Caucasian line (southern Poland, Russia, Ukraine and Slovakia).
# https://animaldiversity.org/accounts/Bison_bonasus/

# animaldiversity doesn't describe them as subspecies, only genetic lines
# The distribution of the lineages is artificial because it's determined by reintroduction programmes, rather than natural dispersal
# That means that the different habitat patches they occupy will not relate to their lineage, just to their overall species distribution

curr.species <- "Bison bonasus"
IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species ~ "bonasus x caucasicus",
                              TRUE ~ subgroup))

## Blastocerus dichotomus ####
# rationale: No taxonomic notes on IUCN, no description of subspecies on wiki or animaldiversity.org

## Canis adustus ####
# rationale: There are seven recognized subspecies of the side-striped jackal:[2]

# L. a. adusta (West Africa to most of Angola) – Sundevall's side-striped jackal
## L. a. bweha (East Africa; Kisumu, Kenya) – Elgon side-striped jackal[15]
# L. a. centralis (Central Africa; Cameroon, near the Uham River)
# L. a. grayi (North Africa; Morocco and Tunisia)
# L. a. kaffensis (Kaffa, southwestern Ethiopia) – Kaffa side-striped jackal
## L. a. lateralis (East Africa; Kenya, Uasin Gishu Plateau, south of Gabon)
## L. a. notatus (East Africa; Kenya, Loita Plains, Rift Valley Province) – Loita side-striped jackal[15]

# Leaving this one be because I can't accurately split the polygons and I'm not 100% on the subspecies
# Kenya points could be bweha, lateralis, or notatus
# unclear what Zimbabwe point is as none are described as being in Southern African

# I think the species where I can't distinguish between subspecies will still be informative for latitudinal gradient,
# but less informative for within-range analysis.
# Might be interesting to compare result from uncleaned polygons with cleaned. Split the streams again after adding this step?

curr.species <- "Canis adustus"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sp.range.polygon) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

## Canis aureus ####
"Recent studies based on mtDNA and morphology have shown that 'Golden Jackals' in Africa are larger in size than those from 
Eurasia and are actually more closely related to the Grey Wolf Canis lupus. African animals hence represent a previously 
overlooked distinct species, the African Wolf, Canis lupaster (see Rueness et al. 2011, Gaubert et al. 2012, Koepfli et al. 
2015, Viranta et al. 2017). However, the putative presence of Golden Jackal in the Sinai Peninsula of Egypt remains unclear 
(see Gaubert et al. 2012, Viranta et al. 2017)."

# rationale: wiki s 7 subspecies
# aureus: Middle East, Iran, Turkmenistan, Afghanistan, Pakistan and Western India
# cruesemanni: Thailand
# ecsedensis: Pannonian Basin, Central Europe
# indicus: India, Nepal, Bangladesh, Bhutan
# moreotica: Southeastern Europe, Moldova, Asia Minor and the Caucasus
# naria: Coastal South West India, Sri Lanka
# syriacus: Israel, Syria,[38] Lebanon,[62] and Jordan

# Have aureus in Iran and probably moreotica in Greece
# subdividing according to : https://doi.org/10.1093/mspecies/sey002
# Following IUCN range description I'm not including the algirensis subspecies
curr.species <- "Canis aureus"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sp.range.polygon) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

# moreotica
sp.range.polygon_try <- raster::disaggregate(sp.range.polygon)
sp.range.polygon_moreotica <- sp.range.polygon_try[2,]

# removing ecsedensis (Pannonian basin/Hungary)
sp.range.polygon_moreotica <- sp.range.polygon_moreotica - countries50[countries50$name == "Hungary",]
sp.range.polygon_moreotica <- raster::disaggregate(sp.range.polygon_moreotica)
sp.range.polygon_moreotica <- sp.range.polygon_moreotica[1,]

# adding Turkey
turkey_range_temp <- countries50[countries50$name == "Turkey",]
turkey_range_temp <- raster::intersect(sp.range.polygon_try[1,], turkey_range_temp)
turkey_range_temp@data <- turkey_range_temp@data[, names(sp.range.polygon_try)]

sp.range.polygon_moreotica <- raster::bind(sp.range.polygon_moreotica, 
                                            turkey_range_temp)

# adding range around azov sea, split at Georgia boundary
azov_range_temp <- Polygon(matrix(c(39.75, 43,
                           34, 43,
                           34, 49,
                           43.75, 49,
                           39.75, 43), 
                        ncol = 2, byrow = TRUE))
azov_range_temp <- SpatialPolygons(list(Polygons(list(azov_range_temp), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
sp.range.polygon_moreotica <- raster::bind(sp.range.polygon_moreotica, 
                                            raster::intersect(azov_range_temp, sp.range.polygon_try[1,]))
sp.range.polygon_moreotica@data$subgroup <- "moreotica"
tm_shape(sp.range.polygon_moreotica) + tm_polygons("subgroup")

# aureus
PolyRm <- raster::bind(turkey_range_temp, 
                       countries50[countries50$name == "Jordan",],
                       countries50[countries50$name == "Syria",],
                       countries50[countries50$name == "India",],
                       azov_range_temp)
sp.range.polygon_aureus <- sp.range.polygon_try[1,] - PolyRm
sp.range.polygon_aureus <- raster::disaggregate(sp.range.polygon_aureus)
sp.range.polygon_aureus <- sp.range.polygon_aureus[1,]
sp.range.polygon_aureus@data$subgroup <- "aureus"
tm_shape(sp.range.polygon_aureus) + tm_polygons("subgroup")

# adding it all up
sp.range.polygon_try <- raster::bind(sp.range.polygon_moreotica, sp.range.polygon_aureus)
sp.range.polygon_otherssp <- sp.range.polygon - sp.range.polygon_try
sp.range.polygon_try <- raster::bind(sp.range.polygon_try, sp.range.polygon_otherssp)
tm_shape(sp.range.polygon_try) + tm_polygons("subgroup")

sp.range.polygon_try <- raster::disaggregate(sp.range.polygon_try)
sp.range.polygon_try@data$SHAPE_Area <- raster::area(sp.range.polygon_try)
sp.range.polygon_try <- sp.range.polygon_try[sp.range.polygon_try$SHAPE_Area > 100000000,]
sp.range.polygon_try <- raster::aggregate(sp.range.polygon_try, by = names(sp.range.polygon_try))

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sp.range.polygon_try)

rm(list = ls(pattern = "^sp.range.polygon"))
rm(list = ls(pattern = "_temp$"))

## Canis latrans ####
# rationale: although subspecies ranges are well described in wiki, there's overlap between them 
# and I don't know exactly where to draw lines. 
# Also range is pretty contiguous (apart from tiburon island). The only subspecies I'm lacking are central american ones

### Canis lupus ####
# for Canis lupus rufus / Canis rufus:
"See Chambers et al. (2012) for a brief review of recent literature concerning the status of this species, 
which they considered a full species, as does this assessment."
# rationale: absolutely loads of subspecies: https://en.wikipedia.org/wiki/Subspecies_of_Canis_lupus
# based on: https://upload.wikimedia.org/wikipedia/commons/1/16/Present_distribution_of_gray_wolf_%28canis_lupus%29_subspecies.png
# Not got: baileyi, arctos, albus, arabs, nubilis
# Got: occidentalis, lycaon, signatus, lupus, italicus, pallipes

# use IUCN taxonomic notes


# removing arctos
sp.range.polygon_try <- sp.range.polygon[-grep("Greenland|Ellesmere|Banks|Melville", sp.range.polygon$island),]

## Canis mesomelas ####
# rationale: two geographically distinct subspecies
curr.species <- "Canis mesomelas"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sp.range.polygon) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

sp.over <- over(sp.gmpd.points, gBuffer(sp.range.polygon, byid = TRUE))
sp.gmpd$subspecies <- sp.over$subspecies
sp.gmpd$HostCorrectedName <- paste(sp.gmpd$HostCorrectedName, sp.over$subspecies, sep = " ")

## Canis simensis ####
# rationale: Two geographically distinct subspecies and I only seem to have citernii
# https://en.wikipedia.org/wiki/Ethiopian_wolf#/media/File:Canis_simensis_subspecies_range.png
curr.species <- "Canis simensis"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sp.range.polygon) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

coords = matrix(c(37, 9,
                  40, 9,
                  40, 14,
                  37, 14,
                  37, 9), 
                ncol = 2, byrow = TRUE)

PolyRm <- Polygon(coords)
PolyRm <- SpatialPolygons(list(Polygons(list(PolyRm), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sp.range.polygon_try <- sp.range.polygon - PolyRm
tm_shape(PolyRm) + tm_polygons() + tm_shape(sp.range.polygon_try) + tm_polygons()
sp.range.polygon_try@data$subgroup <- "citernii"

sp.range.polygon_otherssp <- raster:intersect(sp.range.polygon, PolyRm)

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sp.range.polygon_try)

rm(list = ls(pattern = "^sp.range.polygon"))

## Capra ibex ####
# rationale: no reported subspecies

### Capra pyrenaica ####
# rationale: four subspecies, two extinct. I probably only have hispanica

## Capreolus capreolus ####


## discard pile ####
curr.species <- "Ovibos moschatus"
sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.range.polygon <- IUCN_Mammals[IUCN_Mammals$binomial == curr.species, ]

sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostReportedName == "Viverra civetta", c("Longitude", "Latitude")])

tm_shape(World, bbox = sp.range.polygon@bbox) + tm_polygons() + tm_shape(sp.range.polygon) + tm_polygons("legend") + tm_shape(sp.gmpd.points) + tm_dots()



sp.bbox <- sp.range.polygon@bbox
sp.range.polygon$id <- rownames(sp.range.polygon@data)
sp.range.polygon$poly <- 1:36
sp.range.polygon_fort <- merge(fortify(sp.range.polygon), sp.range.polygon@data, bi = "id")

# upper
coords = matrix(c(33.0, 0,
                  -10, 0,
                  -10, 40,
                  60.0, 40,
                  60.0, 3.7,
                  45.0, 3.7,
                  35.5, 5,
                  33.0, 0), 
                ncol = 2, byrow = TRUE)

# lower?
coords = matrix(c(33.0, 0,
                  35.5, 5,
                  45.0, 3.7,
                  60.0, 3.7,
                  60.0, -30,
                  -10, -30,
                  -10, 0,
                  33.0, 0), 
                ncol = 2, byrow = TRUE)

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  coord_sf(xlim = c(30,50), ylim = c(-5,15), expand = TRUE) +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = sp.range.polygon_fort, 
                aes(x = long, y = lat, group = group, colour = poly, fill = poly),
                size = 0.05) +
  geom_polypath(data = as.data.frame(coords), aes(x = V1, y = V2), fill = NA, colour = "black")

PolyRm <- Polygon(coords)
PolyRm <- SpatialPolygons(list(Polygons(list(PolyRm), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
PolyRm_try <- gUnion(PolyRm, sp.range.polygon[sp.range.polygon$poly == 13, ])

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  coord_sf(xlim = sp.bbox[1,], ylim = sp.bbox[2,], expand = TRUE) +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = sp.range.polygon_fort, 
                aes(x = long, y = lat, group = group, colour = poly, fill = poly),
                size = 0.05) +
  geom_polypath(data = fortify(PolyRm_try), aes(x = long, y = lat), fill = NA, colour = "black")

sp.range.polygon.test <- st_as_sf(sp.range.polygon)
PolyRm.test <- st_as_sf(PolyRm_try)
sp.range.polygon_try <- st_difference(sp.range.polygon.test, PolyRm.test)

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  coord_sf(xlim = sp.bbox[1,], ylim = sp.bbox[2,], expand = TRUE) +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = sp.range.polygon_fort, 
                aes(x = long, y = lat, group = group, colour = poly, fill = poly),
                size = 0.05) +
  geom_polypath(data = (sp.range.polygon_try), aes(x = long, y = lat), fill = "black", colour = "black")



##
#library(mapview)
#geojson.io
sf::sf_use_s2(FALSE)
library(tmap)
tmap_mode("plot")
tm_shape(sp.range.polygon_try) + tm_polygons("binomial")