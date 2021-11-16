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
library(rnaturalearth)      # river data
library(raster)             # 
library(rgdal)              # read shapefiles
library(rgeos)              # Just for gBuffer
library(stringr)            #
library(tidyverse)          # beware of conflicts (mainly with raster)

source(here::here("Functions.R"))

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
nrow(GMPD_Raw_Data) #Beginning with 24323 rows

IUCN_Mammals <- readOGR(here::here("Data/IUCN"), "MAMMALS") #this takes a while

# Basic data cleaning #########################################################################################

## Removing/adjusting unuseable data ##########################################################################

# removing: 
## domestic species - may skew results
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
  filter(!is.na(Prevalence)) %>%
  filter(HostEnvironment != "marine") %>%
  filter(NativeRange != "No" &
           !is.na(NativeRange)) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  filter(!is.na(NumSamples) | !is.na(HostsSampled))
nrow(GMPD_Data) #12063 rows

GMPD_Data <- GMPD_Data %>%
  mutate(HostCorrectedName = case_when(HostCorrectedName == "Alces americanus" ~ "Alces alces",                         # same IUCN polygon
                                       HostCorrectedName == "Felis manul"      ~ "Otocolobus manul",                    # IUCN name differs
                                       HostCorrectedName == "Equus burchellii" ~ "Equus quagga",                        # IUCN name differs
                                       HostCorrectedName == "Taurotragus oryx" ~ "Tragelaphus oryx",                    # IUCN name differs
                                       HostCorrectedName == "Neotragus moschatus" ~ "Nesotragus moschatus",             # IUCN name differs
                                       TRUE                                    ~ HostCorrectedName)) %>%
  dplyr::select(-ParasiteReportedName, -HostReportedName, -HasBinomialName, -NativeRange, -Intensity, -IntensityMeasure, -SampleNotes) %>%   #not used
  distinct()                                                                                                            #remove duplicated rows
nrow(GMPD_Data) #12005

GMPD_Data <- GMPD_Data %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n()>1) %>% 
  ungroup()
nrow(GMPD_Data) #10164

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
nrow(GMPD_Data) #10010

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
nrow(GMPD_Data) #9592

# HostsSampled vs NumSamples
GMPD_Data <- GMPD_Data %>%
  mutate(HostsSampled = case_when(is.na(HostsSampled) ~ NumSamples,
                                  TRUE                ~ HostsSampled)) %>%
  mutate(NumSamples = case_when(is.na(NumSamples) ~ HostsSampled,
                                TRUE              ~ NumSamples)) %>%
  dplyr::select(-NumSamples)
nrow(GMPD_Data) #9592

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
nrow(GMPD_Data) #9545

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
nrow(GMPD_Data) #8869

# Native and non-native ####
## Removing non-native polygons ####

# reducing to relevant subsets of IUCN data to save space 

IUCN_Data <- IUCN_Mammals[IUCN_Mammals$binomial %in% Hostlist, ]
IUCN_Data <- raster::bind(IUCN_Data, IUCN_Mammals[IUCN_Mammals$binomial == "Cervus canadensis",])
IUCN_Data[IUCN_Data$binomial == "Cervus canadensis", "binomial"] <- "Cervus elaphus"

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
nrow(GMPD_Data_res) #7487

GMPD_Data_res <- GMPD_Data_res %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data_res) #7287

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

## Acinonyx jubatus ####
# rationale: There are multiple subspecies, but GMPD data only occupies the range of A. j. jubatus
# therefore want to remove all polygons representing other subspecies. 
# I've removed polygons based on wiki distribution map
# https://en.wikipedia.org/wiki/Cheetah#/media/File:Acinonyx_jubatus_subspecies_range_IUCN_2015.png

sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == "Acinonyx jubatus", ]

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
PolyRm_try <- gUnion(PolyRm, sp.range.polygon[sp.range.polygon$poly == 13, ])

sp.range.polygon_try <- sp.range.polygon - PolyRm_try

sf::sf_use_s2(FALSE)
tmap_mode("plot")
tm_shape(sp.range.polygon) + tm_polygons("legend") #+ tm_shape(PolyRm_try) + tm_polygons()
tm_shape(sp.range.polygon_try) + tm_polygons("legend")

## Aepyceros melampus ####
# rationale: Two subspecies, only the common impala (subsp. melampus) is well represented in GMPD based on wiki and iucn maps
# https://en.wikipedia.org/wiki/Impala#/media/File:Aepyceros_melampus.svg
# IUCN: "In Namibia, the Black-faced Impala is naturally confined to the Kaokoland in the north-west, and neighbouring south-western Angola"

sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == "Aepyceros melampus", ]
coords = matrix(c(17, -20,
                  10, -20,
                  10, -13,
                  17, -13,
                  17, -20), 
                ncol = 2, byrow = TRUE)

PolyRm <- Polygon(coords)
PolyRm <- SpatialPolygons(list(Polygons(list(PolyRm), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

tm_shape(PolyRm, bbox = sp.range.polygon@bbox) + tm_polygons() + tm_shape(sp.range.polygon) + tm_polygons("legend") 

## Alcelaphus buselaphus ####
# rationale: 8 subspecies, GMPD data appears to represent major and cokii
# https://en.wikipedia.org/wiki/Hartebeest#/media/File:Alcelaphus_recent.png
# There is only one sample location for each subspecies so I'll have to drop this one :(
GMPD_Data[which(GMPD_Data$HostCorrectedName == "Alcelaphus buselaphus"),c("Latitude", "Longitude")]

## Alces alces ####
# rationale: There are multiple subspecies, GMPD represents shirasi, gigas, andersoni, and americana in North America
# and alces in Europe/Western Russia, but not buturlini, cameloides, or pfizenmayeri (east of Yenisei river)
# Splitting European polygon around Yenisei should do it. 
# Yenisei data source: https://doi.org/10.1016/j.dib.2018.09.016

sp.range.polygon <- IUCN_Native_Data[IUCN_Native_Data$binomial == "Alces alces", ]

River_Data50 <- ne_load(scale = 50,
                             type = "rivers_lake_centerlines",
                             category = "physical",
                             destdir = here::here("Data/Extras/ne_rivers"))

YeniAnga <- River_Data50[(River_Data50$name == "Yenisey" | River_Data50$name == "Angara"),]

YA_Temp <- disaggregate(YeniAnga)
YA_Temp$ID <- LETTERS[1:12]
tm_shape(YA_Temp) + tm_lines("ID", lwd = 2)
tm_shape(YA_Temp[str_detect(YA_Temp$ID, "A|G|I|K", negate = TRUE),]) + tm_lines("ID", lwd = 2)
YA_Temp <- YA_Temp[-grep("A|G|I|K", YA_Temp$ID),]
YA_coords <- unlist(coordinates(YA_Temp), recursive = FALSE)
YA_coords[[2]] <- YA_coords[[2]][nrow(YA_coords[[2]]):1,]
YA_coords[[8]] <- YA_coords[[8]][nrow(YA_coords[[8]]):1,]
YA_coords <- YA_coords[c(3,1,4,2,5,6,7,8)]
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
tm_shape(sp.range.polygon_try) + tm_polygons()

## Next species ####
# Antidorcas marsupilis







## discard pile ####

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
tmap_mode("view")
tm_shape(sp.range.polygon_try) + tm_polygons("binomial")