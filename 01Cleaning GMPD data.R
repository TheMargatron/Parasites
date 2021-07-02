# Reading and cleaning GMPD data by sample and location
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# takes:
## Raw GMPD data from https://doi.org/10.1002/ecy.1799
## Terrestrial mammal IUCN range polygons https://www.iucnredlist.org/resources/spatial-data-download

# makes: 
## GMPD_Data = cleaned GMPD data
## GMPD_Location_Data = associated location data for GMPD samples
## GMPD_GMPD_Trait_Data = associated Trait information summarised from GMPD
## GMPD_Raw_Data 
## Country_Match

############################################## Libraries ######################################################

library(CoordinateCleaner)  # cleaning geographic data
library(geosphere)          # calculating distances
library(here)               #
library(maps)               # iso 3166 country codes and mapnames
library(raster)             # 
library(rgdal)              # read shapefiles
library(stringr)            #
library(tidyverse)          #  beware of conflicts (mainly with raster)

source(here::here("Functions.R"))

GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
nrow(GMPD_Raw_Data) #Beginning with 24323 rows

IUCN_Mammals <- readOGR(here::here("Data/IUCN"), "MAMMALS") #this takes a while

############################################## Basic data cleaning ############################################

# removing: 
## domestic species - may skew results
## marine species - no associated abiotic data or likely poor relationship with abiotic data
## data that is unuseable in future analyses (missing variables)
## unique host-parasite pairs - uninformative in models

GMPD_Data <- GMPD_Raw_Data %>%
  filter(HostCorrectedName != "Ovis aries" &
           HostCorrectedName != "Bos frontalis" &
           HostCorrectedName != "no binomial name") %>%
  filter(!is.na(Prevalence)) %>%
  filter(HostEnvironment != "marine") %>%
  filter(NativeRange != "No" &
           !is.na(NativeRange)) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  filter(!is.na(NumSamples) | !is.na(HostsSampled))
nrow(GMPD_Data) #12328 rows

GMPD_Data <- GMPD_Data %>%
  mutate(HostCorrectedName = case_when(HostCorrectedName == "Alces americanus" ~ "Alces alces",                         # same IUCN polygon
                                       HostCorrectedName == "Felis manul"      ~ "Otocolobus manul",                    # IUCN name differs
                                       HostCorrectedName == "Equus burchellii" ~ "Equus quagga",                        # IUCN name differs
                                       HostCorrectedName == "Taurotragus oryx" ~ "Tragelaphus oryx",                    # IUCN name differs
                                       TRUE                                    ~ HostCorrectedName)) %>%
  select(-ParasiteReportedName, -HostReportedName, -HasBinomialName, -NativeRange, -Intensity, -IntensityMeasure, -SampleNotes) %>%   #not used
  distinct()                                                                                                            #remove duplicated rows
nrow(GMPD_Data) #12263

GMPD_Data <- GMPD_Data %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n()>1) %>% 
  ungroup()
nrow(GMPD_Data) #10378

## adjust prevalence data that has been reported as a percentage

GMPD_Data <- GMPD_Data %>%
  mutate(Prevalence = case_when(Prevalence >1 ~ Prevalence/100,
                                TRUE          ~ Prevalence))

############################################## Sample cleaning ################################################

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
nrow(GMPD_Data) #10223

# host age and sex NA duplications
GMPD_Data <- GMPD_Data %>%
  group_by_at(vars(-HostAge)) %>%
  fill(HostAge, .direction = "updown") %>% 
  
  group_by_at(vars(-HostSex)) %>%
  fill(HostSex, .direction  = "updown") %>% 
  
  group_by_at(vars(-HostAge, -HostSex)) %>%
  fill(c(HostSex, HostAge), .direction = "updown") %>% 
  distinct()
nrow(GMPD_Data) #9784

# HostsSampled vs NumSamples
GMPD_Data <- GMPD_Data %>%
  mutate(HostsSampled = case_when(is.na(HostsSampled) ~ NumSamples,
                                  TRUE                ~ HostsSampled)) %>%
  mutate(NumSamples = case_when(is.na(NumSamples) ~ HostsSampled,
                                TRUE              ~ NumSamples)) %>%
  filter(HostsSampled == NumSamples)
nrow(GMPD_Data) #8860

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
  select(-sample_temp) %>%
  ungroup()
nrow(GMPD_Data) #8819

############################################## Cleaning by location ###########################################

# for matching LocationName with iso3166 country names
## regex adjustments to iso3166 map names deal with incorrect matches
## not future-proof

Country_Match <- iso3166$mapname
Country_Match[81] <- "France(?!s)"
Country_Match[171] <- "Niger(?!i)"
Country_Match[112] <- "India(?!n)"
Country_Match[159] <- "(?<!inner )Mongolia"
Country_Match[122] <- "(?<!New )Jersey"
Country_Match[152] <- "(?<!New )Mexico"
Country_Match[251] <- "(?<![:alpha:])USA"
Country_Match[126] <- "Kenya(?! border)"
Country_Match[86] <- "(?<![:alpha:])UK(?![:alpha:])"
Country_Match[87] <- "(?<!.)Georgia(?!.)"
Country_Match[24] <- "(?<![:alpha:])Saba"
Country_Match[123] <- "(?<![:alpha:])Jordan"
Country_Match[181] <- "(?<![:alpha:])Oman"


Country_Match <- str_c(Country_Match, collapse = "|")
State_Match <- str_c(state.name, collapse = "|")

# All GMPD location descriptions needing matched to a country, 1984 unique descriptions
GMPD_Location_Data <- GMPD_Data %>%
  select(LocationName) %>%
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

# lacking country name in description or misspelled country name (while already having one correct country)
Country_Data_Temp[which(Country_Data_Temp$LocationName == "Masai Mara National Park, Nairobi National Park, Ngorongoro crater,  Serengeti National Park and Namibia"), "countries"] <- "namibia|kenya|tanzania"
Country_Data_Temp[which(Country_Data_Temp$LocationName == "Czech Republic and Slovatkia"), "countries"] <- "czech republic|slovakia"

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

GMPD_Location_Data <- full_join(State_Data_Temp, Country_Data_Temp, by = c("rowname", "LocationName"))

# This bit is very tailored so won't work for any other data
# Number of rows will increase slightly because of location descriptions for which there are multiple countries. 
# They'll be removed again by coordinate cleaner

GMPD_Location_Data <- GMPD_Location_Data %>%
  mutate(countries = case_when(countries == "lebanon"                                      ~ "usa",
                               countries == "trinidad"                                     ~ "argentina",
                               countries == "italy|austria|switzerland"                    ~ "italy",
                               str_detect(LocationName, "Slovak/Hungary")                  ~ "slovakia",
                               countries == "san marino"                                   ~ "usa",
                               str_detect(LocationName, "La Canada-Flintridge")            ~ "usa",
                               countries == "macedonia"                                    ~ "greece",
                               countries == "spain|spain"                                  ~ "spain",
                               countries == "japan|japan"                                  ~ "japan",
                               countries == "namibia|namibia"                              ~ "namibia",
                               countries == "bolivia|bolivia"                              ~ "bolivia",
                               countries == "saint martin|france"                          ~ "france",
                               TRUE                                                        ~ countries)) %>%
  mutate(countries = case_when(countries != ""                                             ~ countries,
                               str_detect(LocationName, "Montana and British Columbia")    ~ "canada|usa",
                               states    != ""                                             ~ "usa",
                               TRUE                                                        ~ "")) %>%
  mutate(countries = case_when(countries != ""                                             ~ countries,
                               str_detect(LocationName, "Serengeti|Ngorongoro|Selous|Temi|Ruaha") ~ "tanzania",
                               str_detect(LocationName, "Kruger National Park|Natal|Skukuza|Queenstown|Eastern Shores|Karroid|KNP|Rooiwal|Rietvlei|Potchefstroom|Benfontein|Pieter|Transvaa|Sabi|Hluhluwe|Kuruman|Mbiyamiti|Ntomeni|West Coast National Park|Transkei") ~ "south africa",
                               str_detect(LocationName, "Noway|Svalbard|Barents")          ~ "norway",
                               str_detect(LocationName, "Bialowie|Polish|Puszcza")         ~ "poland",
                               str_detect(LocationName, "Swiss|Zurich")                    ~ "switzerland",
                               str_detect(LocationName, "Bale|Sidamo|Urso")                ~ "ethiopia",
                               str_detect(LocationName, "Zimbawe|Hippo|Mana pools|Buffalo Range") ~ "zimbabwe",
                               str_detect(LocationName, "Great Britain|Scotland|England|ENGLAND|(?<!of )Wales|WALES|United Kingdom|Oxford|Somerset|Skye|Gloucester|Angus|Clyde") ~ "uk",
                               str_detect(LocationName, "Ontario|Quebec|Saska|Newfoundland|Nova Scotia|Alberta|British|Yukon(?! F)|Brunswick|Prince Ed|Vancouver|Hudson|Vancourver|Baffin|Northwest Territories|Repulse") ~ "canada",
                               str_detect(LocationName, "Kenia|Masai|Nairobi|Jogi|Bungoma") ~ "kenya",
                               str_detect(LocationName, "Cameroun|Ngaoun")                 ~ "cameroon",
                               str_detect(LocationName, "Yellowstone|Channel Islands|Eastern US|(?<![:alpha;])CA(?!.)|Yosemite|Nebrask|Califonia|Orange Free|Purdue|Marion|Susitna|Orgeon|Rocky|Tallahal|ern Arctic|Arctic Nor") ~ "usa",
                               str_detect(LocationName, "Reykjavik")                       ~ "iceland",
                               str_detect(LocationName, "Tokyo|Hokkaido|Hoakkaido|Sobo|Shimane|Akita") ~ "japan",
                               str_detect(LocationName, "Copenhagen")                      ~ "denmark",
                               str_detect(LocationName, "Parana|Leones|Marajo|Jequitinhonha|Grosso|Pantanal") ~ "brazil",
                               str_detect(LocationName, "Brandenburg|Berlin")              ~ "germany",
                               str_detect(LocationName, "Eslovaquia|Zilina|Slovak")        ~ "slovakia",
                               str_detect(LocationName, "Mexican")                         ~ "mexico",
                               str_detect(LocationName, "New Zeland|Flagstaff")            ~ "new zealand",
                               str_detect(LocationName, "Budakeszi")                       ~ "hungary",
                               str_detect(LocationName, "Moravia|Mim")                     ~ "czech republic",
                               str_detect(LocationName, "Karak")                           ~ "jordan",
                               str_detect(LocationName, "Sorbe|Zaragoza|Malaga|Catalonia|Sierras de|Madrid|Sueve|Jaen|Pallars|Aller|Pyrenees") ~ "spain",
                               str_detect(LocationName, "Zaire")                           ~ "democratic republic of the congo",
                               str_detect(LocationName, "Kaa")                             ~ "bolivia",
                               str_detect(LocationName, "Adiopodoume")                     ~ "ivory coast",
                               str_detect(LocationName, "Italian|Sondrio|Brembana|Belviso") ~ "italy",
                               str_detect(LocationName, "Etosha")                          ~ "namibia",
                               str_detect(LocationName, "Rahasthan")                       ~ "india",
                               str_detect(LocationName, "Ankole|Koja|Jie")                 ~ "uganda",
                               str_detect(LocationName, "Meurthe|French|Bauges|Savoy")     ~ "france",
                               str_detect(LocationName, "Umsalala")                        ~ "sudan",
                               str_detect(LocationName, "Dalmatia")                        ~ "croatia",
                               str_detect(LocationName, "Orkendalen")                      ~ "greenland",
                               str_detect(LocationName, "Kerguelen")                       ~ "french southern and antarctic lands",
                               str_detect(LocationName, "alpine areas")                    ~ "italy|switzerland|france",
                               str_detect(LocationName, "Vojvodina")                       ~ "serbia",
                               str_detect(LocationName, "Danubian wetlands")               ~ "romania",
                               str_detect(LocationName, "Bandia|Saboya")                   ~ "senegal",
                               str_detect(LocationName, "Batie")                           ~ "burkina faso",
                               str_detect(LocationName, "Kaisho, Karagwe")                 ~ "tanzania|uganda|zambia",
                               str_detect(LocationName, "Ndoki|Republic of the Congo")     ~ "republic of congo",
                               TRUE                                                        ~ "")) %>%
  separate(countries, into = c("A","B","C"), sep = "\\|") %>%
  pivot_longer(cols = c("A","B","C"), values_to = "mapname", values_drop_na = TRUE) %>%
  select(-name, -states) %>%
  mutate(mapname = case_when(mapname == "uk"      ~ "uk(?!r)",
                             mapname == "norway"  ~ "norway(?!:bouvet|:svalbard|:jan mayen)",
                             mapname == "finland" ~ "finland(?!:aland)",
                             mapname == "china"   ~ "china(?!:hong kong|:macao)",
                             TRUE                 ~ mapname)) # ignore warning

GMPD_Location_Data <- iso3166 %>%
  select(a3, mapname) %>%
  mutate(mapname = tolower(mapname)) %>%
  right_join(GMPD_Location_Data, by = "mapname") %>%
  rename(countrycode = a3)

GMPD_Data <- full_join(GMPD_Location_Data, GMPD_Data, by = "LocationName")

# CoordinateCleaner tests

GMPD_Data <- GMPD_Data %>%
  clean_coordinates(lon = "Longitude",
                    lat = "Latitude",
                    species = "HostCorrectedName",
                    countries = "countrycode",
                    tests = c("capitals","centroids","institutions", "countries"),
                    range_ref = IUCN_Mammals,
                    value = "clean")
nrow(GMPD_Data) #8157

GMPD_Data <- GMPD_Data %>%
  rename(binomial = HostCorrectedName) %>%
  cc_iucn(IUCN_Mammals,
          lon = "Longitude",
          lat = "Latitude",
          species = "binomial") %>%
  rename(HostCorrectedName = binomial)
nrow(GMPD_Data) #7058

GMPD_Data <- GMPD_Data %>%
  group_by(ParasiteCorrectedName, HostCorrectedName) %>%
  filter(n() > 1) %>% 
  ungroup()
nrow(GMPD_Data) #6750

GMPD_Trait_Data <- GMPD_Data %>%
  select(HostCorrectedName, ParasiteCorrectedName, Group, HostOrder, HostFamily, HostEnvironment, ParType, ParPhylum, ParClass) %>%
  unique()

GMPD_Data <- GMPD_Data %>%
  select(HostCorrectedName, ParasiteCorrectedName, 
         Citation, LocationName, Longitude, Latitude, 
         PopulationType, SamplingBasis, Prevalence, 
         HostsSampled, HostSex, HostAge, NumSamples, SamplingType)

############################################## Restricting by proximity #######################################

# Creating nested data frame

Host_Par_Loc_Nest <- GMPD_Data %>%
  select(HostCorrectedName, ParasiteCorrectedName, Longitude, Latitude) %>%
  group_by(HostCorrectedName, ParasiteCorrectedName) %>%
  nest(Location = c(Longitude, Latitude))
nrow(Host_Par_Loc_Nest) # 1165

# Restricting to those that occupy at least two 60/60 res grid squares

Rastr_60 <- raster(resolution = (60/60))

Host_Par_Loc_Nest <- Host_Par_Loc_Nest %>%
  mutate(Across60 = restrict(Location, rastr = Rastr_60)) %>%
  filter(Across60) %>%
  select(-Across60)
nrow(Host_Par_Loc_Nest) # 886

GMPD_Data <- merge(GMPD_Data, Host_Par_Loc_Nest[c(1, 2)], by = c("HostCorrectedName", "ParasiteCorrectedName"), 
                   sort = FALSE, all.x = FALSE)
nrow(GMPD_Data) # 5868

# Saving relevant subset of IUCN data

Hostlist <- unique(GMPD_Data$HostCorrectedName)
IUCN_Data_List <- lapply(Hostlist, function(host) IUCN_Mammals[IUCN_Mammals$binomial == host, ])
names(IUCN_Data_List) <- Hostlist
IUCN_Data_List <- lapply(IUCN_Data_List, function(host) {host@data <- droplevels(host@data); return(host)})
IUCN_Data <- raster::bind(IUCN_Data_List)

############################################### Distance metrics and range traits #############################

Distances_Data <- range_distances(GMPD_Data, IUCN_Data)
GMPD_Data <- merge(GMPD_Data, Distances_Data[["DistanceMetrics"]], all.x = TRUE)
GMPD_Trait_Data <- merge(Distances_Data[["RangeTraits"]], GMPD_Trait_Data, by = "HostCorrectedName", all = TRUE)

############################################### Write files ###################################################

write.csv(GMPD_Data, file = here::here("Data/Data back ups/GMPD_Data.csv"), row.names = FALSE)
write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)
write.csv(GMPD_Location_Data, file = here::here("Data/Data back ups/GMPD_Location_Data.csv"), row.names = FALSE)
saveRDS(Host_Par_Loc_Nest, file = here::here("Data/Data back ups/Host_Par_Loc_Nest"))
saveRDS(IUCN_Data_List, file = here::here("Data/Data back ups/IUCN_Data_List")) # too large to commit 

rm(list = ls(pattern = "_Temp$"))
