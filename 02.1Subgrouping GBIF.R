# Subgroup method
plot_temp <- function(host, infra = TRUE){
  sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == host,]
  sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == host,]
  if(infra){
    tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gbif) + tm_dots("subgroup")
  } else{
    tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gbif) + tm_dots()
  }
}

# small sample sizes ####
# Some species have very few samples so need revisiting
# https://doi.org/10.1111/ecog.01509
few.gbif <- GBIF_Spatial@data %>% group_by(species) %>% summarise(n = n()) %>% filter(n < 100)
few.gbif <- left_join(few.gbif, summarise(group_by(GBIF_Raw_Data, species), rawn = n()), by = "species")

few.raw <- GBIF_Raw_Data[GBIF_Raw_Data$species %in% few.gbif$species,] %>%
  filter(countryCode != "" & countryCode != "XK" & countryCode != "ZZ") %>%
  mutate(countryCode = countrycode(countryCode, origin = "iso2c", destination = "iso3c"))

few.raw <- clean_coordinates(few.raw,
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

few.gbif <- left_join(few.gbif, summarise(group_by(few.raw, species), ccn = n()), by = "species")

# three points where data was removed are :
## basis of record
## coordinate uncertainty and precision
## "zoo" in locality

# can't really compromise on the latter two, but maybe on basis of record
# definitely don't want fossil specimens, but can look at remaining samples by species

# coordinate uncertainty and precision
few.raw <- few.raw %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.02 | is.na(coordinatePrecision)) 
few.gbif <- left_join(few.gbif, summarise(group_by(few.raw, species), cucpn = n()), by = "species")

# zoo in locality
few.raw <- few.raw %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE)))
few.gbif <- left_join(few.gbif, summarise(group_by(few.raw, species), zoon = n()), by = "species")

# basis of record
few.raw <- few.raw %>%
  filter(basisOfRecord != "FOSSIL_SPECIMEN")
few.gbif <- left_join(few.gbif, summarise(group_by(few.raw, species), born = n()), by = "species")

# converting to spatial
few.raw <- few.raw %>%
  mutate(infraspecificEpithet = case_when(infraspecificEpithet == "" ~ NA_character_,
                                          TRUE ~ infraspecificEpithet)) %>%
  mutate(subgroup = infraspecificEpithet)

few.raw <- as.data.frame(few.raw)

few.raw.spatial <- SpatialPointsDataFrame(coords      = few.raw[, c("decimalLongitude", "decimalLatitude")],
                                          data        = few.raw[, names(few.raw)[!names(few.raw) %in% c("decimalLongitude", "decimalLatitude")]], 
                                          proj4string = CRS(proj4string(IUCN_Native_Data)))

# Acinonyx jubatus ####
curr.species <- "Acinonyx jubatus"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("hecki soemmeringii", "jubatus"), buff = c(1.8, 1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Aepyceros melampus ####
curr.species <- "Aepyceros melampus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]

# separate petersi and remove melampus sample
clip_poly <- Polygon(matrix(c(13.1, -17.4,
                              19.3, -17.4, 
                              19.3, -24.6, 
                              13.1, -24.6, 
                              13.1, -17.4), 
                            ncol = 2, byrow = TRUE))
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

petersi_points <- sp.gbif[clip_poly,]
petersi_points@data[is.na(petersi_points$subgroup), "subgroup"] <- "petersi"
petersi_points <- petersi_points[petersi_points$subgroup != "melampus", ]

# remove southern sa samples from melampus
poly_temp <- terra::buffer(sprp[sprp$subgroup == "melampus",],4)
poly_temp <- poly_temp - clip_poly

melampus_points <- sp.gbif[poly_temp,]
melampus_points@data[is.na(melampus_points$subgroup), "subgroup"] <- "melampus"
melampus_points <- melampus_points[melampus_points$subgroup != "petersi", ]

sp.gbif <- raster::bind(petersi_points, melampus_points)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$HostCorrectedName != curr.species,], sp.gbif)

# Alcelaphus buselaphus ####
curr.species <- "Alcelaphus buselaphus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]

# lichtensteinii
clip_poly <- countries50[countries50$name == "Malawi",]
clip_poly <- terra::buffer(clip_poly, 0)
sprp <- raster::bind(sprp, clip_poly)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "lichtensteinii"

# cokii x lelwel
clip_poly <- Polygon(matrix(c(37.5, 1.5,
                              39, 1.5,
                              39, -0.4,
                              37.5, -0.4,
                              37.5, 1.5), 
                            ncol = 2, byrow = TRUE))
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, clip_poly)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "cokii x lelwel"

# cokii
poly_temp <- gBuffer(sprp[sprp$subgroup == "cokii",], width = 0.2, byid = TRUE)
poly_temp <- poly_temp - sprp[sprp$subgroup == "cokii x lelwel",]
sprp <- raster::bind(sprp[sprp$subgroup != "cokii",], poly_temp)

## major = 2
## swaynei = 0.5
## caama = 0.5
## lichtensteinii = 0.1
## lelwel = 1
## cokii x lelwel = 0
## cokii = 0

buffer_temp <- data.frame(subgroup = c("caama", "cokii", "cokii x lelwel", "lelwel", "lichtensteinii", "major", "swaynei"), 
                          buff = c(0.5, 0, 0, 1, 0.1, 2, 0.5))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Alces alces ####
# buffer = 8
curr.species <- "Alces alces"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
# alces
# paper:A REVIEW OF CIRCUMPOLAR MOOSE POPULATIONS WITH EMPHASIS ON EURASIAN MOOSE DISTRIBUTIONS AND DENSITIES
# paper: Presence of moose (Alces alces) in Southeastern Germany

clip_poly <- Polygon(matrix(c(24.6, 52.9,
                              46.8, 52.9, 
                              46.8, 49, 
                              24.6, 49, 
                              24.6, 52.9), 
                            ncol = 2, byrow = TRUE))
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

clip_points <- matrix(c(82.53, 65.06),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 10000)

sprp <- raster::bind(sprp, clip_poly, clip_points)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "alces"

# buturlini cameloides pfizenmayeri
clip_points <- matrix(c(89.4, 68.8, 
                        170.2, 68.7, 
                        167.2, 61.5, 
                        159, 55.7),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 10000)

sprp <- raster::bind(sprp, clip_points)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "buturlini cameloides pfizenmayeri"

# pip_test
buffer_temp <- data.frame(subgroup = c("americana andersoni gigas shirasi", "alces", "buturlini cameloides pfizenmayeri"), 
                          buff = c(7, 1.5, 3.5))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Antidorcas marsupialis ####
curr.species <- "Antidorcas marsupialis"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif@data[sp.gbif$subgroup %in% c("angolensis", "hofmeyri"), "subgroup"] <- "angolensis hofmeyri"

# angolensis hofmeyri
sprp_angolensis <- gBuffer(sprp[sprp$subgroup == "angolensis hofmeyri",], width = 0.6, byid = TRUE)
sprp_marsupialis <- sprp[sprp$subgroup == "marsupialis",]

clip_points <- matrix(c(31.7, -25.28),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_angolensis <- raster::bind(sprp_angolensis,terra::buffer(clip_points, 100000))
sprp_angolensis@data[is.na(sprp_angolensis$subgroup), "subgroup"] <- "angolensis hofmeyri"
sprp_angolensis <- sprp_angolensis - sprp_marsupialis

angolensis_points <- sp.gbif[sprp_angolensis,]
angolensis_points@data[is.na(angolensis_points$subgroup), "subgroup"] <- "angolensis hofmeyri"

# marsupialis
sprp_marsupialis <- gBuffer(sprp_marsupialis, width = 2.5, byid = TRUE)
sprp_marsupialis <- sprp_marsupialis - sprp_angolensis

marsupialis_points <- sp.gbif[sprp_marsupialis,]
marsupialis_points@data[is.na(marsupialis_points$subgroup), "subgroup"] <- "marsupialis"

sp.gbif <- raster::bind(angolensis_points, marsupialis_points)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Antilocapra americana ####
curr.species <- "Antilocapra americana"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "americana sonoriensis"

# Bison bison ####
curr.species <- "Bison bison"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "bison athabascae"

# Bison bonasus ####
curr.species <- "Bison bonasus"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "bonasus caucasicus"

# Blastocerus dichotomus ####
curr.species <- "Blastocerus dichotomus"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Canis adustus ####
curr.species <- "Canis adustus"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Canis aureus ####
curr.species <- "Canis aureus"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# aureus syriacus
sprp_aureus <- sprp[sprp$subgroup == "aureus syriacus",]
clip_points <- matrix(c(45.3, 47.8,
                        43.8, 49.2),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))

sprp_aureus <- raster::bind(sprp_aureus, terra::buffer(clip_points, 10000))
sprp_aureus@data$subgroup <- "aureus syriacus"
sprp_aureus <- gBuffer(sprp_aureus, width = 0.2, byid = TRUE)

# indicus naria cruesmanni
sprp_indicus <- sprp[sprp$subgroup == "indicus naria cruesmanni",]
sprp_indicus <- gBuffer(sprp_indicus, width = 2.1, byid = TRUE)
sprp_indicus <- sprp_indicus - sprp_aureus

# moreotica ecsedensis
# Not going beyond 2 because I don't know if widespread samples are established populations or migratory
sprp_moreotica <- sprp[sprp$subgroup == "moreotica ecsedensis",]
sprp_moreotica <- gBuffer(sprp_moreotica, width = 2, byid = TRUE)
sprp_moreotica <- sprp_moreotica - sprp_aureus

sprp_try <- raster::bind(sprp_aureus, sprp_indicus, sprp_moreotica)

# pip_test
buffer_temp <- data.frame(subgroup = c("aureus syriacus", "indicus naria cruesmanni", "moreotica ecsedensis"), 
                          buff = c(0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Canis latrans ####
curr.species <- "Canis latrans"
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Canis lupus ####
curr.species <- "Canis lupus"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

sp.gbif <- sp.gbif[is.na(sp.gbif$subgroup) | !sp.gbif$subgroup %in% c("familiaris", "dingo"),]
sp.gbif <- sp.gbif[!str_detect(sp.gbif$verbatimScientificName, "familiaris"),]

sprp_crassodon <- gBuffer(sprp[sprp$subgroup == "crassodon",], width = 0.1, byid = TRUE)
sprp_lycaon <- gBuffer(sprp[sprp$subgroup == "lycaon occidentalis nubilis",], width = 3, byid = TRUE) - sprp_crassodon
sprp_signatus <- gBuffer(sprp[sprp$subgroup == "signatus",], width = 1.4, byid = TRUE)
sprp_pallipes <- gBuffer(sprp[sprp$subgroup == "pallipes chanco",], width = 0.1, byid = TRUE)
sprp_lupus <- gBuffer(sprp[sprp$subgroup == "lupus",], width = 1.06, byid = TRUE) - sprp_pallipes - sprp[sprp$subgroup == "italicus",]
sprp_italicus <- gBuffer(sprp[sprp$subgroup == "italicus",], width = 3, byid = TRUE) - sprp_lupus - sprp_signatus

sprp_try <- raster::bind(sprp[sprp$subgroup %in% c("baileyi", "arctos", "arabs"),],
                         sprp_crassodon, sprp_lycaon, sprp_signatus, sprp_pallipes, sprp_lupus, sprp_italicus)

# pip_test
buffer_temp <- data.frame(subgroup = c("baileyi", "arctos", "crassodon", 
                                       "lycaon occidentalis nubilis", "signatus", "arabs",
                                       "pallipes chanco", "lupus", "italicus"), 
                          buff = c(3, 3, 0, 0, 0, 0, 0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = sp.gbif, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Canis mesomelas ####
curr.species <- "Canis mesomelas"
# map # plot_temp(curr.species)
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("mesomelas", "schmidti"), buff = c(3, 1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Capra ibex ####
curr.species <- "Capra ibex"
plot_temp(curr.species)
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("not used"), buff = c(1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Capra pyrenaica ####
curr.species <- "Capra pyrenaica"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("victoriae", "hispanica"), buff = c(0.45, 0.5)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Capreolus capreolus ####
curr.species <- "Capreolus capreolus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

sprp_caucasicus <- gBuffer(sprp[sprp$subgroup == "caucasicus",], width = 1, byid = TRUE)
sprp_coxi <- gBuffer(sprp[sprp$subgroup == "coxi",], width = 0.5, byid = TRUE)

sprp_capreolus <- gBuffer(sprp[sprp$subgroup == "italicus garganta capreolus",], width = 4, byid = TRUE)
sprp_capreolus <- sprp_capreolus - sprp_caucasicus - sprp_coxi

sprp_try <- raster::bind(sprp_caucasicus, sprp_coxi, sprp_capreolus)

# pip_test
buffer_temp <- data.frame(subgroup = c("caucasicus", "coxi", "italicus garganta capreolus"), 
                          buff = c(0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Capricornis crispus ####
curr.species <- "Capricornis crispus"
# map # plot_temp(curr.species)
# 
# GBIF_Spatial <- GBIF_Spatial[(GBIF_Spatial$species != curr.species | 
#                                 GBIF_Spatial$countryCode != "TWN"),]

GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Cephalophus natalensis ####
curr.species <- "Cephalophus natalensis"
# map # plot_temp(curr.species)

GBIF_Spatial <- GBIF_Spatial[GBIF_Spatial$verbatimScientificName != "Cephalophus harveyi",]

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("natalensis", "robertsi"), buff = c(1, 0.5)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Cerdocyon thous ####
curr.species <- "Cerdocyon thous" 
# map # plot_temp(curr.species)
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("azarae entrerianus", "thous aquilus germanus"), buff = c(3, 3)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Cervus canadensis and elaphus ####
curr.species <- "Cervus elaphus"
# map # plot_temp(curr.species) # crashes R, too many points

gbif_sp_temp <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]

## elaphus (species)
clip_poly <- Polygon(matrix(c(-24, 65,
                              25, 71,
                              34, 53,
                              60, 50,
                              60, 30,
                              -24, 30,
                              -24, 65), 
                            ncol = 2, byrow = TRUE))
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
elaphus_gbif <- gbif_sp_temp[clip_poly,]

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
# atlanticus (subspecies)
clip_points <- matrix(c(17, 68,
                        20, 66,
                        17, 65,
                        17, 63,
                        17, 60),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 200000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "atlanticus"

clip_points <- matrix(c(12, 55.47,
                        12, 55.2),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 50000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "atlanticus"

# hispanicus (subspecies)
clip_points <- matrix(c(-8, 41,
                        -8.5, 39),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- raster::bind(sprp[sprp$subgroup == "hispanicus", ], 
                            terra::buffer(clip_points, 50000))
clip_points <- terra::buffer(clip_points, 0.7)

clip_points <- clip_points - terra::buffer(countries50[countries50$name == "France",], 0.5)
poly_temp <- terra::buffer(countries50[countries50$name %in% c("Spain", "Andorra"),], 0.06)
clip_points <- gUnion(clip_points, poly_temp)

sprp <- raster::bind(sprp[sprp$subgroup != "hispanicus",], clip_points)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "hispanicus"

# italicus (subspecies)
clip_points <- matrix(c(14,41,
                        18,40,
                        14,38),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 25000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "italicus"

# scoticus (subspecies)
clip_points <- matrix(c(-6.5, 58,
                        -7, 57,
                        
                        -9, 54,
                        -8.5, 53,
                        -8, 53,
                        -10, 52,
                        -8.5, 52,
                        -7, 52,
                        
                        -4.5, 52,
                        -4, 53),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 30000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "scoticus"

# elaphus (subspecies)
clip_points <- matrix(c(4, 51.5, 
                        2.5, 51,
                        7.5, 53,
                        
                        8, 44,
                        7, 43.5,
                        6, 43,
                        
                        22, 37,
                        20, 40),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 5000)
clip_points <- raster::bind(sprp[sprp$subgroup == "elaphus",], clip_points)

clip_points <- terra::buffer(clip_points)
clip_points <- raster::bind(clip_points, terra::buffer(countries50[countries50$name == "Ukraine",]))

clip_points <- (clip_points - terra::buffer(sprp[sprp$subgroup == "atlanticus",], 0.3) - 
                  terra::buffer(sprp[sprp$subgroup == "brauneri",], 0.6) - 
                  sprp[sprp$subgroup == "hispanicus",] -
                  terra::buffer(sprp[sprp$subgroup == "italicus",], 0.7)-
                  terra::buffer(sprp[sprp$subgroup == "maral",], 0.3) - 
                  terra::buffer(sprp[sprp$subgroup == "scoticus",], 0.6))

sprp <- raster::bind(sprp[sprp$subgroup != "elaphus",], clip_points)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "elaphus"

# pip_test
buffer_temp <- data.frame(subgroup = c("atlanticus", "barbarus", "brauneri", "corsicanus", "elaphus", "hispanicus", "italicus", "maral", "scoticus"), 
                          buff = c(0.3, 0.5, 0.6, 0.8, 0, 0, 0.7, 0.3, 0.6))
sp.gbif <- pip_test(curr.species, dat = elaphus_gbif, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- GBIF_Spatial[GBIF_Spatial$species != curr.species,]
GBIF_Spatial <- raster::bind(GBIF_Spatial, sp.gbif)

# canadensis
sp.gbif <- pip_test(curr.species, dat = gbif_sp_temp, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == "Cervus canadensis",], 
                    buff = data.frame(subgroup = c("other", "canadensis nannodes roosevelti"), buff = c(5, 15)))

GBIF_Spatial <- raster::bind(GBIF_Spatial, sp.gbif)

# Cervus nippon ####
curr.species <- "Cervus nippon"
# map # plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# aplodontus
clip_points <- matrix(c(140.5, 36),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 5000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "aplodontus"

# nippon
poly_temp <- terra::buffer(sprp[sprp$subgroup == "nippon",], 0.1)
poly_temp <- poly_temp - terra::buffer(sprp[sprp$subgroup == "aplodontus",], 0.1)
sprp <- raster::bind(sprp[sprp$subgroup != "nippon",], poly_temp)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "nippon"

# other (taioanus)
clip_points <- matrix(c(121, 24),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 110000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "other"

# pip_test
buffer_temp <- data.frame(subgroup = c("aplodontus", "nippon", "other", "yesoensis"), 
                          buff = c(0.2, 0, 0.3, 0.1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Chrysocyon brachyurus ####
curr.species <- "Chrysocyon brachyurus"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Civettictis civetta ####
curr.species <- "Civettictis civetta"
# map # plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Conepatus chinga ####
curr.species <- "Conepatus chinga"
plot_temp(curr.species)

# ones further south could be misidentified Conepatus humboldtii 
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 2),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Connochaetes gnou ####
curr.species <- "Connochaetes gnou"
plot_temp(curr.species)

# samples in Namibia may be hybrids with C. taurinus
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 4),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Connochaetes taurinus ####
curr.species <- "Connochaetes taurinus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# albojubatus
clip_points <- matrix(c(38.49, -3.6,
                        38.738527, -1.211118),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 10000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "albojubatus"

# mearnsi
# decided to change western samples labelled as albojubatus assuming they are misidentified 
# mainly bc of location 

clip_points <- matrix(c( 34.379701, -1.607878,
                         36.088075, -0.361138, 
                         36.349000, -0.756626 ),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 10000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "mearnsi"

# pip_test
# assuming samples in coastal SA are artificially managed populations or hybrid with C. gnou
buffer_temp <- data.frame(subgroup = c("albojubatus", "cooksoni", "johnstoni", "mearnsi", "taurinus"), 
                          buff = c(0.2, 0.1, 0, 0.13, 2.5))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Crocuta crocuta ####
curr.species <- "Crocuta crocuta"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 4),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Cynictis penicillata ####
curr.species <- "Cynictis penicillata"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Damaliscus lunatus ####
curr.species <- "Damaliscus lunatus"
plot_temp(curr.species)

# pip_test
# assuming samples in coastal SA are artificially managed populations
buffer_temp <- data.frame(subgroup = c("jimela", "korrigum", "lunatus", "superstes", "tiang", "topi"), 
                          buff = c(2, 0, 1, 0, 0, 1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Damaliscus pygargus ####
curr.species <- "Damaliscus pygargus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# phillipsi
clip_points <- matrix(c(31.871424, -28.788270, 
                        29.907618, -30.293607),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 20000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "phillipsi"

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup == c("phillipsi", "pygargus"), buff = c(0.9, 1.6)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Equus grevyi ####
curr.species <- "Equus grevyi"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Equus quagga ####
curr.species <- "Equus quagga"
plot_temp(curr.species)

# Namibian and South African samples outside of range are likely Equus zebra
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species ==  curr.species,]

sprp <- gBuffer(sprp, width = 1, byid = TRUE)
clip_poly <- Polygon(matrix(c(14.38, -20.19, 
                              19.13, -18.63, 
                              19.13, -25.31, 
                              14.38, -25.31, 
                              14.38, -20.19), 
                            ncol = 2, byrow = TRUE))
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- sprp - clip_poly

sp.gbif <- sp.gbif[sprp,]
sp.gbif@data$subgroup <- "not used"

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Equus zebra ####
curr.species <- "Equus zebra"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("zebra", "hartmannae"), buff = c(1, 1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Felis silvestris ####
curr.species <- "Felis silvestris"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

unique(IUCN_Native_Data@data[IUCN_Native_Data$binomial == "Felis lybica", "subgroup"])
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species &
                          (!GBIF_Spatial$subgroup %in% c("lybica", "ornata", "cafra")|is.na(GBIF_Spatial$subgroup)),]

# silvestris
sprp_silvestris <- gBuffer(sprp[sprp$subgroup == "silvestris",], byid = TRUE, width = 2.5)

# caucasica
sprp_caucasica <- gBuffer(sprp[sprp$subgroup == "caucasica",], byid = TRUE, width = 1)
sprp_caucasica <- sprp_caucasica - sprp_silvestris
sprp <- raster::bind(sprp[sprp$subgroup == "grampia",], sprp_caucasica, sprp_silvestris)

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("caucasica", "grampia", "silvestris"), buff = c(0, 0, 0)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)


# Genetta genetta ####
curr.species <- "Genetta genetta"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# afra
clip_points <- matrix(c(-3, 31),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 30000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "afra"

# senegalensis
clip_points <- matrix(c(2, 8),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 100000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "senegalensis"

# pip_test
buffer_temp <- data.frame(subgroup = c("afra", "dongolana", "felina", "senegalensis"), 
                          buff = c(0, 0.5, 1, 2))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Genetta thierryi ####
curr.species <- "Genetta thierryi"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif.raw <- few.raw.spatial[few.raw.spatial$species == curr.species,]

# tm_shape(sprp) + tm_polygons("subgroup") + 
#   tm_shape(sp.gbif.raw) + tm_dots() + 
#   tm_shape(sp.gbif) + tm_dots(col = "blue") 
# checked samples manually on gbif or on host sites and they all seem legit apart from one in DEU

sp.gbif <- sp.gbif.raw[terra::buffer(sprp, 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Giraffa camelopardalis ####
curr.species <- "Giraffa camelopardalis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# rothschildi
sprp_roths <- gBuffer(sprp[sprp$subgroup == "rothschildi",], width = 0.7, byid = TRUE)
sprp_roths <- sprp_roths - sprp[sprp$subgroup == "reticulata",] - sprp[sprp$subgroup == "tippelskirchi",]

clip_points <- matrix(c(37.5, 0.15),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_roths <- raster::bind(sprp_roths, terra::buffer(clip_points, 20000))
sprp_roths@data[is.na(sprp_roths$subgroup), "subgroup"] <- "rothschildi"

clip_points <- matrix(c(37.7, 0.5),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_roths <- sprp_roths - terra::buffer(clip_points, 20000)

# reticulata and tippelskirchi
sprp_retic <- gBuffer(sprp[sprp$subgroup == "reticulata",], width = 0.7, byid = TRUE)
sprp_retic <- sprp_retic - sprp_roths - terra::buffer(sprp[sprp$subgroup == "tippelskirchi",], 0.1)

sprp_tipp <- gBuffer(sprp[sprp$subgroup == "tippelskirchi",], width = 0.5, byid = TRUE)
sprp_tipp <- sprp_tipp - sprp_roths - sprp_retic

sprp_try <- raster::bind(sprp_roths, sprp_retic, sprp_tipp, sprp[sprp$subgroup == "giraffa",])

# pip_Test
buffer_temp <- data.frame(subgroup = c("giraffa", "reticulata", "rothschildi", "tippelskirchi"), 
                          buff = c(0.5, 0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Herpestes ichneumon ####
curr.species <- "Herpestes ichneumon"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("ichneumon", "sangronizi numidicus widdringtonii"), buff = c(5, 3)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Hippopotamus amphibius ####
curr.species <- "Hippopotamus amphibius"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 4),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Hippotragus niger ####
curr.species <- "Hippotragus niger"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species ,], 0.2),]
sp.gbif@data$subgroup <- "niger kirkii roosevelti"

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Hyaena hyaena ####
curr.species <- "Hyaena hyaena"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Kobus ellipsiprymnus ####
curr.species <- "Kobus ellipsiprymnus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

sprp_defassa <- gBuffer(sprp[sprp$subgroup == "defassa",], width = 1, byid = TRUE)
sprp_defassa <- sprp_defassa - sprp[sprp$subgroup == "ellipsiprymnus",]

sprp_ellipsiprymnus <- gBuffer(sprp[sprp$subgroup == "ellipsiprymnus",], width = 1, byid = TRUE)
sprp_ellipsiprymnus <- sprp_ellipsiprymnus - sprp[sprp$subgroup == "defassa",]

sprp_try <- raster::bind(sprp_ellipsiprymnus, sprp_defassa)

# pip_test
# dropping most points in SA and Namibia because they are likely managed for hunting 
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = data.frame(subgroup = c("defassa", "ellipsiprymnus"), 
                                      buff = c(0, 0)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Kobus leche ####
curr.species <- "Kobus leche"
plot_temp(curr.species)

buffer_temp <- data.frame(subgroup = c("kafuensis", "leche", "smithemani"), 
                          buff = c(0, 0.5, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Leopardus geoffroyi ####
curr.species <- "Leopardus geoffroyi"
plot_temp(curr.species)

GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Leopardus pardalis ####
curr.species <- "Leopardus pardalis"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# mitis
sprp_mitis <- gBuffer(sprp[sprp$subgroup == "mitis",], width = 2, byid = TRUE)
sprp_mitis <- sprp_mitis - terra::buffer(sprp[sprp$subgroup == "pardalis",], 0.1)

# pardalis
clip_points <- matrix(c(-103.3, 29.3),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))

sprp <- raster::bind(sprp_mitis, sprp[sprp$subgroup == "pardalis",], terra::buffer(clip_points, 1000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "pardalis"

# Leopardus tigrinus ####
curr.species <- "Leopardus tigrinus"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 5),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Lontra canadensis ####
curr.species <- "Lontra canadensis"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# mainland
clip_points <- matrix(c(-88, 47),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 30000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "canadensis pacifica sonora lataxina mira"

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("periclyzomae", "canadensis pacifica sonora lataxina mira"), buff = c(0.5, 0.5)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Lutra lutra ####
curr.species <- "Lutra lutra"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
poly_temp <- sprp[sprp$subgroup == "lutra", ]
poly_temp <- gBuffer(poly_temp, byid = TRUE, width = 5)
poly_temp <- poly_temp - terra::buffer(sprp[sprp$subgroup == "angustifrons", ], 1)
poly_temp <- poly_temp - sprp[sprp$subgroup == "meridionalis", ]
sprp <- raster::bind(sprp[sprp$subgroup != "lutra",], poly_temp)

# pip_test
buffer_temp <- data.frame(subgroup = c("chinensis", "aurobrunnea kutab monticolus", "angustifrons", "lutra", "meridionalis"), 
                          buff = c(0.5, 0, 1, 0, 0))

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Lycalopex gymnocercus ####
curr.species <- "Lycalopex gymnocercus"
plot_temp(curr.species)

# ones to the west and south could be misidentified L. culpaeus
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1.5),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Lycaon pictus ####
curr.species <- "Lycaon pictus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif.raw <- few.raw.spatial[few.raw.spatial$species == curr.species,]

# tm_shape(sprp) + tm_polygons("subgroup") +
#   tm_shape(sp.gbif.raw) + tm_dots("basisOfRecord") 
# checked samples manually on gbif or on host sites and they all seem legit 

sp.gbif.raw@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif.raw)

# Lynx canadensis ####
curr.species <- "Lynx canadensis"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Lynx lynx ####
curr.species <- "Lynx lynx"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# lynx
sprp_lynx <- sprp[sprp$subgroup == "lynx",]
sprp_lynx <- gBuffer(sprp_lynx, width = 1.2, byid = TRUE)
sprp_lynx <- sprp_lynx - terra::buffer(sprp[sprp$subgroup == "other",], 0.1)

# other
sprp_other <- sprp[sprp$subgroup == "other",]
sprp_other <- gBuffer(sprp_other, width = 3, byid = TRUE)
sprp_other <- sprp_other - sprp_lynx

sprp <- raster::bind(sprp_other, sprp_lynx, sprp[sprp$subgroup == "balcanicus",])

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("balcanicus", "lynx", "other"), buff = c(1, 0, 0)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Lynx pardinus ####
curr.species <- "Lynx pardinus"

# http://dx.doi.org/10.13140/RG.2.2.12500.94087
# keeping all samples because they are from recent historic range

GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Lynx rufus ####
curr.species <- "Lynx rufus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sprp_fasciatus <- sprp[sprp$subgroup == "fasciatus",]
sprp_rufus <- sprp[sprp$subgroup == "rufus",]
sprp_mexican <- sprp[sprp$subgroup == "mexican group",]

# fasciatus
sprp_fasciatus <- gBuffer(sprp_fasciatus, width = 1, byid = TRUE)
sprp_fasciatus <- sprp_fasciatus - sprp_rufus - sprp_mexican

# mexican
sprp_mexican <- gBuffer(sprp_mexican, width = 1, byid = TRUE)
sprp_mexican <- sprp_mexican - sprp[sprp$subgroup == "fasciatus",] - sprp_rufus

# rufus
sprp_rufus <- gBuffer(sprp_rufus, width = 2, byid = TRUE)
sprp_rufus <- sprp_rufus - sprp[sprp$subgroup == "fasciatus",] - sprp[sprp$subgroup == "mexican group",]

clip_points <- matrix(c(-97.504, 25.377),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_rufus <- sprp_rufus - terra::buffer(clip_points, 10000)

sprp <- raster::bind(sprp_fasciatus, sprp_mexican, sprp_rufus)

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("fasciatus", "mexican group", "rufus"), buff = c(0, 0, 0)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Martes americana ####
curr.species <- "Martes americana"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Martes foina ####
curr.species <- "Martes foina"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# nehringi syriaca
sprp_nehringi <- gBuffer(sprp[sprp$subgroup == "nehringi syriaca",], width = 1, byid = TRUE)
clip_points <- matrix(c(45.496,35.404),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_nehringi <- raster::bind(sprp_nehringi, terra::buffer(clip_points, 100000))
sprp_nehringi@data[is.na(sprp_nehringi$subgroup), "subgroup"] <- "nehringi syriaca"

sprp_nehringi <- sprp_nehringi - sprp[sprp$subgroup %in% c("milleri", "western"),] 

# western
sprp_western <- gBuffer(sprp[sprp$subgroup == "western",], width = 5, byid = TRUE)
sprp_western <- sprp_western - sprp_nehringi - sprp[sprp$subgroup %in% c("milleri", "bunites", "rosanowi"),]

sprp_try <- raster::bind(sprp_nehringi, sprp_western, sprp[sprp$subgroup %in% c("eastern", "milleri", "bunites", "rosanowi"),])

# pip_Test
buffer_temp <- data.frame(subgroup = c("western", "bunites", "milleri", 
                                       "nehringi syriaca", "rosanowi", "eastern"), 
                          buff = c(0, 0, 0, 
                                   0, 0, 1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Martes martes ####
curr.species <- "Martes martes"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Martes melampus ####
curr.species <- "Martes melampus"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 0.5),]
sp.gbif@data$subgroup <- "melampus"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Martes pennanti/Pekania pennanti ####
curr.species <- "Pekania pennanti"
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == "Martes pennanti",]) + tm_polygons("subgroup") +
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species,]) + tm_dots()
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Meles meles ####
curr.species <- "Meles meles"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

sprp_marianensis <- gBuffer(sprp[sprp$subgroup == "marianensis",], width = 0.05, byid = TRUE)
sprp_marianensis <- sprp_marianensis - sprp[sprp$subgroup == "meles milleri heptneri",]

sprp_meles <- gBuffer(sprp[sprp$subgroup == "meles milleri heptneri",], width = 1.5, byid = TRUE)
sprp_meles <- sprp_meles - sprp_marianensis

sprp_try <- raster::bind(sprp_meles, sprp_marianensis)

# pip_Test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = data.frame(subgroup = c("marianensis", "meles milleri heptneri"), buff = c(0, 0)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Melogale moschata/subaruantiaca ####
curr.species <- "Melogale moschata"
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == "Melogale subaurantiaca",]) + tm_polygons("subgroup") +
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species,]) + tm_dots()

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 0.5),]
sp.gbif@data$subgroup <- "subaurantiaca"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Mephitis mephitis ####
curr.species <- "Mephitis mephitis"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 2),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Mustela erminea ####
curr.species <- "Mustela erminea"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# arctica
sprp_arctica <- gBuffer(sprp[sprp$subgroup == "arctica",], width = 2.7, byid = TRUE)
sprp_arctica <- sprp_arctica - sprp[sprp$subgroup == "kadiacensis",]

# hibernica
sprp_hibernica <- gBuffer(sprp[sprp$subgroup == "hibernica",], width = 0.1, byid = TRUE)
clip_points <- matrix(c(-4.5441,54.2216),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_hibernica <- raster::bind(sprp_hibernica, terra::buffer(clip_points, 20000))
sprp_hibernica@data[is.na(sprp_hibernica$subgroup), "subgroup"] <- "hibernica"

# stabilis
sprp_stabilis <- gBuffer(sprp[sprp$subgroup == "stabilis",], width = 1, byid = TRUE)
sprp_stabilis <- sprp_stabilis - sprp[sprp$subgroup == "aestiva",] - sprp_hibernica

# aestiva
sprp_aestiva <- gBuffer(sprp[sprp$subgroup == "aestiva",], width = 1.5, byid = TRUE)
sprp_aestiva <- sprp_aestiva - sprp_stabilis - sprp[sprp$subgroup %in% c("minima", "kaneii and other", "erminea"),]

sprp_try <- raster::bind(sprp[sprp$subgroup %in% c("kadiacensis", "polaris", "kaneii and other", "erminea", "minima"),],
                         sprp_arctica, sprp_hibernica, sprp_stabilis, sprp_aestiva)

# pip_Test
buffer_temp <- data.frame(subgroup = c("arctica", "kadiacensis", "polaris", 
                                       "kaneii and other", "erminea", "hibernica",
                                       "stabilis", "aestiva", "minima"), 
                          buff = c(0, 0, 0, 
                                   0, 0, 0,
                                   0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)


# Mustela lutreola ####
curr.species <- "Mustela lutreola"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Mustela nivalis ####
curr.species <- "Mustela nivalis"
plot_temp(curr.species)

# numidica
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sprp <- raster::bind(sprp, terra::buffer(countries50[countries50$name == "Egypt",], 0))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "numidica"

# eurasian
clip_points <- matrix(c(-6.412437, 57.395929,
                        11.319496, 55.264806, 
                        15.254495, 68.522656),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 100000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "eurasian"

# pip_test
buffer_temp <- data.frame(subgroup = c("american", "numidica", "formosana", "eurasian"), 
                          buff = c(0, 0, 0.3, 0.3))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Mustela putorius ####
curr.species <- "Mustela putorius"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# removing domesic ferret as much as possible
sp.gbif <- sp.gbif[!str_detect(sp.gbif$verbatimScientificName, "furo"),]

# anglia caledoniae
# https://www.vwt.org.uk/wp-content/uploads/2016/04/Polecat-Report-2016.pdf
# roughly following distribution from above as it's described as verifiable records or true polecats
# essentially only including mainland England, Wales, and Scotland samples

sprp_anglia <- raster::disaggregate(countries50[countries50$name == "United Kingdom",])
clip_points <- matrix(c(-1.692,53.958),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_anglia <- sprp_anglia[clip_points,]
sprp_anglia@data <- sprp@data[sprp$subgroup == "anglia caledoniae",]
sprp_anglia <- gBuffer(sprp_anglia, width = 0.25, byid = TRUE)

# rothschildi
sprp_rothschildi <- gBuffer(sprp[sprp$subgroup == "rothschildi",], width = 0.1, byid = TRUE)

# putorius and european
sprp_putorius <- gBuffer(sprp[sprp$subgroup == "putorius and european",], width = 1.3, byid = TRUE)
sprp_putorius <- sprp_putorius - sprp_rothschildi - sprp_anglia

sprp_try <- raster::bind(sprp_anglia, sprp_rothschildi, sprp_putorius)

# pip_test
buffer_temp <- data.frame(subgroup = c("anglia caledoniae", "rothschildi", "putorius and european"), 
                          buff = c(0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = sp.gbif, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Mustela vison/Neovison vison ####
curr.species <- "Mustela vison"
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == "Neovison vison",]) + tm_polygons("subgroup") +
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species,]) + tm_dots()

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == "Neovison vison",], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Nyctereutes procyonoides ####
curr.species <- "Nyctereutes procyonoides"
plot_temp(curr.species)

buffer_temp <- data.frame(subgroup = c("albus", "other", "viverrinus"), 
                          buff = c(0.5, 0.5, 0.3))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Odocoileus hemionus ####
curr.species <- "Odocoileus hemionus"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 6),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Odocoileus virginianus ####
curr.species <- "Odocoileus virginianus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# northern
states_temp <- c("Utah", "California", "Idaho", "British Columbia", "New Brunswick", "Newfoundland and Labrador")
clip_poly <- terra::buffer(states50[states50$name %in% states_temp,], 0)
sprp <- raster::bind(sprp, clip_poly)
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "northern"

# southern
clip_points <- matrix(c(-78.9, -1.2),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 100000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "southern"

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("northern", "southern"), buff = c(1, 1.4)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Otocyon megalotis ####
curr.species <- "Otocyon megalotis"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("megalotis", "virgatus"), buff = c(1, 1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Ourebia ourebi ####
curr.species <- "Ourebia ourebi"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "oribi"

# Ovibos moschatus ####
curr.species <- "Ovibos moschatus"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 9),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Ovis canadensis ####
curr.species <- "Ovis canadensis"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Ovis dalli ####
curr.species <- "Ovis dalli"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 2),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Ozotoceros bezoarticus ####
curr.species <- "Ozotoceros bezoarticus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif.raw <- few.raw.spatial[few.raw.spatial$species == curr.species,]

# tm_shape(sprp) + tm_polygons("subgroup") +
#   tm_shape(sp.gbif.raw) + tm_dots("basisOfRecord")
# no weird outliers

# adding to celer to catch point
clip_points <- matrix(c(-62.108,-32.974),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 10000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "celer"

# pip_test
buffer_temp <- data.frame(subgroup = c("leucogaster", "bezoarticus", "arerunguaensis", "uruguayensis", "celer"), 
                          buff = c(2.5, 3.5, 0, 0, 2.5))

sp.gbif <- pip_test(curr.species, dat = sp.gbif.raw, 
                    range.polygon = sprp, 
                    buff = buffer_temp)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Panthera leo ####
curr.species <- "Panthera leo"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("leo", "persica"), buff = c(5, 1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Panthera onca ####
curr.species <- "Panthera onca"
plot_temp(curr.species)

# central
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
clip_points <- matrix(c(-79, -2.5,
                        -79.5, -0.5),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 10000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "central"

# pip_test
buffer_temp <- data.frame(subgroup = c("central", "northern", "southern"), 
                          buff = c(0.5, 0.5, 0.5))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Panthera pardus ####
curr.species <- "Panthera pardus"
plot_temp(curr.species)

buffer_temp <- data.frame(subgroup = c("delacouri", "fusca", "kotiya", "nimr", "pardus"), 
                          buff = c(0, 0, 0, 0.5, 0.5))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Pecari tajacu ####
curr.species <- "Pecari tajacu"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
# northern
poly_temp <- sprp[sprp$subgroup == "northern",]
poly_temp <- gBuffer(poly_temp, width = 4, byid = TRUE)
poly_temp <- poly_temp - terra::buffer(sprp[sprp$subgroup == "southern",])

sprp <- raster::bind(sprp[sprp$subgroup != "northern",], poly_temp)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("northern", "southern"), buff = c(0, 1)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Pelea capreolus ####
curr.species <- "Pelea capreolus"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 2),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Philantomba monticola ####
curr.species <- "Philantomba monticola"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Procapra gutturosa ####
curr.species <- "Procapra gutturosa"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif.raw <- few.raw.spatial[few.raw.spatial$species == curr.species,]

# tm_shape(sprp) + tm_polygons("subgroup") +
#   tm_shape(sp.gbif.raw) + tm_dots("basisOfRecord") 
# no weird outliers 

sp.gbif.raw@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif.raw)

# Procyon lotor ####
curr.species <- "Procyon lotor"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# IUCN polygon is missing Vancouver island
sprp_vancouverensis <- raster::disaggregate(countries50[countries50$name == "Canada",])
clip_points <- matrix(c(-125.5640,49.6982),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_vancouverensis <- sprp_vancouverensis[clip_points,]
sprp_vancouverensis@data <- sprp@data[sprp$subgroup == "grinnelli",]
sprp_vancouverensis <- gBuffer(sprp_vancouverensis, width = 0.15, byid = TRUE)

clip_points <- matrix(c(-123.0869,48.5379),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_vancouverensis <- raster::bind(sprp_vancouverensis, terra::buffer(clip_points, 15000))
sprp_vancouverensis@data$subgroup <- "vancouverensis"

# inesperatus 
sprp_inesperatus <- gBuffer(sprp[sprp$subgroup == "inesperatus",], width = 0.1, byid = TRUE)

# mainland
sprp_mainland <- gBuffer(sprp[sprp$subgroup == "mainland",], width = 2.1, byid = TRUE)
clip_points <- matrix(c(-78.388,26.622,
                        -78.207,7.881),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 100000)

sprp_mainland <- sprp_mainland - sprp_vancouverensis - sprp_inesperatus - clip_points - sprp[sprp$subgroup == "grinnelli",]

sprp_try <- raster::bind(sprp_vancouverensis, sprp_inesperatus, sprp_mainland, sprp[sprp$subgroup == "grinnelli",])

# pip_test
buffer_temp <- data.frame(subgroup = c("vancouverensis", "inesperatus", "mainland", "grinnelli"), 
                          buff = c(0, 0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Puma concolor ####
curr.species <- "Puma concolor"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
# couguar
poly_temp <- sprp[sprp$subgroup == "couguar",]
poly_temp <- gBuffer(poly_temp, width = 2.5, byid = TRUE)
poly_temp <- poly_temp - terra::buffer(sprp[sprp$subgroup == "concolor",])

sprp <- raster::bind(sprp[sprp$subgroup != "couguar",], poly_temp)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("couguar", "concolor"), buff = c(0, 1)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Rangifer tarandus ####
curr.species <- "Rangifer tarandus"
plot_temp(curr.species)
# avoiding samples outside of iucn range in Scandinavia because of semi-domestication

# pip_test
buffer_temp <- data.frame(subgroup = c("american", "other", "tarandus"), 
                          buff = c(3, 1.8, 0.8))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Raphicerus campestris ####
curr.species <- "Raphicerus campestris"
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("campestris", "neumanni"), buff = c(1, 1)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Redunca arundinum ####
curr.species <- "Redunca arundinum"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 0.7),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Redunca fulvorufula ####
curr.species <- "Redunca fulvorufula"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif.raw <- few.raw.spatial[few.raw.spatial$species == curr.species,]

# tm_shape(sprp) + tm_polygons("subgroup") +
#   tm_shape(sp.gbif.raw) + tm_dots("basisOfRecord") 
# no weird outliers 

sp.gbif <- pip_test(curr.species, dat = sp.gbif.raw, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("chanleri", "fulvorufula"), buff = c(1.5, 0.5)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Rupicapra pyrenaica ####
curr.species <- "Rupicapra pyrenaica"
plot_temp(curr.species)

buffer_temp <- data.frame(subgroup = c("ornata", "parva", "pyrenaica"), buff = c(0.5, 0.5, 1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Rupicapra rupicapra ####
curr.species <- "Rupicapra rupicapra"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# rupicapra cartusiana
sprp_rupicapra <- gBuffer(sprp[sprp$subgroup == "rupicapra cartusiana",], width = 2.5, byid = TRUE)
clip_points <- matrix(c(7.384,50.520,
                        16.684,49.334,
                        16.847,50.211),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)
sprp_rupicapra <- sprp_rupicapra - clip_points - terra::buffer(sprp[sprp$subgroup == "balcanica",], 0.2)

# balcanica
sprp_balcanica <- gBuffer(sprp[sprp$subgroup == "balcanica",], width = 0.7, byid = TRUE)
sprp_balcanica <- sprp_balcanica - sprp_rupicapra

sprp_try <- raster::bind(sprp_rupicapra, sprp_balcanica, sprp[sprp$subgroup %in% c("asiatica", "carpatica", "caucasica", "tatrica"),])

# pip_test
buffer_temp <- data.frame(subgroup = c("asiatica", "balcanica", "carpatica", "caucasica", "rupicapra cartusiana", "tatrica"), 
                          buff = c(0.1, 0, 0.7, 0.1, 0, 0.1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Sus scrofa ####
curr.species <- "Sus scrofa"
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

sp.gbif <- sp.gbif[!str_detect(sp.gbif$verbatimScientificName, "domestic"),]
sp.gbif <- sp.gbif[sp.gbif$establishmentMeans != "INTRODUCED",]

# introduced to Sicily and only domestic on Mallorca

# vittatus
sprp_vittatus <- gBuffer(sprp[sprp$subgroup == "vittatus",], width = 0.2, byid = TRUE)
clip_points <- matrix(c(119.525,-8.598,
                        126.161,-7.765),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)
sprp_vittatus <- raster::bind(sprp_vittatus, clip_points)
sprp_vittatus@data$subgroup <- "vittatus"

# meridionalis
sprp_meridionalis <- gBuffer(sprp[sprp$island %in% c("Sardinia", "Corsica"),], width = 0.1, byid = TRUE)

spain_temp <- terra::buffer(countries50[countries50$name == "Spain",], 0.1)
clip_poly <- matrix(c(-6.2880,36.7198,
                      -4.4866,37.5771,
                      -3.8487,37.8578,
                      -3.1598,38.0415,
                      -3.02319,38.22421,
                      -1.6319,38.7475,
                      -0.7722,37.8905,
                      -0.6954,37.6261,
                      -2.1120,36.5546,
                      -5.6015,35.9596,
                      -6.0499,36.1747,
                      -6.2880,36.7198),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

sprp_meridionalis <- raster::bind(sprp_meridionalis, raster::intersect(spain_temp, clip_poly))
sprp_meridionalis@data$subgroup <- "meridionalis"

# algira
sprp_algira <- gBuffer(sprp[sprp$subgroup == "algira",], width = 0.5, byid = TRUE)
sprp_algira <- sprp_algira - sprp_meridionalis

# indian
sprp_indian <- gBuffer(sprp[sprp$subgroup == "indian",], width = 1, byid = TRUE)
sprp_indian <- sprp_indian - sprp_vittatus

# majori
sprp_majori <- gBuffer(sprp[sprp$subgroup == "majori",], width = 0.1, byid = TRUE)
clip_points <- matrix(c(10.394,42.755,
                        16.436,41.011, 
                        16.592,40.628, 
                        17.032,40.514),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)

sprp_majori <- raster::bind(sprp_majori, clip_points)
sprp_majori@data$subgroup <- "majori"
sprp_majori <- sprp_majori - sprp[sprp$subgroup == "western",]

# western
# large buffer to include recent range expansion into Scotland and because of high density of samples around range borders
sprp_western <- gBuffer(sprp[sprp$subgroup == "western",], width = 6, byid = TRUE)
clip_points <- matrix(c(3.253,39.793,
                        13.466,37.859,
                        18.728,57.589,
                        14.704,37.622),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 80000)
sprp_western <- sprp_western - clip_points - sprp_algira - sprp_indian - sprp[sprp$subgroup == "eastern",] - sprp_majori - sprp_meridionalis

sprp_try <- raster::bind(sprp_algira, sprp[sprp$subgroup == "eastern",], sprp_indian, 
                         sprp_majori, sprp_meridionalis, sprp_vittatus, sprp_western)

# pip_test
buffer_temp <- data.frame(subgroup = c("algira", "eastern", "indian",
                                       "majori", "meridionalis", "vittatus", "western"), 
                          buff = c(0, 0.6, 0, 
                                   0, 0, 0, 0))
sp.gbif <- pip_test(curr.species, dat = sp.gbif, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)


writeOGR(sp.gbif[,names(sp.gbif)[!names(sp.gbif) %in% "gbifID"]], here::here("Data/Big species"), paste(curr.species, "dots"), driver = "ESRI Shapefile")
writeOGR(sprp_try, here::here("Data/Big species"), paste(curr.species, "try"), driver = "ESRI Shapefile")

# Sylvicapra grimmia ####
curr.species <- "Sylvicapra grimmia"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Syncerus caffer ####
curr.species <- "Syncerus caffer"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]

# caffer
sprp_caffer <- sprp[sprp$subgroup == "caffer", ]
sprp_caffer <- gBuffer(sprp_caffer, width = 0.7, byid = TRUE)

clip_points <- matrix(c(30.071875, 0.806031),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 10000)

sprp_caffer <- sprp_caffer - clip_points

sp.gbif_caffer <- sp.gbif[sprp_caffer,]
sp.gbif_caffer@data$subgroup <- "caffer"

# nanus
sprp_nanus <- sprp[sprp$subgroup == "nanus",]
sprp_nanus <- sprp_nanus - sprp_caffer

sprp_nanus <- raster::bind(sprp_nanus, clip_points)
sprp_nanus@data[is.na(sprp_nanus$subgroup), "subgroup"] <- "nanus"

sp.gbif_nanus <- sp.gbif[sprp_nanus,]
sp.gbif_nanus@data$subgroup <- "nanus"

# brachyceros
sprp_brachyceros <- gBuffer(sprp[sprp$subgroup == "brachyceros",], width = 1.5, byid = TRUE)
sprp_brachyceros <- sprp_brachyceros - sprp_nanus

sp.gbif_brachyceros <- sp.gbif[sprp_brachyceros,]
sp.gbif_brachyceros@data$subgroup <- "brachyceros"

# aequinoctialis
# [of ssp. aequinoctialis] This subspecies is sometimes considered to be the same as [brachyceros]
# so am happy to overwrite
# also seems more likely that samples in very close vicinity of nanus sample and close to the border of the nanus range polygon are more likely to be nanus than aequinoctialis because polygons are rough estimates  

sp.gbif_aequinoctialis <- sp.gbif[sprp[sprp$subgroup == "aequinoctialis",],]
sp.gbif_aequinoctialis@data <- sp.gbif_aequinoctialis@data %>%
  mutate(subgroup = case_when(is.na(subgroup) ~ "nanus",
                              subgroup == "brachyceros" ~ "aequinoctialis",
                              TRUE ~ subgroup))

sp.gbif <- raster::bind(sp.gbif_aequinoctialis, sp.gbif_brachyceros, sp.gbif_caffer, sp.gbif_nanus)

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Taxidea taxus ####
curr.species <- "Taxidea taxus"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Tragelaphus angasii ####
curr.species <- "Tragelaphus angasii"
plot_temp(curr.species)
# restricted buffer to avoid samples in private game reserves

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(sprp, 0.1),]
sp.gbif@data$subgroup <- "not used"

# Tragelaphus oryx/Taurotragus oryx ####
curr.species <- "Taurotragus oryx"
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == "Tragelaphus oryx",]) + tm_polygons("subgroup") +
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species,]) + tm_dots()

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == "Tragelaphus oryx",], 2),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Tragelaphus scriptus ####
curr.species <- "Tragelaphus scriptus"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Tragelaphus spekii ####
curr.species <- "Tragelaphus spekii"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Tragelaphus strepsiceros ####
curr.species <- "Tragelaphus strepsiceros"
plot_temp(curr.species)

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("chora", "cottoni", "strepsiceros"), buff = c(0.1, 0, 0.6)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Urocyon cinereoargenteus ####
curr.species <- "Urocyon cinereoargenteus"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# northern
clip_points <- matrix(c(-68.8, 44.4),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 10000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "northern"

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("venezuelae", "northern"), buff = c(0.6, 1.1)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Ursus americanus ####
curr.species <- "Ursus americanus"
plot_temp(curr.species)
# using different method to avoid time taken to buffer mainland polygon

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif@data$subgroup <- "mainland"

# vancouveri
sprp_vancouveri <- terra::buffer(sprp[sprp$subgroup == "vancouveri",], 0.1)
sp.gbif_vancouveri <- sp.gbif[sprp_vancouveri, ]
sp.gbif_vancouveri@data$subgroup <- "vancouveri"
sp.gbif <- sp.gbif[gDisjoint(sprp_vancouveri, sp.gbif, byid = TRUE)[,1],]

# perniger
sprp_perniger <- terra::buffer(sprp[sprp$subgroup == "perniger",], 0.02)
sp.gbif_perniger <- sp.gbif[sprp_perniger, ]
sp.gbif_perniger@data$subgroup <- "perniger"
sp.gbif <- sp.gbif[gDisjoint(sprp_perniger, sp.gbif, byid = TRUE)[,1],]

# hamiltoni = 0
sprp_hamiltoni <- sprp[sprp$subgroup == "hamiltoni",]
sp.gbif_hamiltoni <- sp.gbif[sprp_hamiltoni,]
sp.gbif_hamiltoni@data$subgroup <- "hamiltoni"
sp.gbif <- sp.gbif[gDisjoint(sprp_hamiltoni, sp.gbif, byid = TRUE)[,1],]

# carlottae pugnax
sprp_carlottae <- terra::buffer(sprp[sprp$subgroup == "carlottae pugnax",], 0.1)
sp.gbif_carlottae <- sp.gbif[sprp_carlottae,]
sp.gbif_carlottae@data$subgroup <- "carlottae pugnax"
sp.gbif <- sp.gbif[gDisjoint(sprp_carlottae, sp.gbif, byid = TRUE)[,1],]

sp.gbif <- raster::bind(sp.gbif, sp.gbif_carlottae, sp.gbif_hamiltoni, sp.gbif_perniger, sp.gbif_vancouveri)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Ursus arctos ####
curr.species <- "Ursus arctos"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# arctos = 0
sprp_arctos <- gBuffer(sprp[sprp$subgroup == "arctos",], witdh = 0.6, byid = TRUE)
sprp_arctos <- sprp_arctos - sprp[sprp$subgroup == "collaris beringianus lasiotus",]

clip_points <- matrix(c(53.059303, 68.059823,
                        28.521339, 44.403096),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_arctos <- raster::bind(sprp_arctos, terra::buffer(clip_points, 100000))
sprp_arctos@data[is.na(sprp_arctos$subgroup), "subgroup"] <- "arctos"

sprp <- raster::bind(sprp[sprp$subgroup != "arctos",], sprp_arctos)

# collaris beringianus lasiotus = 0
clip_points <- matrix(c(109.922607, 51.681679,
                        127.667930, 37.551208,
                        137.436492, 54.868479),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 100000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "collaris beringianus lasiotus"

# crowtheri = NA
# horribilis = 0

sprp_horribilis <- gBuffer(sprp[sprp$subgroup == "horribilis",], width = 0.6, byid = TRUE)
sprp_horribilis <- sprp_horribilis - terra::buffer(sprp[sprp$subgroup == "sitkensis",], 0.1)
sprp <- raster::bind(sprp[sprp$subgroup != "horribilis",], sprp_horribilis)

# marsicanus = 0.1
# middendorfi = 0
# pruinosus isabellinus gobiensis = 0.1
# pyrenaicus = 0.5
# sitkensis = 0.1
# syriacus = 0.1
# syriacus samples sometimes given arctos epithet "The Syrian brown bear (Ursus arctos syriacus or Ursus arctos arctos)" (wiki)
# so am happy to overwrite arctos

# ungavaensis = NA

# pip_test
buffer_temp <- data.frame(subgroup = c("arctos", "collaris beringianus lasiotus", "horribilis", "marsicanus", "middendorfi", "pruinosis isabellinus gobiensis", "pyrenaicus", "sitkensis", "syriacus"), 
                          buff = c(0, 0, 0, 0.1, 0, 0.1, 0.5, 0.1, 0.1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Ursus maritimus ####
curr.species <- "Ursus maritimus"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

# Vulpes lagopus ####
curr.species <- "Vulpes lagopus"
plot_temp(curr.species)

# lagopus
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
clip_points <- matrix(c(-166.67, 53.85,
                        -145.75, 60.5,
                        15.79, 69,
                        178.84, 71.12),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 10000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "lagopus"

# pip_test
buffer_temp <- data.frame(subgroup = c("beringensis", "foragoapusis", "fuliginosus", "lagopus", "pribilofensis"), 
                          buff = c(0.1, 0.5, 0.5, 1, 0.5))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Vulpes macrotis ####
curr.species <- "Vulpes macrotis"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
# macrotis
poly_temp <- sprp[sprp$subgroup == "macrotis",]
poly_temp <- gBuffer(poly_temp, width = 1.5, byid = TRUE)
poly_temp <- poly_temp - sprp[sprp$subgroup == "mutica",]

sprp_try <- raster::bind(sprp[sprp$subgroup != "macrotis",], poly_temp)

# stray velox sample

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species | IUCN_Native_Data$binomial == "Vulpes velox",]) + tm_polygons() + 
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species | GBIF_Spatial$species == "Vulpes velox",]) + tm_dots("verbatimScientificName")

clip_points <- matrix(c(-104.08, 40.99),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_try <- raster::bind(sprp_try, terra::buffer(clip_points, 50000))
sprp_try@data[is.na(sprp_try$subgroup), "subgroup"] <- "velox"

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = data.frame(subgroup = c("macrotis", "mutica", "velox"), 
                                      buff = c(0, 0, 0)))
sp.gbif@data[sp.gbif$subgroup == "velox", "species"] <- "Vulpes velox"

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Vulpes velox ####
curr.species <- "Vulpes velox"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif.raw <- few.raw.spatial[few.raw.spatial$species == curr.species,]
sp.gbif.raw <- raster::bind(sp.gbif.raw, sp.gbif[sp.gbif$verbatimScientificName == "Vulpes macrotis",])

tm_shape(sprp) + tm_polygons("subgroup") +
  tm_shape(sp.gbif.raw) + tm_dots("verbatimScientificName") 

# Vulpes velox
clip_points <- matrix(c(-94.607,39.036,
                        -114.094,51.117,
                        -115.220,55.076,
                        -101.305,32.283),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(gBuffer(sprp, width = 0.1, byid = TRUE), terra::buffer(clip_points, 80000))
sprp@data$subgroup <- "velox"
sp.gbif <- sp.gbif.raw[sprp,]
sp.gbif@data$subgroup <- "not used"

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Vulpes macrotis
sprp_macrotis <- IUCN_Native_Data[IUCN_Native_Data$binomial == "Vulpes macrotis",]
clip_points <- matrix(c(-107.0814,25.0714,
                        -97.239,15.905,
                        -93.746,30.434),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 80000)
sprp_macrotis <- raster::bind(gBuffer(sprp_macrotis, width = 1.5, byid = TRUE), clip_points)
sprp_macrotis <- sprp_macrotis - sprp
sprp_macrotis@data$subgroup <- "macrotis"
sp.macrotis <- sp.gbif.raw[sprp_macrotis,]
sp.macrotis@data$species<- "Vulpes macrotis"

# macrotis mutica
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.macrotis) + tm_dots()

sprp_try <- raster::bind(sprp_try, clip_points)
sprp_try@data[is.na(sprp_try$subgroup), "subgroup"] <- "macrotis"

# pip_test
sp.macrotis <- pip_test("Vulpes macrotis", dat = sp.macrotis, 
                    range.polygon = sprp_try, 
                    buff = data.frame(subgroup = c("macrotis", "mutica"), buff = c(0, 0)))

GBIF_Spatial <- raster::bind(GBIF_Spatial, sp.macrotis)

# Vulpes vulpes ####
curr.species <- "Vulpes vulpes" # crashes R
sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

sp.gbif <- sp.gbif[terra::buffer(sprp, 5),]

# schrencki
sprp_schrencki <- gBuffer(sprp[sprp$subgroup == "schrencki",], width = 0.1, byid = TRUE)
clip_poly <- matrix(c(141.1, 54.3,
                      144.4, 54.3,
                      144.4, 45.6,
                      141.1, 45.6,
                      141.1, 54.3),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_schrencki <- raster::bind(sprp_schrencki, clip_poly)
sprp_schrencki@data$subgroup <- "schrencki"

# japonica 
sprp_japonica <- gBuffer(sprp[sprp$subgroup == "japonica",], width = 0.1, byid = TRUE)
clip_points <- matrix(c(132.660,33.594),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp_japonica <- raster::bind(sprp_japonica, terra::buffer(clip_points, 50000))
sprp_japonica@data$subgroup <- "japonica"

# niloticus 
sprp_niloticus <- gBuffer(sprp[sprp$subgroup == "niloticus",], width = 1.5, byid = TRUE)

# barbara atlantica
sprp_barbara <- gBuffer(sprp[sprp$subgroup == "barbara atlantica",], width = 0.5, byid = TRUE)
clip_points <- matrix(c(132.660,33.594),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)
sprp_barbara <- sprp_barbara - clip_points - spain_temp

clip_points <- matrix(c(10.825,33.803),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)
sprp_barbara <- raster::bind(sprp_barbara, clip_points)
sprp_barbara@data$subgroup <- "barbara atlantica"

# silacea
sprp_silacea <- gBuffer(sprp[sprp$subgroup == "silacea",], width = 0.1, byid = TRUE)
sprp_silacea <- sprp_silacea - sprp[sprp$subgroup == "crucifera",]

# vulpes
sprp_vulpes <- gBuffer(sprp[sprp$subgroup == "vulpes",], width = 0.05, byid = TRUE)
clip_points <- matrix(c(16.532,56.505,
                        16.870,56.994,
                        18.325,57.208,
                        18.664,57.722,
                        18.488,59.289,
                        22.414,60.180,
                        49.649,69.289,
                        23.7344,70.7012,
                        18.828,69.735,
                        18.072,69.471,
                        17.447,69.314,
                        16.414,68.750,
                        15.958,69.122,
                        15.586,68.546,
                        14.943,68.570,
                        14.174,68.186,
                        13.430,68.078,
                        12.493,65.964,
                        11.148,64.895,
                        9.0581,63.5176,
                        8.3504,63.2105,
                        7.6531,63.0648,
                        6.3209,62.4456,
                        5.5871,60.5669,
                        5.2488,59.2243,
                        11.6132,58.0899),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)
sprp_vulpes <- raster::bind(sprp_vulpes, clip_points)
sprp_vulpes@data$subgroup <- "vulpes"
sprp_vulpes <- sprp_vulpes - sprp[sprp$subgroup %in% c("crucifera", "unknown"),]

# unknown
sprp_unknown <- gBuffer(sprp[sprp$subgroup == "unknown",], width = 0.4, byid = TRUE)
clip_points <- matrix(c(32.8807,34.8661),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 60000)
sprp_unknown <- raster::bind(sprp_unknown, clip_points)
sprp_unknown@data$subgroup <- "unknown"
sprp_unknown <- sprp_unknown - sprp_niloticus - sprp[sprp$subgroup %in% c("crucifera", "vulpes"),]

# crucifera
sprp_crucifera <- gBuffer(sprp[sprp$subgroup == "crucifera",], width = 0.5, byid = TRUE)
clip_points <- matrix(c(-4.627,54.184,
                        -6.568,57.435,
                        10.5307,55.2159,
                        10.9414,54.9260,
                        11.5935,54.7570,
                        12.3585,54.9824,
                        11.9076,54.9905,
                        11.5171,55.4977,
                        12.1773,55.7795,
                        22.4904,58.4443,
                        22.6140,58.8556),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points <- terra::buffer(clip_points, 50000)
sprp_crucifera <- raster::bind(sprp_crucifera, clip_points)
sprp_crucifera@data$subgroup <- "crucifera"
sprp_crucifera <- sprp_crucifera - sprp_unknown - sprp_silacea - sprp_vulpes

sprp_try <- raster::bind(sprp_barbara, sprp_crucifera, sprp_japonica, sprp_niloticus,
                         sprp_schrencki, sprp_silacea, sprp_unknown, sprp_vulpes, 
                         sprp[sprp$subgroup %in% c("ichnusae", "splendidissima"),])

# pip_test
buffer_temp <- data.frame(subgroup = c("barbara atlantica", "crucifera", "ichnusae", 
                                       "japonica", "niloticus", "schrencki",
                                       "silacea", "splendidissima", "unknown", "vulpes"), 
                          buff = c(0, 0, 0.1, 
                                   0, 0, 0,
                                   0, 0.1, 0, 0))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp_try, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Finishing up ####

GBIF_Subgroups <- as.data.frame(GBIF_Spatial)

# tidying
rm(list = ls(pattern = "^sprp"))
rm(list = ls(pattern = "_temp$"))
rm(turkey_rangepol, sp.gmpd.points, list = ls(pattern = "^clip_"))
rm(curr.species, list = ls(pattern = "^sp."), 
   list = ls(pattern = "^gbif_"), list = ls(pattern = "^few."))

