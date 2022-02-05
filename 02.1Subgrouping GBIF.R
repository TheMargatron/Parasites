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

## Canis aureus ####
curr.species <- "Canis aureus"
plot_temp(curr.species)

## Canis latrans ####
curr.species <- "Canis latrans"
# map # plot_temp(curr.species)
# unsure of island subspecies
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

## Canis lupus ####
curr.species <- "Canis lupus"
plot_temp(curr.species)

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

## Capreolus capreolus ####
curr.species <- "Capreolus capreolus"
gbif_sp_temp <- GBIF_Spatial[GBIF_Spatial$species == curr.species & !is.na(GBIF_Spatial$subgroup),]
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(gbif_sp_temp) + tm_dots("subgroup")

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
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,],]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Giraffa camelopardalis ####
curr.species <- "Giraffa camelopardalis"
plot_temp(curr.species)

# giraffa = 0.5

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

## Kobus ellipsiprymnus ####
curr.species <- "Kobus ellipsiprymnus"
plot_temp(curr.species)

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
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

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

## Lynx pardinus ####
curr.species <- "Lynx pardinus"
plot_temp(curr.species)

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

## Martes foina ####
curr.species <- "Martes foina"
plot_temp(curr.species)

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

## Meles meles ####
curr.species <- "Meles meles"
plot_temp(curr.species) # lots of data, crashes R

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

## Mustela erminea ####
curr.species <- "Mustela erminea"
plot_temp(curr.species) # lots of samples, crashes R

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

buffer_temp <- data.frame(subgroup = c("american", "numidica", "formosana", "eurasian"), 
                          buff = c(0, 0, 0.3, 0.3))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Mustela putorius ####
curr.species <- "Mustela putorius"
plot_temp(curr.species)

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
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("leucogaster", "bezoarticus", "arerunguaeansis"), buff = c(2.5, 2, 0)))

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
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

## Procyon lotor ####
curr.species <- "Procyon lotor"
plot_temp(curr.species)

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
plot_temp(curr.species)

sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("chanleri", "fulvorufula"), buff = c(1.5, 0)))
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Rupicapra pyrenaica ####
curr.species <- "Rupicapra pyrenaica"
plot_temp(curr.species)

buffer_temp <- data.frame(subgroup = c("ornata", "parva", "pyrenaica"), buff = c(0.5, 0.5, 1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Rupicapra rupicapra ####
# Need to rerun 01 scripts so that Aubrac is kept 
curr.species <- "Rupicapra rupicapra"
plot_temp(curr.species)

sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]

# pip_test
buffer_temp <- data.frame(subgroup = c("asiatica", "balcanica", "carpatica", "caucasica", "rupicapra cartusiana", "tatrica"), 
                          buff = c(0.1, 0.7, 0.7, 0.1, x, 0.1))
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = buffer_temp)
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Sus scrofa ####
curr.species <- "Sus scrofa"
plot_temp(curr.species)

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

sprp <- raster::bind(sprp[sprp$subgroup != "macrotis",], poly_temp)

# stray velox sample

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species | IUCN_Native_Data$binomial == "Vulpes velox",]) + tm_polygons() + 
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species | GBIF_Spatial$species == "Vulpes velox",]) + tm_dots("verbatimScientificName")

clip_points <- matrix(c(-104.08, 40.99),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS(proj4string(IUCN_Native_Data)))
sprp <- raster::bind(sprp, terra::buffer(clip_points, 50000))
sprp@data[is.na(sprp$subgroup), "subgroup"] <- "velox"

# pip_test
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = sprp, 
                    buff = data.frame(subgroup = c("macrotis", "mutica", "velox"), 
                                      buff = c(0, 0, 0)))
sp.gbif@data[sp.gbif$subgroup == "velox", "species"] <- "Vulpes velox"

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

# Vulpes velox ####
curr.species <- "Vulpes velox"
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons() + 
  tm_shape(GBIF_Spatial[GBIF_Spatial$species == curr.species,]) + tm_dots("verbatimScientificName")
# looks likes some are actually macrotis

GBIF_Spatial@data <- GBIF_Spatial@data %>%
  mutate(species = case_when(verbatimScientificName == "Vulpes velox subsp. macrotis" ~ "Vulpes macrotis",
                             TRUE ~ species)) %>%
  mutate(subgroup = case_when(verbatimScientificName == "Vulpes velox subsp. macrotis" ~ "macrotis",
                             TRUE ~ subgroup)) %>%
  mutate(subgroup = case_when(species == curr.species ~ "not used",
                              TRUE ~ subgroup))

## Vulpes vulpes ####
curr.species <- "Vulpes vulpes" # crashes R
plot_temp(curr.species)
