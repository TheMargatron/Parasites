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

## Aepyceros melampus ####
curr.species <- "Aepyceros melampus"
plot_temp(curr.species)

# sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
#                     range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
#                     buff = data.frame(subgroup = c("hecki soemmeringii", "jubatus"), buff = c(1.8, 1)))
# 
# GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$HostCorrectedName != curr.species,], sp.gbif)

## Alcelaphus buselaphus ####
curr.species <- "Alcelaphus buselaphus"
plot_temp(curr.species)

## Alces alces ####
# buffer = 8
curr.species <- "Alces alces"
plot_temp(curr.species)

## Antidorcas marsupialis ####
curr.species <- "Antidorcas marsupialis"
plot_temp(curr.species)

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

# Canis latrans ####
curr.species <- "Canis latrans"
# map # plot_temp(curr.species)
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

## Capra pyrenaica ####
curr.species <- "Capra pyrenaica"
plot_temp(curr.species)

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

## Connochaetes taurinus ####
curr.species <- "Connochaetes taurinus"
plot_temp(curr.species)

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

## Damaliscus lunatus ####
curr.species <- "Damaliscus lunatus"
plot_temp(curr.species)

## Damaliscus pygargus ####
curr.species <- "Damaliscus pygargus"
plot_temp(curr.species)

# Equus grevyi ####
curr.species <- "Equus grevyi"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Equus quagga ####
curr.species <- "Equus quagga"
plot_temp(curr.species)

## Equus zebra ####
curr.species <- "Equus zebra"
plot_temp(curr.species)
sp.gbif <- pip_test(curr.species, dat = GBIF_Spatial, 
                    range.polygon = IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 
                    buff = data.frame(subgroup = c("zebra", "hartmannae"), buff = c(1, 1)))

GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Felis silvestris ####
curr.species <- "Felis silvestris"
plot_temp(curr.species)

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

## Hippotragus niger ####
curr.species <- "Hippotragus niger"
plot_temp(curr.species)

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

## Leopardus pardalis ####
curr.species <- "Leopardus pardalis"
plot_temp(curr.species)

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

## Lynx lynx ####
curr.species <- "Lynx lynx"
plot_temp(curr.species)

## Lynx pardinus ####
curr.species <- "Lynx pardinus"
plot_temp(curr.species)

## Lynx rufus ####
curr.species <- "Lynx rufus"
plot_temp(curr.species)

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

# Mustela erminea ####
curr.species <- "Mustela erminea"
plot_temp(curr.species) # lots of samples, crashes R

# Mustela lutreola ####
curr.species <- "Mustela lutreola"
plot_temp(curr.species)

sp.gbif <- GBIF_Spatial[GBIF_Spatial$species == curr.species,]
sp.gbif <- sp.gbif[terra::buffer(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,], 1),]
sp.gbif@data$subgroup <- "not used"
GBIF_Spatial <- raster::bind(GBIF_Spatial[GBIF_Spatial$species != curr.species,], sp.gbif)

## Mustela nivalis ####
curr.species <- "Mustela nivalis"
plot_temp(curr.species)

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

## Odocoileus virginianus ####
curr.species <- "Odocoileus virginianus"
plot_temp(curr.species)

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

## Ovibos moschatus ####
curr.species <- "Ovibos moschatus"
plot_temp(curr.species)

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

## Rangifer tarandus ####
curr.species <- "Rangifer tarandus"
plot_temp(curr.species)

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
curr.species <- "Rupicapra rupicapra"
plot_temp(curr.species)

## Sus scrofa ####
curr.species <- "Sus scrofa"
plot_temp(curr.species)

# Sylvicapra grimmia ####
curr.species <- "Sylvicapra grimmia"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

## Syncerus caffer ####
curr.species <- "Syncerus caffer"
plot_temp(curr.species)

# Taxidea taxus ####
curr.species <- "Taxidea taxus"
plot_temp(curr.species)
GBIF_Spatial@data[GBIF_Spatial$species == curr.species, "subgroup"] <- "not used"

## Tragelaphus angasii ####
curr.species <- "Tragelaphus angasii"
plot_temp(curr.species)

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

## Tragelaphus strepsiceros ####
curr.species <- "Tragelaphus strepsiceros"
plot_temp(curr.species)

## Urocyon cinereoargenteus ####
curr.species <- "Urocyon cinereoargenteus"
plot_temp(curr.species)

## Ursus americanus ####
curr.species <- "Ursus americanus"
plot_temp(curr.species)

# northern
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]
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

## Ursus arctos ####
curr.species <- "Ursus arctos"
plot_temp(curr.species)

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
