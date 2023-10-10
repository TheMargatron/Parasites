# Label IUCN and GMPD data with subgroups
# Written by Margaret Bolton mb804(at)exeter.ac.uk


# Prep ####

# going by subgroup rather than subspecies because of geographic groupings
IUCN_Native_Data$subgroup <- IUCN_Native_Data$subspecies

IUCN_Native_Data <- IUCN_Native_Data %>%
  dplyr::mutate(subgroup = case_when(stringr::str_detect(subgroup, " ssp. ") ~ stringr::str_split(subgroup, "ssp. ", simplify = TRUE)[,2],
                                     TRUE ~ subgroup))

# ignoring invalid geometries introduced by S2 for now 
# currently (27/09/23) st_make_valid does nothing
sf::sf_use_s2(FALSE)

GMPD_Data$subgroup <- GMPD_Data$HostReportedSubspecies

GMPD_Spatial <- sf::st_as_sf(GMPD_Data,
                             coords = c("Longitude", "Latitude"),
                             crs = Projection_String)

# subsp_list <- GMPD_Data %>% 
#   dplyr::select(HostCorrectedName) %>%
#   dplyr::distinct() %>% 
#   dplyr::arrange(HostCorrectedName)
# 
# 
# write.csv(subsp_list, here::here("Data/Extras/Subs_temp.csv"), row.names = FALSE)

## Acinonyx jubatus ####
# There are multiple subspecies that are quite geographically separated, 
# GMPD data only occupies the range of A. j. jubatus
# I've removed polygons based on wiki distribution map
# https://en.wikipedia.org/wiki/Cheetah#/media/File:Acinonyx_jubatus_subspecies_range_IUCN_2015.png

species_temp <- "Acinonyx jubatus"
range_temp <- IUCN_Native_Valid[IUCN_Native_Valid$binomial == species_temp, ]

#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# upper
clip_poly_temp <- matrix(c(33.0, 0,
                           -10, 0,
                           -10, 40,
                           60.0, 40,
                           60.0, 3.7,
                           45.0, 3.7,
                           35.5, 5,
                           33.0, 0), 
                         ncol = 2, byrow = TRUE) %>%
  list() %>% 
  st_polygon() %>% 
  st_sfc(crs = Projection_String) %>% 
  st_sf(crs = Projection_String)

clip_points_temp <- matrix(c(38.7, 4.2),
                           ncol = 2, byrow = TRUE) %>% 
  as.data.frame() %>% 
  st_as_sf(crs = Projection_String,
    coords = c(1,2))

clip_poly_test <- range_temp %>% 
  st_intersection(clip_points_temp) 

clip_poly_temp <- gUnion(clip_poly_temp, range_temp[clip_points_temp, ])

range_temp_jubatus <- range_temp - clip_poly_temp
range_temp_jubatus@data$subgroup <- "jubatus"
range_temp_otherssp <- raster::intersect(range_temp, clip_poly_temp)
range_temp_try <- raster::bind(range_temp_jubatus, range_temp_otherssp)
range_temp_try@data[is.na(range_temp_try@data$subgroup), "subgroup"] <- "hecki soemmeringii"

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# No gmpd subgroup info to begin with so am happy overwriting
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "jubatus"

rm(list = ls(pattern = "^range_temp"))

## Aepyceros melampus ####
# Two subspecies, both represented in GMPD 
# Also black-faced impala seems geographically distinct 
# https://en.wikipedia.org/wiki/Impala#/media/File:Aepyceros_melampus.svg
# IUCN: "In Namibia, the Black-faced Impala is naturally confined to the Kaokoland in the north-west, and neighbouring south-western Angola"

species_temp <- "Aepyceros melampus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

# All petersi samples already have subgroup data
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp & 
                    is.na(GMPD_Spatial@data$subgroup), "subgroup"] <- "melampus"

## Alcelaphus buselaphus ####
# 8 subspecies, GMPD data appears to represent major and cokii, and they are geographically distinct
# https://en.wikipedia.org/wiki/Hartebeest#/media/File:Alcelaphus_recent.png
# There is only one sample location for each subspecies 

species_temp <- "Alcelaphus buselaphus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

# All cokii samples already have subgroup data
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp & 
                    is.na(GMPD_Spatial@data$subgroup), "subgroup"] <- "major"

## Alces alces ####
# There are multiple subspecies, GMPD represents shirasi, gigas, andersoni, and americana in North America
# and alces in Europe/Western Russia, but not buturlini, cameloides, or pfizenmayeri (east of Yenisei river)
# Splitting European polygon around Yenisei should do it. 
# Yenisei data source: https://doi.org/10.1016/j.dib.2018.09.016

species_temp <- "Alces alces"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

# Used combination of Yenisei and Angara 
YA_temp <- River_Data50[River_Data50$name %in% c("Yenisey", "Angara"),]

YA_temp <- disaggregate(YA_temp)
YA_temp$ID <- LETTERS[1:nrow(YA_temp)]
#map# tm_shape(YA_temp) + tm_lines("ID", lwd = 2)
drop_ID <- "A|G|I|K"
#map# tm_shape(YA_temp[str_detect(YA_temp$ID, drop_ID, negate = TRUE),]) + tm_lines("ID", lwd = 2)

YA_temp <- YA_temp[-grep(drop_ID, YA_temp$ID),] # to avoid self intersections
YA_coords <- unlist(coordinates(YA_temp), recursive = FALSE)
names(YA_coords) <- YA_temp$ID
YA_coords[["C"]] <- YA_coords[["C"]][nrow(YA_coords[["C"]]) - 1:1,]
YA_coords[["L"]] <- YA_coords[["L"]][nrow(YA_coords[["L"]]):1,]
YA_coords <- YA_coords[c("D","B","E","C","F","H","J","L")]
YA_coords <- do.call(rbind, YA_coords)

YA_curr_coords <- rbind(YA_coords,
                   matrix(c(YA_coords[nrow(YA_coords), 1], (range_temp@bbox[2,2] + 1),
                            (range_temp@bbox[1,2] + 1), (range_temp@bbox[2,2] + 1),
                            (range_temp@bbox[1,2] + 1), (range_temp@bbox[2,1] - 1),
                            YA_coords[1,1], (range_temp@bbox[2,1] - 1),
                            YA_coords[1,1], YA_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(YA_curr_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

#map# tm_shape(clip_poly_temp) + tm_polygons() + tm_shape(range_temp) + tm_polygons(alpha = 0)

# alces subspecies
range_temp_try <- range_temp - clip_poly_temp

poly_temp <- Polygon(matrix(c(0, 75,
                          110, 75,
                          110, 30,
                          0, 30,
                          0, 75), 
                        ncol = 2, byrow = TRUE))
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_alces <- raster::intersect(range_temp_try, poly_temp)
range_temp_alces@data$subgroup <- "alces"

# americanus subspecies group (americana, andersoni, gigas, shirasi)
range_temp_americanus <- range_temp_try - poly_temp
range_temp_americanus@data$subgroup <- "americana andersoni gigas shirasi"

# russian/asian subspecies (buturlini, cameloides, pfizenmayeri)
range_temp_otherssp <- raster::intersect(range_temp, clip_poly_temp)
range_temp_otherssp@data$subgroup <- "buturlini cameloides pfizenmayeri"

range_temp_try <- raster::bind(range_temp_alces, range_temp_americanus, range_temp_otherssp)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("alces", "americana andersoni gigas shirasi"), buff = c(0.5, 4.5)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Antidorcas marsupilis ####
# There are three recognised subspecies (wiki) 
# I only seem to have A. m. marsupialis in GMPD data whose range lies south of Orange river (wiki) 
# In the Eastern part of this range it appears to be restricted by the Vaal rather than the Orange as wiki describes it as extending up to Kimberley which is North of Orange
# "trekbokken" historically occurred across the Orange river
# https://www.ewt.org.za/wp-content/uploads/2019/02/3.-Springbok-Antidorcas-marsupialis_LC.pdf
# However, such large scale migrations no longer occur 
"Anderson C, Schultze E, Codron D, Bissett C, Gaylard A, Child MF. 2016. A conservation
assessment of Antidorcas marsupialis. In Child MF, Roxburgh L, Do Linh San E, Raimondo D, Davies-Mostert HT, editors.
The Red List of Mammals of South Africa, Swaziland and Lesotho. South African National Biodiversity Institute and
Endangered Wildlife Trust, South Africa."
# 
"While the subspecies distinction is debated, morphometric data reveal a difference in size between 
Springbok occurring on either sides of the Orange and Vaal Rivers, which is evidence for maintaining 
subspecies status (Peters & Brink 1992). "

species_temp <- "Antidorcas marsupialis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

OV_temp <- River_Data50[(River_Data50$name == "Orange"| River_Data50$name == "Vaal"),]

#map# tm_shape(range_temp) + tm_polygons() + tm_shape(OV_temp) + tm_lines() + tm_shape(sp.gmpd.points) + tm_dots()

OV_temp <- disaggregate(OV_temp)
OV_temp$ID <- LETTERS[1:nrow(OV_temp)]
#map# tm_shape(OV_temp) + tm_lines("ID", lwd = 2)
drop_ID <- "G|A|E|B|D" ##### added D
#map# tm_shape(OV_temp[str_detect(OV_temp$ID, drop_ID, negate = TRUE),]) + tm_lines("ID", lwd = 2)
OV_temp <- OV_temp[-grep(drop_ID, OV_temp$ID),] # to avoid self intersections

OV_edit <- coordinates(OV_temp[OV_temp$ID == "F",])[[1]][[1]]
split_point <- nearestPointOnLine(OV_edit, 
                                  tail(coordinates(OV_temp[OV_temp$ID == "H",])[[1]][[1]], 1L))

split_row <- which(apply(OV_edit, 1, function(x) all(x == split_point)))
OV_edit <- OV_edit[split_row:nrow(OV_edit),]
#map# tm_shape(OV_temp) +tm_lines() + tm_shape(SpatialLines(list(Lines(Line(OV_edit), ID = "6"))))  + tm_lines( col = "red")

OV_coords <- unlist(coordinates(OV_temp), recursive = FALSE)
names(OV_coords) <- OV_temp$ID
OV_coords[["F"]] <- OV_edit

# Turning them round
OV_coords <- OV_coords[c("J","I","C","H","F")] #### removed D between J and I
OV_coords <- do.call(rbind, OV_coords)

OV_coords <- rbind(OV_coords,
                   matrix(c((range_temp@bbox[1,1] - 1), OV_coords[nrow(OV_coords), 2],
                            (range_temp@bbox[1,1] - 1), (range_temp@bbox[2,2] + 1),
                            (range_temp@bbox[1,2] + 1), (range_temp@bbox[2,2] + 1),
                            (range_temp@bbox[1,2] + 1), OV_coords[1,2], 
                            OV_coords[1,1], OV_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(OV_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

#map# tm_shape(clip_poly_temp) + tm_polygons()

range_temp_marsupialis <- range_temp - clip_poly_temp
range_temp_marsupialis@data$subgroup <- "marsupialis"
range_temp_otherssp <- raster::intersect(range_temp, clip_poly_temp)
range_temp_try <-  raster::bind(range_temp_marsupialis, range_temp_otherssp)

range_temp_try@data[is.na(range_temp_try@data$subgroup), "subgroup"] <- "angolensis hofmeyri"
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "marsupialis"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "^OV_"))

## Antilocapra americana ####
# three subspecies, but not all represented by GMPD (missing peninsularis)
# "Mitochondrial DNA analyses since the early 1990s support the idea of clines within a wide-ranging species rather than separate subspecies (O’Gara and Yoakum 2004)."
# but peninsularis is quite geographically distinct from other subpopulations
species_temp <- "Antilocapra americana"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

range_temp_try <- raster::disaggregate(range_temp)

range_temp_try@data$subgroup <- "americana sonoriensis"
range_temp_try@data[17,"subgroup"] <- "peninsularis"

range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "americana sonoriensis"

rm(list = ls(pattern = "^range_temp"))

## Axis axis ####
# Only have one sample point :(
species_temp <- "Axis axis"
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]

#map# unique(sp.gmpd.points@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Bison bison ####
# "There are two recognized subspecies in North America: Bison bison bison and B. b. athabascae."

"A herd of hybrid plains bison x wood bison lived wild in the Yukon, Canada. The Wood bison is a 
distinct subspecies that almost became extinct in the 20th century. In an attempt to save the Plains 
bison subspecies, between 1925 and 1928, thousands of Plains bison were released into Wood Buffalo 
Park (a preserve for the Wood bison subspecies). They readily interbred and produced a 12,000 strong 
herd by 1934. The Wood bison was nearly hybridized into extinction. A small genetically pure herd was 
recovered from an isolated area in 1959 and is now being kept isolated from introduced Plains bison. 
Unfortunately, recent genetic testing seems to indicate that these supposedly pure wood bison are 
also hybridized with the plains subspecies, though the majority of the genetic makeup is that of the 
'Wood Buffalo.'" #https://en.wikipedia.org/wiki/Bovid_hybrid

species_temp <- "Bison bison"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

IUCN_Native_Data@data[IUCN_Native_Data$binomial == species_temp, "subgroup"] <- "bison athabascae"
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "bison athabascae"

## Bison bonasus ####
# "Three subspecies of European bison existed in the recent past, but only one, 
# the nominate subspecies (B. b. bonasus), survives today" (wiki)
# two subspecies are widely recognized as the Lowland Bison (Bison bonasus bonasus) and 
# the Caucasian Bison (Bison bonasus caucasicus) (Kowalczyk and Plumb 2020)" (IUCN)

# "European bison herds, scattered across Central and Eastern Europe, represent two genetic lines"
# the lowland line (Poland, Belarus, and Lithuania) and the lowland-Caucasian line (southern Poland, Russia, Ukraine and Slovakia).
# https://animaldiversity.org/accounts/Bison_bonasus/

# The distribution of the lineages is artificial because it's determined by reintroduction programmes, rather than natural dispersal
# That means that the different habitat patches they occupy will not relate to their lineage, just to their overall species distribution

species_temp <- "Bison bonasus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

IUCN_Native_Data@data[IUCN_Native_Data$binomial == species_temp, "subgroup"] <- "bonasus caucasicus"
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "bonasus caucasicus"

## Blastocerus dichotomus ####
# No taxonomic notes on IUCN, no description of subspecies on wiki or msotw
#map# species_temp <- "Blastocerus dichotomus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Canis adustus ####
# There are seven recognized subspecies of the side-striped jackal:[2]

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
# Also range is continuous across Africa with no clear subpopulations

# I think the species where I can't distinguish between subspecies will still be informative for latitudinal 
# gradient, but less informative for within-range analysis.
# Might be interesting to compare result from uncleaned polygons with cleaned. 
# Split the streams again after adding this step?

#map# species_temp <- "Canis adustus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Canis aureus ####
"Recent studies based on mtDNA and morphology have shown that 'Golden Jackals' in Africa are larger in size than those from 
Eurasia and are actually more closely related to the Grey Wolf Canis lupus. African animals hence represent a previously 
overlooked distinct species, the African Wolf, Canis lupaster (see Rueness et al. 2011, Gaubert et al. 2012, Koepfli et al. 
2015, Viranta et al. 2017). However, the putative presence of Golden Jackal in the Sinai Peninsula of Egypt remains unclear 
(see Gaubert et al. 2012, Viranta et al. 2017)." #IUCN

# wiki s 7 subspecies
# aureus: Middle East, Iran, Turkmenistan, Afghanistan, Pakistan and Western India
# cruesemanni: Thailand
# ecsedensis: Pannonian Basin, Central Europe
# indicus: India, Nepal, Bangladesh, Bhutan
# moreotica: Southeastern Europe, Moldova, Asia Minor and the Caucasus
# naria: Coastal South West India, Sri Lanka
# syriacus: Israel, Syria,[38] Lebanon,[62] and Jordan

# Have aureus in Iran and probably moreotica in Greece
# Could subdivide according to: https://doi.org/10.1093/mspecies/sey002
# or: https://en.wikipedia.org/wiki/Sri_Lankan_jackal#/media/File:Canis_aureus_subspecies_range.png
# However: 
"In Greece and Dalmatia, C. aureus is documented from the Holocene (Sommer and Benecke 2005; 
Malez 1984 in Rutkowski et al. 2015) and the ancient Mediterranean populations have persisted 
and merged with jackals coming from Asia (Fabbri et al. 2014; Rutkowski et al. 2015)."
# Indicates good dispersal ability and admixture of subspecies populations 
# little reason to expect they are genetically isolated, but they will be locally adapted to some degree
species_temp <- "Canis aureus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp,]

# indicus naria and cruesmanni
clip_poly_temp <- countries50[countries50$name %in% c("India", "Nepal", "Bhutan", 
                                                 "Bangladesh", "Myanmar", "Lao PDR",
                                                 "Thailand", "Cambodia", "Vietnam",
                                                 "Sri Lanka"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.145)

poly_temp <- Polygon(matrix(c(68.5, 24.6,
                              70.7, 24.6,
                              70.72, 24.48,
                              70.59, 24.43,
                              69.73, 23.85,
                              68.5, 23.85,
                              68.5, 24.6), 
                            ncol = 2, byrow = TRUE))
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

clip_poly_temp <- clip_poly_temp - poly_temp

range_temp_indicus <- raster::intersect(range_temp, clip_poly_temp)
range_temp_indicus@data$subgroup <- "indicus naria cruesmanni"

# moreotica and ecsedensis
poly_temp <- countries50[countries50$name %in% c("Turkey", "Greece", "Bulgaria", 
                                                 "Albania", "Macedonia", "Romania",
                                                 "Serbia", "Montenegro", "Bosnia and Herz.",
                                                 "Hungary", "Croatia", "Slovenia", 
                                                 "Italy", "Austria", "Kosovo"),]
poly_temp <- terra::buffer(poly_temp, 0.1)
poly_temp <- gUnion(poly_temp, terra::buffer(countries50[countries50$name == "Moldova",], 1.5))

range_temp_moreotica <- raster::intersect(range_temp, poly_temp)
range_temp_moreotica@data$subgroup <- "moreotica ecsedensis"

# aureus and syriacus
range_temp_aureus <- range_temp - clip_poly_temp - poly_temp
range_temp_aureus@data$subgroup <- "aureus syriacus"

# adding it all up
range_temp_try <- raster::bind(range_temp_indicus, range_temp_moreotica, range_temp_aureus)

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("aureus syriacus", "moreotica ecsedensis"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

# Following IUCN range description I'm not including the algirensis subspecies (sample in Tunisia)
species_temp <- "Canis aureus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

sp.gmpd.points <- sp.gmpd.points[range_temp,]
GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Canis latrans ####
# although subspecies ranges are well described in wiki, there's overlap between them 
# and I don't know exactly where to draw lines. 
# Also range is pretty continuous (apart from tiburon island, which isn't in polygons anyway). 
# The only subspecies I'm lacking are central american ones which don't appear geographically 
# isolated from descriptions or geographic boundaries

#map# species_temp <- "Canis latrans"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Canis lupus ####
# for Canis lupus rufus / Canis rufus:
"See Chambers et al. (2012) for a brief review of recent literature concerning the status of this species, 
which they considered a full species, as does this assessment."

# Researching individual subspecies from IUCN taxonomic notes 
## other general sources: 
# https://en.wikipedia.org/wiki/Subspecies_of_Canis_lupus and 
# https://en.wikipedia.org/wiki/List_of_gray_wolf_populations_by_country

species_temp <- "Canis lupus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

# North America
# mainland population is generally continuous apart from baileyi in Arizona/New Mexico
# arctos is geographically isolated and not represented
## remaining subspecies are lycaon, occidentalis, and nubilis
## uncertain about geographic boundaries between species so will leave as one
clip_poly_temp <- countries50[countries50$continent == "North America",]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.4) 

range_temp_america <- raster::intersect(range_temp, clip_poly_temp)
range_temp_america@data$subgroup <- "lycaon occidentalis nubilis"

range_temp_america@data <- range_temp_america@data %>%
  mutate(subgroup = case_when(str_detect(island, "Greenland|Ellesmere|Banks|Melville") ~ "arctos",
                              island == "Baffin" ~ "manningi",
                              island == "Vancouver" ~ "crassodon",
                              TRUE ~ subgroup))

# baileyi
range_temp_baileyi <- range_temp_america[states50[states50$name %in% c("Arizona", "New Mexico"),],]
range_temp_baileyi@data$subgroup <- "baileyi"  

# adding america together
range_temp_america <- range_temp_america - states50[states50$name %in% c("Arizona", "New Mexico"),]
range_temp_america <- raster::bind(range_temp_america, range_temp_baileyi)

#map# tm_shape(range_temp_america) + tm_polygons("subgroup")

# Europe:
# signatus is isolated (https://doi.org/10.1111/mec.14824) (Iberian population)
# italicus is isolated (https://doi.org/10.1016/j.mambio.2017.01.005) and does not hybridise much with domestic dogs (https://doi.org/10.1007/BF03194151)

range_temp_lupus <- range_temp - clip_poly_temp
range_temp_lupus@data$subgroup <- "lupus"

# signatus
clip_poly_temp <- countries50[countries50$name %in% c("Spain", "Portugal"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.2)
range_temp_signatus <- raster::intersect(range_temp_lupus, clip_poly_temp)

clip_poly_temp <- countries50[countries50$name == "Andorra",]
clip_poly_temp <- terra::buffer(clip_poly_temp, 1)
range_temp_signatus <- range_temp_signatus - clip_poly_temp
range_temp_signatus@data$subgroup <- "signatus"
range_temp_lupus <- range_temp_lupus - range_temp_signatus

# italicus
clip_poly_temp <- countries50[countries50$name %in% c("France", "Italy", "Andorra"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.08)
# clip_poly_temp <- clip_poly_temp - countries50[countries50$name %in% c("Belgium", "Germany",
#                                                              "Austria", "Switzerland"),]
poly_temp <- terra::buffer(countries50[countries50$name == "Andorra",], 1)
clip_poly_temp <- gUnion(clip_poly_temp, poly_temp) 

poly_temp <- matrix(c(8.919, 45.885,
                      8.919, 47.33,
                      14.64, 47.33,
                      14.64, 44.46,
                      10.56, 45.3,
                      8.919, 45.885), 
                    ncol = 2, byrow = TRUE)
poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - poly_temp

range_temp_italicus <- raster::intersect(range_temp_lupus, clip_poly_temp)
range_temp_italicus@data$subgroup <- "italicus"
range_temp_lupus <- range_temp_lupus - range_temp_italicus

# Asia:
## very few asian subspecies represented within gmpd, and those that are do not appear to be isolated
## Subpopulations are also less well documented generally 
# isolated arabs population in Arabian peninsula (Yemen, Oman, Southern Saudi Arabia)
# isolated arabs and pallipes population in Sinai peninsula, Israel, Jordan, Lebanon, Southern Syria
# isolated pallipes population in Northwest India

# pallipes and chanco
clip_poly_temp <- matrix(c(29, 40.8,
                      29, 41.03,
                      29.05, 41.06, 
                      29.065, 41.09,
                      29.08, 41.13,
                      29.07, 41.16,
                      29.16, 41.26,
                      35, 43,
                      40.01, 43.38,
                      40.01, 43.41, 
                      52, 39, 
                      54.02, 35.67, 
                      62.76, 34.43, 
                      70.45, 37.29, 
                      75.16, 42, 
                      80.21, 42.13, 
                      98.49, 48.37, 
                      
                      100.54, 51.594,
                      103.723, 51.706,
                      104.4, 51.6,
                      105.8, 51.98,
                      106.16, 52.36,
                      108.006, 53.159,
                      109.676, 55.689,
                      112.093, 56.23,
                      113.016, 56.205,
                      114.422, 56.254,
                      116.531, 56.352,
                      117.718, 56.57,
                      119.168, 57.004,
                      121.453, 56.98,
                      123.87, 56.401,
                      126.375, 55.515,
                      128.902, 54.692,
                      131.649, 54.039,
                      133.671, 54.3,
                      
                      135.7, 54.9,
                      147, 54.9,
                      147, 6.7,
                      27, 6.7, # south
                      27, 36.5,
                      25, 39.7,
                      26.3, 40.05,
                      26.4, 40.14,
                      26.4, 40.2,
                      26.8, 40.45,
                      29, 40.8),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

poly_temp <- terra::buffer(countries50[countries50$name %in% c("Georgia", "Azerbaijan", "Iran", 
                                                               "Afghanistan", "Tajikistan", "Kyrgyzstan", 
                                                               "China", "Mongolia"),], 0.09)
clip_poly_temp <- gUnion(clip_poly_temp, poly_temp)

poly_temp <- matrix(c(70.508, 40.961,
                     70.954, 41.255,
                     71.684, 41.683,
                     73.651, 40.928,
                     71.162, 39.582,
                     70.069, 40.445,
                     70.508, 40.961), 
                   ncol = 2, byrow = TRUE)

poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- gUnion(clip_poly_temp, poly_temp)

range_temp_pallipes <- raster::intersect(range_temp_lupus, clip_poly_temp)
range_temp_pallipes@data$subgroup <- "pallipes chanco"
range_temp_lupus <- range_temp_lupus - range_temp_pallipes

# arabs
clip_poly_temp <- countries50[countries50$name %in% c("Saudi Arabia", "Bahrain", "Qatar", "United Arab Emirates", "Oman", "Yemen"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.5)
clip_poly_temp <- clip_poly_temp - terra::buffer(countries50[countries50$name %in% c("Egypt", "Israel", "Jordan", 
                                                              "Iraq", "Kuwait"),], 0.05)
range_temp_arabs <- raster::intersect(range_temp_pallipes, clip_poly_temp)
range_temp_arabs@data$subgroup <- "arabs"
range_temp_pallipes <- range_temp_pallipes - range_temp_arabs

# adding all up 
range_temp_try <- raster::bind(range_temp_america, range_temp_lupus, range_temp_signatus, range_temp_italicus, range_temp_pallipes, range_temp_arabs)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
buff_temp <- data.frame(subgroup = c("signatus", "italicus", "lupus", "pallipes chanco", "lycaon occidentalis nubilis"),
                                 buff = c(1.5, 0, 0, 0, 1))
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = buff_temp)

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))

## Canis mesomelas ####
# two geographically distinct subspecies, already done :)
species_temp <- "Canis mesomelas"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("mesomelas", "schmidti"), buff = c(0.5, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Canis simensis ####
# Two geographically distinct subspecies and I only seem to have citernii
# https://en.wikipedia.org/wiki/Ethiopian_wolf#/media/File:Canis_simensis_subspecies_range.png
species_temp <- "Canis simensis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

clip_poly_temp = matrix(c(37, 9,
                  40, 9,
                  40, 14,
                  37, 14,
                  37, 9), 
                ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- range_temp - clip_poly_temp
range_temp_try@data$subgroup <- "citernii"

range_temp_otherssp <- raster::intersect(range_temp, clip_poly_temp)
range_temp_otherssp@data$subgroup <- "simensis"

range_temp_try <- raster::bind(range_temp_try, range_temp_otherssp)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "citernii"

rm(list = ls(pattern = "^range_temp"))

## Capra ibex ####
# no reported subspecies in iucn, wiki, or msotw
#map# species_temp <- "Capra ibex"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Capra pyrenaica ####
# four subspecies, two extinct. I probably only have hispanica which is geographically distinct from victoriae
# subdivided based on IUCN description of ranges
species_temp <- "Capra pyrenaica"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

range_temp_try <- raster::disaggregate(range_temp)
range_temp_try@data$temp <- letters[1:nrow(range_temp_try@data)]
#map# tm_shape(range_temp_try) + tm_polygons("temp")

split_hispanica <- c("a", "b", "c", "d","e", "f", "g", "i", "m", "n")
range_temp_try@data <- range_temp_try@data %>%
  mutate(subgroup = case_when(temp %in% split_hispanica ~ "hispanica",
                              TRUE ~ "victoriae")) %>% 
  dplyr::select(-temp)

range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "hispanica"

## Capreolus capreolus ####
# The IUCN notes the following five confirmed subspecies: italicus, garganta, capreolus, caucasicus, and coxi
# Definitely have italicus and garganta
# uncertain about capreolus but I expect so, plus it's not geographically isolated
# Don't appear to have coxi and caucasicus and they are geographically separated

species_temp <- "Capreolus capreolus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

range_temp_try <- raster::disaggregate(range_temp)

# coxi
clip_points_temp <- matrix(c(34.7, 41.4,
                        27.6, 38.3,
                        36.4, 37.1,
                        36.4, 34.5),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_coxi <- range_temp_try[clip_points_temp,]
range_temp_coxi@data$subgroup <- "coxi"

# caucasicus
clip_points_temp <- matrix(c(44.2, 41.4),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_caucasicus <- range_temp_try[clip_points_temp,]
range_temp_caucasicus@data$subgroup <- "caucasicus"

# italicus garganta capreolus
clip_points_temp <- matrix(c(54.045, 36.89),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_remove <- range_temp_try[clip_points_temp,]

range_temp_capreolus <- range_temp_try - range_temp_caucasicus - range_temp_coxi - range_temp_remove
range_temp_capreolus@data$subgroup <- "italicus garganta capreolus"

range_temp_try <- raster::bind(range_temp_coxi, range_temp_caucasicus, range_temp_capreolus)
range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "italicus garganta capreolus"

rm(list = ls(pattern = "^range_temp"))

## Capricornis crispus ####
# No subspecies noted by IUCN, wiki, msotw
#map# species_temp <- "Capricornis crispus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Cephalophus natalensis ####
# Two subspecies have been named: C. n. natalensis and C. n. robertsi (north of the Limpopo river)
# I only seem to have natalensis

species_temp <- "Cephalophus natalensis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

L_temp <- River_Data50[(River_Data50$name == "Limpopo"),]
L_coords <- unlist(coordinates(L_temp), recursive = FALSE)[[1]]
L_coords <- rbind(L_coords,
                matrix(c(33.53, -25.2,
                         33.53, (range_temp@bbox[2,1] - 1),
                         (range_temp@bbox[1,1] - 1), (range_temp@bbox[2,1] - 1),
                         (range_temp@bbox[1,1] - 1), L_coords[1,2],
                         L_coords[1,]), 
                       ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(L_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_natalensis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_natalensis@data$subgroup <- "natalensis"

range_temp_robertsi <- range_temp - clip_poly_temp
range_temp_robertsi@data$subgroup <- "robertsi"

range_temp_try <- raster::bind(range_temp_natalensis, range_temp_robertsi)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "natalensis"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "^L_"))
rm(list = ls(pattern = "_temp$"))

## Cerdocyon thous ####
# wiki lists 5 subspecies: thous, azarae, entrerianus, aquilus, germanus
# I seem to have azarae and entrerianus and the others are geographically separate

species_temp <- "Cerdocyon thous"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# azarae entrerianus
range_temp_try <- raster::disaggregate(range_temp)
clip_points_temp <- matrix(c(-49.5, -18.6),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_azarae <- range_temp_try[clip_points_temp,]
range_temp_azarae@data$subgroup <- "azarae entrerianus"

# thous aquilus germanus
range_temp_try <- range_temp - range_temp_azarae
range_temp_try@data$subgroup <- "thous aquilus germanus"

range_temp_try <- raster::bind(range_temp_try, range_temp_azarae)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "azarae entrerianus"

rm(list = ls(pattern = "^range_temp"))

## Cervus canadensis ####
# I only have north american subspecies: canadensis, nannodes, roosevelti
# Not clear how to separate subspecies even though there are subpopulations so just going by continent
"Here we recognise the subspecies and their distributions as follows:

C. c. canadensis – N America
C. c. alashanicus – N China
C. c. nannodes - California
C. c. roosevelti - Vancouver Island, Washington state and Oregon
C. c. sibiricus – NE Kazakhstan and N Xinjiang to S Siberia and N Mongolia
C. c. xanthopygus – SE Siberia, Russian Far East, Ussuriland, Manchuria
C. c. macneilli – Lydekker 1909 (Central and SW China (N Qinghai, Gansu, Shaanxi, W Sichuan and E Xizang)
C. c. wallichii - G. Cuvier 1823 (SW China (SE Xizang), Bhutan)"

species_temp <- "Cervus canadensis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

clip_poly_temp <- matrix(c(-125, 61,
                      -77,61,
                      -77,25,
                      -125,25,
                      -125,61),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_americas <- raster::intersect(range_temp, clip_poly_temp)
range_temp_americas@data$subgroup <- "canadensis nannodes roosevelti"

range_temp_otherssp <- range_temp - clip_poly_temp
range_temp_otherssp@data$subgroup <- "other"
range_temp_try <- raster::bind(range_temp_americas, range_temp_otherssp)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial$HostCorrectedName == species_temp, "subgroup"] <- "canadensis nannodes roosevelti"

# Extra GMPD bit
# Some samples recorded as Cervus elaphus 
"Until recently, biologists considered the red deer and elk or wapiti (C. canadensis) the same species, 
forming a continuous distribution throughout temperate Eurasia and North America. This belief was based 
largely on the fully fertile hybrids that can be produced under captive conditions.[27][28][29]"
elaphus_temp <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == "Cervus elaphus", ]
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(elaphus_temp) + tm_dots("subgroup")

range_temp_temp <- gBuffer(range_temp_try[range_temp_try$subgroup != "other",], byid = TRUE, width = 4)
clip_points_temp <- matrix(c(-148.5, 64),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- SpatialPolygonsDataFrame(terra::buffer(clip_points_temp, 100000), 
                                        data = range_temp_temp@data,
                                        match.ID = FALSE)

range_temp_temp <- raster::bind(range_temp_temp, clip_points_temp)

elaphus_temp <- elaphus_temp[range_temp_temp,]
elaphus_temp@data$HostCorrectedName <- species_temp
elaphus_temp@data$subgroup <- "canadensis nannodes roosevelti"
GMPD_Spatial <- raster::bind(GMPD_Spatial, elaphus_temp)

range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Cervus elaphus ####
# I seem to have: elaphus, italicus, and montanus
# but not: barbarus, corsicanus, maral, brauneri
# Other sources (msotw and wiki) also record scoticus in British Isles, and atlanticus in Norway/South Sweden
"Several subspecies of Western Red Deer have been recognized with their ranges as follows:
C. e. elaphus: Ireland, Great Britain, continental Europe
C. e. barbarus: Atlas Mountains (Algeria, Tunisia)
C. e. corsicanus: Corsica (extinct, reintroduced in 1985), Sardinia
C. e. maral: Anatolia,
C. e. italicus: Italy (Ferrara)
C. e. brauneri Crimea (Russia)
C. e. montanus (syn. Carpathicus) Carpathian mountains" #IUCN

# Not sure which subspecies are around Georgia/East Turkey/Iran, or south of Finland (in Russia) from IUCN description
# But based on wiki they are maral and elaphus, respectively

"However, recent analyses call into question the veracity of the morphology-based subspecific taxonomy of C. elaphus 
more generally. Rather, three genetic lineages were clearly differentiated corresponding approximately with geographical 
factors. One lineage was distributed in central-western Europe (and also Ukraine), one from eastern Europe to the 
Middle East, and a third which corresponded to C. e. barbarus and C. e. corsicanus from North Africa and Sardinia, 
the latter two nomials are considered synonyms and were recommended considered as a subspecies"

# will group montanus in with elaphus as it is contiguous, has no clear boundaries, and is not always treated as separate subspecies from elaphus 
# but keep scoticus and atlanticus separate as they're isolated.
# Though there have been European introductions to England

species_temp <- "Cervus elaphus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# British Isles (scoticus)
buffer_temp <- countries50[countries50$name %in% c("United Kingdom", "Ireland"),]
buffer_temp <- terra::buffer(buffer_temp, 0.3)
range_temp_scoticus <- raster::intersect(range_temp, buffer_temp)
range_temp_scoticus@data$subgroup <- "scoticus"

# North Africa (barbarus)
buffer_temp <- countries50[countries50$name == "Tunisia",]
buffer_temp <- terra::buffer(buffer_temp, 0.3)
buffer_temp <- buffer_temp + countries50[countries50$name == "Morocco",]
range_temp_barbarus <- raster::intersect(range_temp, buffer_temp)
range_temp_barbarus@data$subgroup <- "barbarus"

# Iberian peninsula (hispanicus)
buffer_temp <- countries50[countries50$name %in% c("Spain", "Portugal", "Andorra"),]
buffer_temp <- terra::buffer(buffer_temp, 0.06)
range_temp_hispanicus <- raster::intersect(range_temp, buffer_temp)
range_temp_hispanicus@data$subgroup <- "hispanicus"

# Corsica and Sardinia (corsicanus)
range_temp_corsicanus <- range_temp[range_temp$island %in% c("Corsica", "Sardinia"),]
range_temp_corsicanus@data$subgroup <- "corsicanus"

# Norway and Sweden (atlanticus)
buffer_temp <- countries50[countries50$name %in% c("Norway", "Sweden"),]
buffer_temp <- terra::buffer(buffer_temp, 0.3)
range_temp_atlanticus <- raster::intersect(range_temp, buffer_temp)
range_temp_atlanticus@data$subgroup <- "atlanticus"

# Italy (italicus)
clip_poly_temp <- matrix(c(18, 40,
                      16, 39,
                      7, 45,
                      13, 45,
                      18, 40),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_italicus <- raster::intersect(range_temp, clip_poly_temp)
range_temp_italicus@data$subgroup <- "italicus"

# Crimea (brauneri)
clip_poly_temp <- matrix(c(32, 46,
                      36, 46,
                      36, 44,
                      32, 44,
                      32, 46),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_brauneri <- raster::intersect(range_temp, clip_poly_temp)
range_temp_brauneri@data$subgroup <- "brauneri"

# Anatolia, Caucasus, Iran (maral)
clip_poly_temp <- matrix(c(28.6, 40.7,
                      37, 46,
                      57, 46,
                      57, 35,
                      30, 35,
                      25.5, 39.5,
                      27.5, 40.7,
                      28.6, 40.7),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_maral <- raster::intersect(range_temp, clip_poly_temp)
range_temp_maral@data$subgroup <- "maral"

# Most of Europe (elaphus), including Carpathian mountains (montanus)
range_temp_try <- raster::bind(range_temp_scoticus, range_temp_barbarus, range_temp_hispanicus, 
                         range_temp_corsicanus, range_temp_atlanticus, range_temp_italicus, 
                         range_temp_brauneri, range_temp_maral)
range_temp_try <- range_temp_try[, 1:29]
range_temp_elaphus <- range_temp - range_temp_try
range_temp_elaphus@data$subgroup <- "elaphus"

range_temp_try <- raster::bind(range_temp_try, range_temp_elaphus)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# original data has some hispanicus samples for which info is retained, and hippelaphus which is replaced with elaphus
clip_points_temp <- matrix(c(7, 44.5),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- SpatialPolygonsDataFrame(terra::buffer(clip_points_temp, 100000), 
                                        data = range_temp_try@data[range_temp_try$subgroup == "elaphus" & range_temp_try$origin == 1,],
                                        match.ID = FALSE)

range_temp_try <- raster::bind(range_temp_try, clip_points_temp)

buff_temp <- data.frame(subgroup = c("scoticus", "hispanicus", "elaphus", "italicus", "atlanticus"), 
                                 buff = c(0.1, 0, 0, 0, 0.5))
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = buff_temp)

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Cervus nippon ####
# Sources:
# https://doi.org/10.1007/s10344-005-0011-5            Groves 2005
## Groves describes the North/South differences between japanese subspecies as sufficient for treatment as separate species
## So although there is no apparent barrier between subspecies, I'm splitting them according to descriptions in Groves and Goodman
# https://doi.org/10.1046/j.1365-294X.2001.01277.x     Goodman et al 2001
# IUCN describes species status of yesoensis and hortulorum as uncertain so maintain as one species under nippon

# splitting at the Hyogo prefecture Eastern boundary based on descriptions in Groves and Goodman 
# only have samples on main Japanese islands so excluding small island subspecies as well as mainland asia subspecies

species_temp <- "Cervus nippon"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

clip_poly_temp <- matrix(c(134.87, 35.8,
                      134.87, 35.66,
                      134.85, 35.59,
                      134.92, 35.54,
                      134.93, 35.51,
                      135.04, 35.54,
                      135.05, 35.41,
                      134.99, 35.39,
                      134.94, 35.41,
                      134.92, 35.31,
                      135.39, 35.1,
                      135.35, 34.96,
                      135.47, 34.92,
                      135.42, 34.82,
                      135.46, 34.74,
                      135.41, 34.69,
                      134.8, 33.8,
                      134.8, 30.9,
                      129.5, 30.9,
                      129.5, 33.3,
                      130.8, 34.5,
                      132.8, 35.7,
                      134.87, 35.8),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
#map# tm_shape(clip_poly_temp) + tm_polygons(alpha = 0) + tm_shape(range_temp) + tm_polygons(alpha = 0)

range_temp_nippon <- raster::intersect(range_temp, clip_poly_temp)
range_temp_nippon@data$subgroup <- "nippon"

range_temp_aplodontus <- range_temp[range_temp$island == "Honshu" & !is.na(range_temp$island),] - clip_poly_temp - countries50[countries50$name %in% c("Russia", "China"),]
range_temp_aplodontus@data$subgroup <- "aplodontus"

range_temp_yesoensis <- range_temp[range_temp$island == "Hokkaido" & !is.na(range_temp$island),]
range_temp_yesoensis@data$subgroup <- "yesoensis"

range_temp_try <- raster::bind(range_temp_nippon, range_temp_aplodontus, range_temp_yesoensis)
range_temp_otherssp <- range_temp - range_temp_try
range_temp_try <- raster::bind(range_temp_otherssp, range_temp_try)
range_temp_try@data[is.na(range_temp_try@data$subgroup), "subgroup"] <- "other"
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# original subgroup designations correspond to polygons :)
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("nippon", "aplodontus", "yesoensis"), buff = c(0.5, 0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))

## Chrysocyon brachyurus ####
# no subspecies described by IUCN, wiki, or msotw
#map# species_temp <- "Chrysocyon brachyurus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Civettictis civetta #### 
# IUCN doesn't describe any subspecies but they are recognised on wiki and msotw
# Subspecific range descriptions lack detail, are merely type specimen locations
# Range appears continuous and hard to subdivide based on biogeographic barriers

species_temp <- "Civettictis civetta"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Conepatus chinga ####
# No subspecies described by IUCN or wiki, and range descriptions not available on msotw
# Also range has no clearly isolated populations
#map# species_temp <- "Conepatus chinga"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Connochaetes gnou ####
# No subspecies described by IUCN, wiki, or msotw
#map# species_temp <- "Connochaetes gnou"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Connochaetes taurinus ####
# already included in species range polygons <3
species_temp <- "Connochaetes taurinus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("mearnsi", "albojubatus", "taurinus"), buff = c(0, 0, 1)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Crocuta crocuta ####
# neither IUCN nor msotw list any subspecies and wiki states:
# "all the variation seen in the then recognised subspecies could also be found in a single population"
#map# species_temp <- "Crocuta crocuta"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Cynictis penicillata ####
# no subspecies listed by IUCN, no range descriptions on wiki, and the 12 listed on msotw have no range descriptions
# Also range appears continuous
#map# species_temp <- "Cynictis penicillata"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Damaliscus lunatus ####
# already included in species range polygons <3
species_temp <- "Damaliscus lunatus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# GMPD subgroup assignment
# Zambian sample already assigned lunatus which also corresponds with GBIF records
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("lunatus", "jimela"), buff = c(4, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Damaliscus pygargus ####
# One subspecies already included, pygargus is the other described by IUCN
species_temp <- "Damaliscus pygargus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == species_temp & is.na(subspecies) ~ "pygargus",
                              TRUE ~ subgroup))
#map# tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]) + tm_polygons("subgroup") 

# GMPD subgroup assignment
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("phillipsi", "pygargus"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Equus grevyi ####
# No subspecies on msotw
# "However, Groves and Bell (2004) concluded that the species is indeed monotypic."
#map# species_temp <- "Equus grevyi"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Equus quagga ####
# No subspecies on msotw
# "The molecular data represented a genetic cline" so it's monotypic
species_temp <- "Equus quagga"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Equus zebra ####
# "We continue to recognize Mountain Zebra as a single species comprising two subspecies."
# polygons already labelled but with rogue hartmannae polygons within zebra range
# this is despite IUCN describing them as allopatric subspecies
species_temp <- "Equus zebra"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# correcting polygons
clip_poly_temp <- Polygon(matrix(c(16.9, -29.5,  
                              28, -29.5,  
                              28, -35, 
                              16.9, -35, 
                              16.9, -29.5), 
                            ncol = 2, byrow = TRUE))
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_try <- raster::disaggregate(range_temp)
range_temp_try <- range_temp_try[clip_poly_temp,]
range_temp_try@data$subgroup <- "zebra"

range_temp_try <- raster::bind(range_temp - range_temp_try, range_temp_try)

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("hartmannae", "zebra"), buff = c(0.5, 0.5)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Felis silvestris and lybica ####
# "A revised taxonomy of the Felidae. The final report of the Cat Classification Task Force of the IUCN/SSC Cat Specialist Group. Cat News Special Issue 11, 80 pp."
# two subspecies (silvestris and caucasica) 
# and possibly a third which I will include because it's an island population (grampia)
# I seem to also have Felis libyca
species_temp <- "Felis silvestris"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Felis silvestris
# grampia
clip_poly_temp <- countries50[countries50$name == "United Kingdom",]
clip_poly_temp <- raster::buffer(clip_poly_temp, 0.2)
range_temp_grampia <- raster::intersect(range_temp, clip_poly_temp)
range_temp_grampia@data$subgroup <- "grampia"

# silvestris
turkey_rangepol <- countries50[countries50$name == "Turkey",]
turkey_rangepol <- raster::disaggregate(turkey_rangepol)
clip_poly_temp <- countries50[countries50$continent == "Europe" & !countries50$name %in% c("United Kingdom", "Ireland", "Russia"),]
clip_poly_temp <- raster::bind(clip_poly_temp, turkey_rangepol[3,])
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.1)
clip_poly_temp <- gUnion(clip_poly_temp, terra::buffer(countries50[countries50$name %in% c("Ukraine", "Croatia"),], 0.3))

range_temp_silvestris <- raster::intersect(range_temp, clip_poly_temp)
range_temp_silvestris@data$subgroup <- "silvestris"

# caucasica
clip_poly_temp <- matrix(c(36.7, 44.6,
                      37.4, 46.6,
                      49.7, 43.6,
                      49, 40.5,
                      46.5, 38.9,
                      46.2, 38.9,
                      44.7, 39.8,
                      42.7, 37.3,
                      36, 35.7,
                      26.2, 36.3,
                      26.3, 40.05,
                      36.7, 44.6),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_caucasica <- raster::intersect(range_temp, clip_poly_temp)
range_temp_caucasica@data$subgroup <- "caucasica"

range_temp_try <- raster::bind(range_temp_grampia, range_temp_silvestris, range_temp_caucasica)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# Sorting out GMPD 
silvestris_temp <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("grampia", "silvestris"), buff = c(0, 0.7)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], silvestris_temp)

## Felis lybica
range_temp <- range_temp - range_temp_try

clip_poly_temp <- matrix(c(7.9, 43.2,
                      10, 43.2,
                      10, 38.5,
                      7.9, 38.5,
                      7.9, 43.2),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_silvestris <- raster::disaggregate(range_temp_silvestris)
range_temp <- raster::bind(range_temp, range_temp_silvestris[clip_poly_temp,])
range_temp@data$binomial <- "Felis lybica"

# ornata
clip_poly_temp <- matrix(c(42.6, 48.6,
                      107.6, 48.6,
                      107.6, 17.5,
                      62.5, 14.8, 
                      56.65, 26.6,
                      50.2, 24.4,
                      47.2, 32.1,
                      42.6, 35.5,
                      42.6, 48.6),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_ornata <- raster::intersect(range_temp, clip_poly_temp)
range_temp_ornata@data$subgroup <- "ornata"

# cafra
clip_poly_temp <- matrix(c(10,2,
                      29.5, -6,
                      40.4, -10.5,
                      40.8, -10.2,
                      42, -35,
                      10, -35,
                      10,2),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - countries50[countries50$name == "Tanzania",]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.01)

range_temp_cafra <- raster::intersect(range_temp, clip_poly_temp)
range_temp_cafra@data$subgroup <- "cafra"

# lybica
range_temp_lybica <- range_temp - range_temp_ornata - range_temp_cafra
range_temp_lybica@data$subgroup <- "lybica"

range_temp_try <- raster::bind(range_temp_ornata, range_temp_cafra, range_temp_lybica)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data, range_temp_try)

# Sorting out GMPD
lybica_temp <- sp.gmpd.points[!is.na(sp.gmpd.points@data$subgroup),]
lybica_temp@data$HostCorrectedName <- "Felis lybica"
GMPD_Spatial <- raster::bind(GMPD_Spatial, lybica_temp)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Galictis cuja ####
# from wiki:
"Four subspecies are recognised:

Galictis cuja cuja – southwestern Bolivia, western Argentina, central Chile
Galictis cuja furax – southern Brazil, northeastern Argentina, Uruguay, and Paraguay
Galictis cuja huronax – south-central Bolivia, eastern Argentina
Galictis cuja luteola – extreme southern Peru, western Bolivia and northern Chile"

# However, range appears continuous and I have no knowledge of barriers between subpopulations
#map# species_temp <- "Galictis cuja"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Genetta genetta ####
# wiki:
"G. g. genetta (Linnaeus), 1758 — Spain, Portugal and France
G. g. afra (Cuvier), 1825 — North Africa[32]
G. g. senegalensis (Fischer), 1829 — sub-Saharan Africa[33]
G. g. dongolana (Hemprich and Ehrenberg), 1832 — Arabia[34]"

# IUCN: 
"This assessment includes the South African Small-spotted Genet (Genetta felina (Thunberg, 1811))"

# I have genetta in Spain, but excluding because introduced
# I don't have afra or dongolana
# Southern African population is G. felina

species_temp <- "Genetta genetta"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

clip_poly_temp <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia", "Libya"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.1)
range_temp_afra <- raster::intersect(range_temp, clip_poly_temp)
range_temp_afra@data$subgroup <- "afra"

clip_poly_temp <- countries50[countries50$name %in% c("South Africa", "Botswana", "Namibia", "Angola", "Zimbabwe"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 1.7)
range_temp_felina <- raster::intersect(range_temp, clip_poly_temp)
range_temp_felina@data$subgroup <- "felina"

clip_poly_temp <- countries50[countries50$name %in% c("Saudi Arabia", "Yemen", "Oman"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.2)
range_temp_dongolana <- raster::intersect(range_temp, clip_poly_temp)
range_temp_dongolana@data$subgroup <- "dongolana"

range_temp_try <- raster::bind(range_temp_afra, range_temp_felina, range_temp_dongolana)
range_temp_senegalensis <- range_temp - range_temp_try
range_temp_senegalensis@data$subgroup <- "senegalensis"

range_temp_try <- raster::bind(range_temp_try, range_temp_senegalensis)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

sp.gmpd.points <- sp.gmpd.points[!is.na(sp.gmpd.points@data$subgroup),]
GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Genetta thierryi ####
# No subspecies listed by IUCN, wiki, or msotw
# plus it has a small continuous range
#map# species_temp <- "Genetta thierryi"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Giraffa camelopardalis ####
# already recorded <3
species_temp <- "Giraffa camelopardalis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# GMPD subgroup assignment
# samples which appear to be within range of rothschildi are already recorded as reticulata
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("reticulata", "angolensis", "giraffa"), buff = c(1, 0, 1)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Herpestes ichneumon ####
# none listed on IUCN, but 11 listed by msotw (location notes from wiki)

# SUBSPECIES ichneumon      from the area of the Nile River in Egypt
# SUBSPECIES angolensis     a male specimen from Quissange in Angola
# SUBSPECIES cafra          based on a description of a specimen from the Cape of Good Hope
# SUBSPECIES centralis      two specimens from Beni, Democratic Republic of the Congo
# SUBSPECIES funestus       a specimen from Naivasha in British East Africa (Kenya)
# SUBSPECIES mababiensis    a specimen from Mababe in northern Bechuanaland (Botswana)
## SUBSPECIES numidicus      two individuals from Algiers in Algeria
# SUBSPECIES parvidens      three specimens collected near the lower Congo River in Congo Free State (roughly DRC)
# SUBSPECIES sabiensis      specimen from Sabi Sand Game Reserve in Southern Africa
## SUBSPECIES sangronizi     specimen from Mogador in Morocco
## SUBSPECIES widdringtonii  specimen from Sierra Morena in Spain

# I might have nominate subspecies, but hard to tell
# Much of range is continuous, apart from Iberian Peninsula and Northwest Africa which I will exclude
# Not possible to justify any further subdivision of range based on available information

species_temp <- "Herpestes ichneumon"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

# Northwest Africa (sangronizi and numidicus) and Iberian peninsula (widdringtonii)
clip_poly_temp <- matrix(c(-16, 44,
                      18, 44,
                      18, 19,
                      -16, 19,
                      -16, 44),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_otherssp <- raster::intersect(range_temp, clip_poly_temp)
range_temp_otherssp@data$subgroup <- "sangronizi numidicus widdringtonii"
range_temp_ichneumon <- range_temp - clip_poly_temp
range_temp_ichneumon@data$subgroup <- "ichneumon"

range_temp_try <- raster::bind(range_temp_otherssp, range_temp_ichneumon)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "ichneumon"

rm(list = ls(pattern = "^range_temp"))

## Hippopotamus amphibius ####
# three subspecies in msotw, 
# though it's possible tschadensis and constrictus are missing because of lack of info rather than anything else
"Genetic analyses have tested the existence of three of these putative subspecies. A study examining 
mitochondrial DNA from skin biopsies taken from 13 sampling locations, considered genetic diversity and 
structure among hippo populations across the continent. The authors found low, but significant, genetic 
differentiation among H. a. amphibius, H. a. capensis, and H. a. kiboko. Neither H. a. tschadensis nor 
H. a. constrictus has been tested.[9][10]" #(wiki)
# I likely only have capensis, described in wiki as from Zambia to South Africa.
# However, it's really hard to tell if there's admixture because their population is so subdivided
# Also subspecies ranges described in wiki have a lot of overlap

#map# species_temp <- "Hippopotamus amphibius"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Hippotragus niger ####
# Only variani is geographically separate, and descriptions of ranges do not specify biogeographical boundaries
"Four subspecies are usually recognized: H. n. niger, H. n. kirkii, H. n. roosevelti and the isolated Giant Sable 
(H. n. variani) from Angola. As for many other antelope species, the validity and precise distribution of most of the 
described subspecies are uncertain. An extensive study of the geographical genetic structure of Hippotragus niger 
identified three genetic subdivisions representing a Kenya and east Tanzania clade (H. n. roosevelti), a west Tanzania 
clade (H. n. kirkii), and a southern African clade (H. n. niger) (Pitra et al. 2002)." #(IUCN)

species_temp <- "Hippotragus niger"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data[IUCN_Native_Data$binomial == species_temp & is.na(IUCN_Native_Data$subgroup), "subgroup"] <- "niger kirkii roosevelti"
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "niger kirkii roosevelti"

## Hyaena hyaena ####
# As of 2005,[3] no subspecies are recognised. (wiki, from msotw)
#map# species_temp <- "Hyaena hyaena"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Kobus ellipsiprymnus ####
# already labelled <3
species_temp <- "Kobus ellipsiprymnus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("defassa", "ellipsiprymnus"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Kobus kob ####
# already done <3
species_temp <- "Kobus kob"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "kob"

## Kobus leche ####
# dooone <£
#map# species_temp <- "Kobus leche"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Leopardus geoffroyi ####
# according to cat specialist group there's only one ssp
#map# species_temp <- "Leopardus geoffroyi"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Leopardus pardalis ####
# cat specialist group say there's two: pardalis (north) and mitis (south)
# "borders between subspecies are speculative" and I have one right on the border
"The morphological differentiation between Central and South American forms is clear and supported partly 
by molecular data as well as a clear biogeographical barrier, the Andes."
# I'll follow their apparent division at the border between Panama and Costa Rica

species_temp <- "Leopardus pardalis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

buffer_temp <- countries50[countries50$continent == "North America" & countries50$name != "Panama",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with

clip_points_temp <- matrix(c(-87.6, 13.3),
                      ncol = 2, byrow = TRUE) # patching up gaps
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
buffer_temp <- gUnion(buffer_temp, terra::buffer(clip_points_temp, 10000))

range_temp_pardalis <- raster::intersect(range_temp, buffer_temp)
range_temp_pardalis@data$subgroup <- "pardalis"
range_temp_mitis <- range_temp - buffer_temp
range_temp_mitis@data$subgroup <- "mitis"

range_temp_try <- raster::bind(range_temp_pardalis, range_temp_mitis)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# extending range_temp_try because of overlapping buffers
clip_points_temp <- matrix(c(-99.5, 32),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- SpatialPolygonsDataFrame(terra::buffer(clip_points_temp, 100000), 
                                        data = range_temp_pardalis@data)

range_temp_try <- raster::bind(range_temp_try, clip_points_temp)

sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try,
                                buff = data.frame(subgroup = c("mitis", "pardalis"), buff = c(0, 0.5)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Leopardus tigrinus ####
# "Until then L. tigrinus is recognised as having two subspecies: tigrinis (south) and oncilla" (IUCN)
# They are geographically distinct and I only have one
species_temp <- "Leopardus tigrinus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

range_temp_try <- raster::disaggregate(range_temp)
range_temp_try@data$subgroup <- c("tigrinus", "oncilla")
range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "tigrinus"

## Lontra canadensis ####
# (wiki)
"L. c. canadensis (Schreber, 1777) – (eastern Canada, U.S., Newfoundland)
L. c. kodiacensis (Goldman, 1935) – (Kodiak Island, Alaska)
L. c. lataxina (Cuvier, 1823) – (U.S.)
L. c. mira (Goldman, 1935) – (Alaska, British Columbia)
L. c. pacifica (J. A. Allen, 1898) – (Alaska, Canada, northern U.S., south to central California, northern Nevada, and northeastern Utah)
L. c. periclyzomae (Elliot, 1905) – (Queen Charlotte Islands, British Columbia)
L. c. sonora (Rhoads, 1898) – (U.S., Mexico)"

# I seem to have canadensis, pacifica, sonora, and lataxina
# Not kodiacensis, mira, periclyzomae
# Only kodiacensis and periclyzomae appear to be isolated
# Also, GBIF distribution of points and wiki distribution map show a more limited range, excluding Hawaiian islands and most of Nunavut
# I'm going to follow suit by removing Hawaii and islands North of Hudson Bay

species_temp <- "Lontra canadensis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

clip_poly_temp <- matrix(c(-78.5, 62.7,
                      -72.7, 62.9,
                      -63.4, 60.8,
                      -48.5, 48,
                      -79.4, 23.9,
                      -169.3, 23.9,
                      -169.3, 75.3,
                      -114.35, 68.97,
                      -113.2, 68.34,
                      -109.37, 68.37,
                      -106.75, 69.16,
                      -101.9, 68.24,
                      -95.79, 71.94,
                      -96.32, 74.27,
                      -88.5, 74.22,
                      -91.93, 71.83,
                      -87.09, 69.35,
                      -85.4, 69.95,
                      -83.1, 69.81,
                      -82.56, 69.71,
                      -81.34, 69.75, 
                      -78.5, 62.7),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- raster::intersect(range_temp, clip_poly_temp)

range_temp_kodiacensis <- range_temp_try[range_temp_try$island %in% c("Kodiak", "Afognak", "Shuyak"),]
range_temp_kodiacensis@data$subgroup <- "kodiacensis"

range_temp_periclyzomae <- range_temp_try[range_temp_try$island %in% c("Langara", "Graham", "Moresby", "Louise"),]
range_temp_periclyzomae@data$subgroup <- "periclyzomae"

range_temp_canadensis <- range_temp_try[!(range_temp_try$island %in% c("Kodiak", "Afognak", "Shuyak", "Langara", "Graham", "Moresby", "Louise")),]
range_temp_canadensis@data$subgroup <- "canadensis pacifica sonora lataxina mira"

range_temp_try <- raster::bind(range_temp_kodiacensis, range_temp_periclyzomae, range_temp_canadensis)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "canadensis pacifica sonora lataxina mira"

rm(list = ls(pattern = "^range_temp"))

## Lutra lutra ####
# following:
# https://doi.org/10.1093/mspecies/sew011
# I appear to have lutra only which appears isolated from other mainland asian subspecies by mountains at 
# Tajikistan/Kyrgyzstan border, and similarly around Turkey's Eastern border (though less distinct here and 
# use of border is approximate)
# Too hard to split aurobrunnea, kutab, and monticolus from each other (doesn't matter anyway though)

species_temp <- "Lutra lutra"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# angustifrons
clip_poly_temp <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.1)
range_temp_angustifrons <- raster::intersect(range_temp, clip_poly_temp)
range_temp_angustifrons@data$subgroup <- "angustifrons"

# nair
clip_poly_temp <- countries50[countries50$name == "Sri Lanka",]
clip_poly_temp <- terra::buffer(clip_poly_temp, 5)
range_temp_nair <- raster::intersect(range_temp, clip_poly_temp)
range_temp_nair@data$subgroup <- "nair"
##map# tm_shape(clip_poly_temp) + tm_polygons(alpha = 0) + tm_shape(range_temp) + tm_polygons(alpha = 0)

# barang
clip_poly_temp <- matrix(c(97.4, 23.8,
                      98.2, 24.3, 
                      108, 24.6,
                      108.05, 21.55,
                      108.05, 18.8,
                      110.2, 13.2, 
                      110.2, -7.9,
                      94.3, -7.9,
                      94.8, 23.1,
                      97.4, 23.8),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - countries50[countries50$name == "China",]

range_temp_barang <- raster::intersect(range_temp, clip_poly_temp)
range_temp_barang@data$subgroup <- "barang"

# hainana
range_temp_hainana <- range_temp[range_temp$island == "Hainan" & !is.na(range_temp$island),]
range_temp_hainana@data$subgroup <- "hainana"

# chinensis
clip_poly_temp <- matrix(c(97.4, 23.8,
                      97.55, 28.5, 
                      97.55, 37.8,
                      123.1, 37.8,
                      123.1, 20.2,
                      97.4, 20.2, 
                      97.4, 23.8),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - range_temp_barang - countries50[countries50$name == "Myanmar",]

range_temp_chinensis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_chinensis@data$subgroup <- "chinensis"

# aurobrunnea x kutab x monticolus
clip_poly_temp <- matrix(c(101, 38.2,
                      101, 21.2,
                      71.8, 21.2,
                      71.8, 38.2,
                      101, 38.2),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- gUnion(clip_poly_temp, countries50[countries50$name == "Tajikistan",]) - range_temp_chinensis - range_temp_barang

poly_temp <- matrix(c(69.45, 39.7,
                      68.3, 38.1,
                      66.1, 40.1,
                      70.4, 41.5, 
                      71.6, 39.7,
                      69.45, 39.7),
                    ncol = 2, byrow = TRUE)
poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - poly_temp

range_temp_akm <- raster::intersect(range_temp, clip_poly_temp)
range_temp_akm@data$subgroup <- "aurobrunnea kutab monticolus"

# seistanica
clip_poly_temp <- matrix(c(59.6, 28.7,
                      56.6, 38.5,
                      60.7, 44.6, 
                      69.2, 44.6,
                      69.3, 40,
                      70, 38.9,
                      72.3, 38.5, 
                      76.5, 28.7,
                      59.6, 28.7),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

clip_points_temp <- matrix(c(73.85, 34.88),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- terra::buffer(clip_points_temp, 100000)

clip_poly_temp <- clip_poly_temp - range_temp_akm - countries50[countries50$name %in% c("Kyrgyzstan", "Iran"),] - clip_points_temp

range_temp_seistanica <- raster::intersect(range_temp, clip_poly_temp)
range_temp_seistanica@data$subgroup <- "seistanica"

# meridionalis
clip_poly_temp <- matrix(c(44.1, 39.31,
                      44.05, 39.34, 
                      41.2, 39.3,
                      41.2, 43.5,
                      49.8, 43.5, 
                      51.7, 38.1,
                      60.2, 38.6, 
                      66.9, 32.4,
                      56.1, 26.3, 
                      40.5, 27.7,
                      46.3, 36.6,
                      44.1, 39.31),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

clip_points_temp <- matrix(c(62.38, 35.2),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- terra::buffer(clip_points_temp, 50000)

clip_poly_temp <- clip_poly_temp - countries50[countries50$name %in% c("Turkey", "Turkmenistan", "Afghanistan"),] - clip_points_temp

clip_points_temp <- matrix(c(53.92, 37.35, 
                        58.27, 37.64),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- terra::buffer(clip_points_temp, 50000)

clip_poly_temp <- gUnion(clip_poly_temp, clip_points_temp)

range_temp_meridionalis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_meridionalis@data$subgroup <- "meridionalis"

# lutra
range_temp_try <- raster::bind(range_temp_angustifrons, range_temp_meridionalis, range_temp_seistanica, range_temp_nair, range_temp_akm, range_temp_barang, range_temp_chinensis, range_temp_hainana)
clip_poly_temp <- matrix(c(57.52, 32.25,
                      63.3, 32.25,
                      63.3, 28.26,
                      57.52, 28.26,
                      57.52, 32.25),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_lutra <- range_temp - range_temp_try - clip_poly_temp
range_temp_lutra@data$subgroup <- "lutra"

range_temp_try <- raster::bind(range_temp_lutra, range_temp_try)
#map# tm_shape(range_temp_try) +tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "lutra"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Lycalopex culpaeus ####
# only one sample point :o
species_temp <- "Lycalopex culpaeus"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Lycalopex fulvipes ####
# only one sample
species_temp <- "Lycalopex fulvipes"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Lycalopex gymnocercus ####
# (wiki)
"Five subspecies are currently recognised, although the geographic range of each is unclear, and the type 
localities of three of them lie outside the present-day range of the species:[1][4]"
# Also has an apparently continuous range and subspecies range descriptions are not 

#map# species_temp <- "Lycalopex gymnocercus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Lycaon pictus ####
# five subspecies in wiki, though:
"these subspecific designations are not universally accepted"
"extensive intermixing has occurred between East African and Southern African populations in the past"
"Southern African and northeastern African populations, with a transition zone"
"West African wild dog population may possess a unique haplotype, thus possibly constituting a truly distinct subspecies"
"The West African wild dog used to be widespread from western to central Africa, from Senegal to Nigeria. "
# Seems as though there has been lots of mixing between the subspecies, apart from manguensis
# I at least have pictus and lupinus based on gbif, probably also somalicus
# Too hard to know what's what and where there is or has been contact between subpopulations

#map# species_temp <- "Lycaon pictus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Lynx canadensis ####
# cat group says "Therefore we conclude that Lynx canadensis is a monotypic species"
#map# species_temp <- "Lynx canadensis"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Lynx lynx ####
# On the basis of current evidence we propose the following six subspecies" (cat group)
# I only have ssp. lynx which is described as:
"Distribution: Scandinavia, Finland, Baltic States, Belarus, European part of Russia E to the Yenissei River."
# Also see their map
species_temp <- "Lynx lynx"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# using Yenisei to try and replicate cat map
# I am assuming that the Yenisei acts as a biogeographic boundary between subspecies
# Although the range does extend South beyond the extent of the river

Yeni_temp <- River_Data50[River_Data50$name == "Yenisey",]
Yeni_temp <- disaggregate(Yeni_temp)

Yeni_coords <- unlist(coordinates(Yeni_temp), recursive = FALSE)[[1]]
Yeni_coords <- Yeni_coords[nrow(Yeni_coords):1,]

Yeni_lynx_coords <- rbind(Yeni_coords,
                         matrix(c(Yeni_coords[nrow(Yeni_coords), 1], (range_temp@bbox[2,2] + 1),
                                  (range_temp@bbox[1,2] + 1), (range_temp@bbox[2,2] + 1),
                                  (range_temp@bbox[1,2] + 1), (range_temp@bbox[2,1] - 1),
                                  1.4, (range_temp@bbox[2,1] - 1),
                                  1.4, 56,
                                  38, 47,
                                  80, 49,
                                  86, 47.5,
                                  100, 42,
                                  109, 42,
                                  96, 48,
                                  Yeni_coords[1,1], Yeni_coords[1,2]), 
                                ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(Yeni_lynx_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string=CRS(proj4string(IUCN_Native_Data)))

#map# tm_shape(clip_poly_temp) + tm_polygons(alpha = 0) + tm_shape(range_temp) + tm_polygons(alpha = 0)

range_temp_otherssp <- raster::intersect(range_temp, clip_poly_temp)
range_temp_otherssp@data[is.na(range_temp_otherssp$subgroup), "subgroup"] <- "other"
range_temp_lynx <- range_temp - clip_poly_temp
range_temp_lynx@data$subgroup <- "lynx"

range_temp_try <- raster::bind(range_temp_otherssp, range_temp_lynx)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "lynx"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "^Yeni_"))

## Lynx pardinus ####
# This is a monotypic species - cat group
#map# species_temp <- "Lynx pardinus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Lynx rufus ####
# cat group splits into 2 on either side of the great plains which acts as a barrier
# They also group mexican populations into a rough subspecies
species_temp <- "Lynx rufus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# Going for quite a crude split trying to replicate cat map

buffer_temp <- countries50[countries50$name == "Mexico",]
buffer_temp <- terra::buffer(buffer_temp, 0.1)
# #map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(buffer_temp) + tm_polygons(alpha = 0)

range_temp_mexican <- raster::intersect(range_temp, buffer_temp)
range_temp_mexican@data$subgroup <- "mexican group"

range_temp_north <- range_temp - buffer_temp

clip_poly_temp <- matrix(c(-101, 29.5,
                      -102, 32,
                      -101.5, 35.5,
                      -100, 39,
                      -101.5, 45,
                      -106, 49,
                      -110, 51,
                      -58, 51,
                      -58, 23,
                      -101, 23,
                      -101, 29),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

clip_poly_temp <- clip_poly_temp - buffer_temp

range_temp_rufus <- raster::intersect(range_temp_north, clip_poly_temp)
range_temp_rufus@data$subgroup <- "rufus"

range_temp_fasciatus <- range_temp_north - clip_poly_temp
range_temp_fasciatus@data$subgroup <- "fasciatus"

range_temp_try <- raster::bind(range_temp_rufus, range_temp_mexican, range_temp_fasciatus)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# extending range_temp_try because of overlapping buffers
# floridanus not recognised
clip_points_temp <- matrix(c(-73, 53),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- SpatialPolygonsDataFrame(terra::buffer(clip_points_temp, 100000), 
                                        data = range_temp_rufus@data)

range_temp_try <- raster::bind(range_temp_try, clip_points_temp)

sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try,
                           buff = data.frame(subgroup = c("fasciatus", "rufus"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Martes americana ####
# thirteen subspecies noted by msotw but no range descriptions, and I have points across the range
# "Martes americana, which is now considered to comprise two subspecies-groups (americana and caurina)" - IUCN
# "found to intergrade in Montana and British Columbia" - IUCN
# Also don't know of any island or otherwise isolated subspecies
species_temp <- "Martes americana"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Martes foina ####
"Some of the island populations are morphologically quite distinct, although the taxonomic significance 
of this is not yet clear"  #IUCN
# based on wiki descriptions of eleven subspecies, I have:
# foina, mediterranea
# but not: bosniaca, bunites, intermedia, kozlovi, milleri, nehringi, rosanowi, syriaca, toufoeus

# bosniaca, nehringi, rosanowi, and syriaca are contiguous with foina 
# Unaware of physical barriers for bosniaca
# But nehringi separated by Caucasus, rosanowi isolated on Crimea, and Syriaca only contiguous via nehringi
# intermedia, bunites, kozlovi, milleri, toufoeus are geographically separate in IUCN polygons
# bunites and milleri are island populations

# bunites: only crete is included in polygon, not skopelos, erimomilos, karpathos, samothrace, seriphos, kythnos, or naxos

species_temp <- "Martes foina"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# bunites and milleri
range_temp_try <- raster::disaggregate(range_temp)
range_temp_try@data$subgroup <- "western"
range_temp_try@data <- range_temp_try@data %>%
  mutate(subgroup = case_when(island == "Crete" ~ "bunites",
                             island == "Rhodes" ~ "milleri",
                             TRUE ~ subgroup))

# central and east asian subspecies
clip_poly_temp <- matrix(c(41, 31,
                      92, 62,
                      121, 39, 
                      95, 15,
                      41, 31),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_eastern <- raster::intersect(range_temp_try, clip_poly_temp)
range_temp_eastern@data$subgroup <- "eastern"
range_temp_try <- range_temp_try - clip_poly_temp
range_temp_try <- raster::bind(range_temp_try, range_temp_eastern)

# nehringi and syriaca #
clip_poly_temp <- matrix(c(29, 40.8,
                      29, 41.03,
                      29.05, 41.06, 
                      29.065, 41.09,
                      29.08, 41.13,
                      29.07, 41.16,
                      29.16, 41.26,
                      35, 43,
                      40.01, 43.38,
                      40.01, 43.41, 
                      52, 39,
                      35, 30,
                      28.4, 36.5, 
                      27, 36.5,
                      25, 39.7,
                      26.3, 40.05,
                      26.4, 40.14,
                      26.4, 40.2,
                      26.8, 40.45,
                      29, 40.8),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- gUnion(clip_poly_temp, countries50[countries50$name %in% c("Georgia", "Azerbaijan"),])

range_temp_anatolia <- raster::intersect(range_temp_try, clip_poly_temp)
range_temp_anatolia@data$subgroup <- "nehringi syriaca"
range_temp_try <- range_temp_try - clip_poly_temp
range_temp_try <- raster::bind(range_temp_try, range_temp_anatolia)

# rosanowi
clip_poly_temp <- matrix(c(33.606, 46.132,
                      33.637, 46.140,
                      33.614, 46.226,
                      33.648, 46.229,
                      33.659, 46.223,
                      33.666, 46.222,
                      33.744, 46.184,
                      33.828, 46.206,
                      34.051, 46.107,
                      34.075, 46.118,
                      34.244, 46.051,
                      34.357, 46.063,
                      34.442, 45.964,
                      34.507, 45.944,
                      34.558, 45.992,
                      34.801, 45.898,
                      34.808, 45.798,
                      34.957, 45.757,
                      36.769, 46.181, 
                      36.659, 45.350,
                      36.478, 45.264,
                      36.516, 44.3,
                      32.3, 44.3,
                      32.3, 45.8,
                      33.56, 45.99,
                      33.606, 46.132),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_rosanowi <- raster::intersect(range_temp_try, clip_poly_temp)
range_temp_rosanowi@data$subgroup <- "rosanowi"
range_temp_try <- range_temp_try - clip_poly_temp
range_temp_try <- raster::bind(range_temp_try, range_temp_rosanowi)

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "western"

rm(list = ls(pattern = "^range_temp"))

## Martes martes #####
# has eight subspecies listed by msotw, but no range descriptions available that I can find
#map# species_temp <- "Martes martes"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Martes melampus ####
# I only have melampus
"The two confirmed subspecies of Japanese marten are:
M. m. melampus lives on several of the Japanese islands.
M. m. tsuensis is found on Tsushima Island.[3]" # wiki

species_temp <- "Martes melampus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("island") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == species_temp & (is.na(island) | island == "Kamishima") ~ "tsuensis",
                              binomial == species_temp ~ "melampus",
                              TRUE ~ subgroup))

#map# tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp,]) + tm_polygons("subgroup")

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "melampus"

## Martes pennanti ####
# "in general, the fisher is recognized to be a monotypic species with no extant subspecies.[11]" (wiki)
#map# species_temp <- "Martes pennanti"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Meles meles ####
# based on wiki description I have subspecies meles and marianensis
# not heptneri, though I can't remove the range of this one as its description overlaps with meles
# don't seem to have milleri, but uncertain of boundaries between meles and milleri
# Need to remove the range of M. canescens by splitting at the Caucasus mountains as it's now a separate species (msotw)

species_temp <- "Meles meles"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# removing M. canescens
clip_poly_temp <- matrix(c(29, 40.8,
                      29, 41.03,
                      29.05, 41.06, 
                      29.065, 41.09,
                      29.08, 41.13,
                      29.07, 41.16,
                      29.16, 41.26,
                      35, 43,
                      40.01, 43.38,
                      40.01, 43.41, 
                      52, 39, 
                      71, 39,
                      71, 29,
                      
                      27, 29,
                      27, 36.5,
                      25, 39.7,
                      26.3, 40.05,
                      26.4, 40.14,
                      26.4, 40.2,
                      26.8, 40.45,
                      29, 40.8),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

poly_temp <- terra::buffer(countries50[countries50$name %in% c("Georgia", "Azerbaijan"),], 0.15)
clip_poly_temp <- gUnion(clip_poly_temp, poly_temp)

range_temp_try <- range_temp - clip_poly_temp
#map# tm_shape(range_temp_try) + tm_polygons(alpha = 0) + tm_shape(range_temp) + tm_polygons(alpha = 0)

# marianensis
clip_poly_temp <- countries50[countries50$name %in% c("Spain", "Portugal"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.15)
range_temp_marianensis <- raster::intersect(range_temp_try, clip_poly_temp)
range_temp_marianensis@data$subgroup <- "marianensis"

# meles
range_temp_try <- range_temp_try - range_temp_marianensis
range_temp_try@data$subgroup <- "meles milleri heptneri"
range_temp_try <- raster::bind(range_temp_try, range_temp_marianensis)

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("marianensis", "meles milleri heptneri"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Melogale moschata/subaurantiaca ####
# Melogale subaurantiaca was previously considered subspecies of M. moschata
# now considered a full species occupying Taiwan
# All my samples appear to be M subaurantiaca so will adjust both polygons and gmpd

species_temp <- "Melogale moschata"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data <- GMPD_Spatial@data %>% 
  mutate(HostCorrectedName = case_when(HostCorrectedName == species_temp ~ "Melogale subaurantiaca",
                                       TRUE ~ HostCorrectedName)) %>%
  mutate(HostReportedSubspecies = case_when(HostCorrectedName == "Melogale subaurantiaca" ~ NA_character_,
                                            TRUE ~ HostReportedSubspecies))

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(binomial = case_when(binomial == species_temp & island == "Taiwan" ~ "Melogale subaurantiaca",
                              TRUE ~ binomial))

#map# species_temp <- "Melogale subaurantiaca"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Mephitis mephitis ####
# I have good sample coverage over most subspecies and they have a continuous range over North America
# Not aware of any island or otherwise isolated subspecies
#map# species_temp <- "Mephitis mephitis"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Mungos mungo ####
# only one sample location 
species_temp <- "Mungos mungo"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Mustela erminea ####
# Following wiki description of 21 subspecies ranges, where available
# North America: kadiacensis, arctica, polaris
# I only have arctica
# range description of arctica does not include US or Eastern Canada (polygons do), 
# but this area and range is contiguous with arctica anyway

# Europe: hibernica, ricinae, stabilis, minima, erminea, aestiva
# I have stabilis and aestiva
# IUCN polygon does not include Caucasus where teberdina is described
# same goes for ricina in Hebrides
# splitting at waterways at inland boundary of Kola peninsula for erminea,
# at Swiss borders for minima 
# (Alps on Southern border, Jura mountains on Northern border, and lakes Constance and Geneva at Northeast and Southwest)
# at the Urals for remaining unrepresented subspecies 

# Asia: tobolica, ferghanae, mongolica, lymani, nippon, karaginensis, kaneii
# Uncertain how to split south of the Yenisei, but none of these are represented anyway
species_temp <- "Mustela erminea"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# North America
# kadaiacensis
range_temp_kadiacensis <- range_temp[range_temp$island == "Kodiak" & !is.na(range_temp$island),]
range_temp_kadiacensis@data$subgroup <- "kadiacensis"

# polaris
range_temp_polaris <- range_temp[range_temp$island == "Greenland" & !is.na(range_temp$island),]
range_temp_polaris@data$subgroup <- "polaris"

# arctica
clip_poly_temp <- countries50[countries50$name %in% c("Canada", "United States"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.3)
range_temp_arctica <- raster::intersect(range_temp, clip_poly_temp)
range_temp_arctica <- range_temp_arctica[!range_temp_arctica$island %in% c("Kodiak", "Greenland"),]
range_temp_arctica@data$subgroup <- "arctica"

# Europe
# hibernica
poly_temp <- countries50[countries50$name == "United Kingdom",]
poly_temp <- raster::disaggregate(poly_temp)
clip_poly_temp <- countries50[countries50$name %in% c("Ireland", "Isle of Man"),]
clip_poly_temp <- raster::bind(clip_poly_temp, poly_temp[3,])
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.15)

range_temp_hibernica <- raster::intersect(range_temp, clip_poly_temp)
range_temp_hibernica@data$subgroup <- "hibernica"

# stabilis
clip_poly_temp <- poly_temp - clip_poly_temp
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.3)

clip_points_temp <- matrix(c(-1.63, 59.53),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- terra::buffer(clip_points_temp, 100000)
clip_poly_temp <- gUnion(clip_poly_temp, clip_points_temp)

range_temp_stabilis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_stabilis@data$subgroup <- "stabilis"

# aestiva
clip_poly_temp <- matrix(c(67.4, 68.8,
                      66.1, 68,
                      65.7, 67.2,
                      63.5, 66.5,
                      62.8, 65.8, 
                      61.1, 65,
                      59.9, 65.1,
                      58.9, 59.4,
                      59.7, 55.3, 
                      57.8, 54.5, 
                      57.1, 52.5,
                      57.2, 50.6,
                      51.3, 45.8,
                      -13.5, 39,
                      -13.5, 71.3, 
                      64.1, 71.3,
                      67.4, 68.8),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

clip_poly_temp <- clip_poly_temp - terra::buffer(countries50[countries50$name == "Kazakhstan",], 0.15) - range_temp_hibernica - range_temp_stabilis

range_temp_aestiva <- raster::intersect(range_temp, clip_poly_temp)
range_temp_aestiva@data$subgroup <- "aestiva"

# erminea
clip_poly_temp <- matrix(c(32.419, 67.130,
                      32.423, 67.132,
                      32.424, 67.138,
                      32.422, 67.143,
                      32.423, 67.148,
                      32.426, 67.149,
                      32.431, 67.155,
                      32.435, 67.164,
                      32.453, 67.168,
                      32.455, 67.173,
                      32.46, 67.175,
                      32.464, 67.18,
                      32.468, 67.181,
                      32.466, 67.183,
                      32.471, 67.184,
                      32.47, 67.188,
                      32.462, 67.193,
                      32.463, 67.198,
                      32.452, 67.208,
                      32.452, 67.211,
                      32.444, 67.218,
                      32.437, 67.221,
                      32.428, 67.233,
                      32.429, 67.236,
                      32.427, 67.238,
                      32.466, 67.276,
                      32.463, 67.28,
                      32.487, 67.29,
                      32.49, 67.293,
                      32.502, 67.299,
                      32.505, 67.304,
                      32.498, 67.309,
                      32.513, 67.316,
                      32.499, 67.355,
                      32.519, 67.369,
                      32.526, 67.381,
                      32.545, 67.388,
                      32.553, 67.392,
                      32.592, 67.404,
                      32.5, 67.5,
                      32.8, 67.54,
                      32.99, 67.57,
                      33.051, 67.603,
                      33.057, 67.611,
                      33.055, 67.626,
                      33.112, 67.662,
                      33.154, 67.699,
                      33.14, 67.741,
                      33.178, 67.766,
                      33.178, 67.864,
                      33.304, 67.962,
                      33.344, 67.995,
                      33.324, 68.013,
                      33.344, 68.043,
                      33.317, 68.072,
                      33.325, 68.08,
                      33.328, 68.09,
                      33.329, 68.091,
                      33.324, 68.094,
                      33.318, 68.094,
                      33.318, 68.095,
                      33.321, 68.098,
                      33.323, 68.099,
                      33.319, 68.104,
                      33.322, 68.106,
                      33.33, 68.109,
                      33.331, 68.109, 
                      33.334,68.107,
                      33.336, 68.106,
                      33.342, 68.109,
                      33.341, 68.112,
                      33.347, 68.119,
                      33.353, 68.113,
                      33.348, 68.21,
                      33.329, 68.218,
                      33.313, 68.221,
                      33.286, 68.221,
                      33.26, 68.214,
                      33.224, 68.261,
                      33.223, 68.325,
                      33.234, 68.331,
                      33.234, 68.336,
                      33.24, 68.34,
                      33.246, 68.342,
                      33.261, 68.351,
                      33.264, 68.352,
                      33.268, 68.353,
                      33.272, 68.353,
                      33.274, 68.352,
                      33.281, 68.354,
                      33.28, 68.355,
                      33.283, 68.356,
                      33.29, 68.357,
                      33.296, 68.362,
                      33.308, 68.365,
                      33.319, 68.375,
                      33.331, 68.433,
                      33.299, 68.443,
                      33.291, 68.451,
                      33.287, 68.452,
                      33.267, 68.45,
                      33.255, 68.451,
                      33.245, 68.449,
                      33.227, 68.448,
                      33.218, 68.449,
                      33.218, 68.45,
                      33.213, 68.45,
                      33.219, 68.454,
                      33.205, 68.461,
                      33.205, 68.463,
                      33.192, 68.471,
                      33.192, 68.473,
                      33.187, 68.474,
                      33.179, 68.481,
                      33.18, 68.482,
                      33.174, 68.485,
                      33.17, 68.486,
                      33.168, 68.49,
                      33.156, 68.494,
                      33.149, 68.509,
                      33.154, 68.509,
                      33.157, 68.511,
                      33.161, 68.512,
                      33.165, 68.515,
                      33.189, 68.516,
                      33.193, 68.522,
                      33.183, 68.525,
                      33.183, 68.526,
                      33.176, 68.529,
                      33.182, 68.531,
                      33.183, 68.536,
                      33.178, 68.541,
                      33.173, 68.544,
                      33.169, 68.55,
                      33.177, 68.557,
                      33.173, 68.559,
                      33.174, 68.562,
                      33.178, 68.564,
                      33.181, 68.57,
                      33.18, 68.572,
                      33.176, 68.572,
                      33.174, 68.574,
                      33.177, 68.577,
                      33.171, 68.578,
                      33.165, 68.578,
                      33.167, 68.581,
                      33.171, 68.582,
                      33.175, 68.586,
                      33.183, 68.589,
                      33.186, 68.591,
                      33.191, 68.593,
                      33.198, 68.598,
                      33.205, 68.602,
                      33.214, 68.604,
                      33.215, 68.606,
                      33.222, 68.607,
                      33.232, 68.613,
                      33.233, 68.616,
                      33.238, 68.62,
                      33.234, 68.623,
                      33.236, 68.625,
                      33.232, 68.63,
                      33.23, 68.63,
                      33.225, 68.632,
                      33.226, 68.635,
                      33.221, 68.639,
                      33.217, 68.64,
                      33.214, 68.643,
                      33.208, 68.643,
                      33.202, 68.646,
                      33.196, 68.647,
                      33.192, 68.65,
                      33.182, 68.654,
                      33.16, 68.658,
                      33.15, 68.666,
                      33.136, 68.674,
                      33.136, 68.678,
                      33.138, 68.679,
                      33.132, 68.697,
                      33.127, 68.698,
                      33.117, 68.707,
                      33.116, 68.71,
                      33.113, 68.712,
                      33.107, 68.715,
                      33.105, 68.718,
                      33.1, 68.72, 
                      33.097, 68.722,
                      33.098, 68.727,
                      33.096, 68.734,
                      33.102, 68.736,
                      33.116, 68.735,
                      33.124, 68.736,
                      33.135, 68.737,
                      33.145, 68.74,
                      33.153, 68.745,
                      33.152, 68.749,
                      33.156, 68.75,
                      33.161, 68.75,
                      33.164, 68.754,
                      33.157, 68.762,
                      33.152, 68.763,
                      33.151, 68.766,
                      33.147, 68.77,
                      33.146, 68.775,
                      33.147, 68.776, 
                      33.145, 68.782, 
                      33.128, 68.786,
                      33.115, 68.786,
                      33.108, 68.789,
                      33.1, 68.798,
                      33.088, 68.802,
                      33.097, 68.809,
                      33.095, 68.813,
                      33.078, 68.819,
                      33.078, 68.823,
                      33.083, 68.83,
                      33.092, 68.835,
                      33.081, 68.844,
                      33.068, 68.849,
                      33.046, 68.848,
                      33.028, 68.851,
                      33.024, 68.856,
                      33.033, 68.86,
                      33.038, 68.864,
                      33.038, 68.866,
                      33.03, 68.868,
                      33.03, 68.871,
                      33.034, 68.874,
                      33.031, 68.879,
                      33.035, 68.885,
                      33.028, 68.888,
                      33.043, 68.906,
                      33.015, 68.946,
                      33.023, 68.96,
                      33.045, 68.977,
                      33.036, 69.009,
                      33.067, 69.056,
                      33.153, 69.067,
                      33.285, 69.079,
                      33.406, 69.098,
                      33.454, 69.134,
                      33.53, 69.161,
                      33.539, 69.251,
                      33.516, 69.306,
                      33.7, 69.5,
                      42, 69,
                      42, 67,
                      40.8, 66.2,
                      38.5, 65.7,
                      34.2, 66.4,
                      33.3, 66.72,
                      32.65, 67.02,
                      32.386, 67.117,
                      32.419, 67.130),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_erminea <- raster::intersect(range_temp_aestiva, clip_poly_temp)
range_temp_erminea@data$subgroup <- "erminea"
range_temp_aestiva <- range_temp_aestiva - clip_poly_temp
range_temp_aestiva <- raster::bind(range_temp_aestiva, range_temp_erminea)

# minima
clip_poly_temp <- countries50[countries50$name == "Switzerland",]
range_temp_minima <- raster::intersect(range_temp_aestiva, clip_poly_temp)
range_temp_minima@data$subgroup <- "minima"
range_temp_minima@data <- range_temp_minima@data[,1:29]

range_temp_aestiva <- range_temp_aestiva - clip_poly_temp
range_temp_aestiva <- raster::bind(range_temp_aestiva, range_temp_minima)

#map# tm_shape(range_temp_aestiva) + tm_polygons("subgroup")

# remaining subspecies
range_temp_try <- raster::bind(range_temp_kadiacensis, range_temp_polaris, range_temp_arctica, range_temp_stabilis, range_temp_hibernica, range_temp_aestiva)
range_temp_otherssp <- range_temp - range_temp_try
range_temp_otherssp@data$subgroup <- "kaneii and other"
range_temp_try <- raster::bind(range_temp_try, range_temp_otherssp)

#map# tm_shape(range_temp_try) +tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("arctica", "stabilis", "aestiva"), buff = c(0, 0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Mustela lutreola ####
# seven subspecies, I only seem to have the french mink (biedermanni)
# uncertain of where to split subspecies from descriptions on wiki
#map# species_temp <- "Mustela lutreola"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Mustela nivalis ####
# Basing division on geography, rather than morphological size based categories (from wiki)
# North American subspecies: allegheniensis, campestris, eskimo, rixosa
# European/Asian: nivalis, boccamela, caucasica, heptneri, namiyei, numidica, pallida, pygmaea, russelliana, vulgaris
species_temp <- "Mustela nivalis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

clip_poly_temp <- countries50[countries50$continent == "North America",]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.8)

range_temp_americas <- raster::intersect(range_temp, clip_poly_temp)
range_temp_americas@data$subgroup <- "american"

range_temp_eurasian <- range_temp - clip_poly_temp
range_temp_eurasian@data$subgroup <- "eurasian"

# numidica
range_temp_numidica <- raster::disaggregate(range_temp_eurasian)
clip_points_temp <- matrix(c(-3, 34),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_numidica <- range_temp_numidica[clip_points_temp, ]
range_temp_numidica@data$subgroup <- "numidica"

# namiyei
range_temp_eurasian@data[range_temp_eurasian$island == "Honshu" & !is.na(range_temp_eurasian$island), "subgroup"] <- "namiyei"

# formosana (taiwan)
# http://dx.doi.org/10.3106/041.035.0305
range_temp_formosana <- raster::disaggregate(range_temp_eurasian)
clip_points_temp <- matrix(c(121, 23.5),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_formosana <- range_temp_formosana[clip_points_temp, ]
range_temp_formosana@data$subgroup <- "formosana"

range_temp_eurasian <- range_temp_eurasian - range_temp_numidica - range_temp_formosana

range_temp_try <- raster::bind(range_temp_americas, range_temp_eurasian, range_temp_numidica, range_temp_formosana)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "eurasian"

rm(list = ls(pattern = "^range_temp"))

## Mustela putorius ####
# msotw has seven subspecies, and based on wiki ranges I have putorius and furo
# not aureola, or mosquensis but these are contiguous with putorius so I won't remove
# not anglia or caledoniae and I can remove as they are isolated
# not rothschildi which is restricted to Dobruja region, inferred to be isolated by Danube and Balkan mountains

species_temp <- "Mustela putorius"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

range_temp_try <- raster::disaggregate(range_temp)

# UK
buffer_temp <- countries50[countries50$name == "United Kingdom",]
buffer_temp <- terra::buffer(buffer_temp, 0.3)

range_temp_ukssp <- raster::intersect(range_temp_try, buffer_temp)
range_temp_ukssp@data$subgroup <- "anglia caledoniae"

# Europe: putorius, aureola, mosquensis, rothschildi, furo
range_temp_putorius <- range_temp_try - buffer_temp
range_temp_putorius@data$subgroup <- "putorius and european"

# rothschildi
DB_temp <- River_Data50[River_Data50$name %in% c("Danube", "Bratul Chillia"),]

DB_temp <- disaggregate(DB_temp)
DB_temp$ID <- LETTERS[1:nrow(DB_temp)]
#map# tm_shape(DB_temp) + tm_lines("ID", lwd = 2)

DB_coords <- rbind(DB_temp@lines[[2]]@Lines[[1]]@coords[116:144,], 
                   DB_temp@lines[[1]]@Lines[[1]]@coords[14:1,])

DB_roth_coords <- rbind(DB_coords,
                        matrix(c(29.78, 45.14,
                                 29.78, 44.61,
                                 28.76, 42.77,
                                 27.874, 42.837,
                                 #28.62, 43.26,
                                 DB_coords[1,]), 
                               ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(DB_roth_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

poly_temp <- matrix(c(29.78, 45.14,
                      29.78, 44.61, 
                      28.76, 42.77,
                      27.874, 42.837,
                      27.687, 42.881,
                      27.484, 42.881,
                      27.27, 42.981,
                      26.564, 42.905,
                      26.366, 42.906,
                      26.004, 42.791,
                      25.631, 42.766,
                      25.059, 42.758,
                      24.125, 42.784,
                      23.68, 42.857,
                      23.461, 43.151,
                      23.109, 43.114,
                      23.015, 43.198,
                      22.3, 43.18,
                      22.3, 44.248,
                      28.2, 45.64, 
                      29.5, 45.5,
                      29.88, 45.3), 
                    ncol = 2, byrow = TRUE)

poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
poly_temp <- raster::intersect(poly_temp, countries50[countries50$name == "Bulgaria",])

clip_poly_temp <- gUnion(poly_temp, clip_poly_temp)
range_temp_rothschildi <- raster::intersect(range_temp_putorius, clip_poly_temp)
range_temp_rothschildi@data$subgroup <- "rothschildi"

range_temp_putorius <- range_temp_putorius - range_temp_rothschildi
range_temp_putorius <- raster::bind(range_temp_putorius, range_temp_rothschildi)

range_temp_try <- raster::bind(range_temp_putorius, range_temp_ukssp)
range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
#map# tm_shape(range_temp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
# one sample location assigned originally assigned furo is now in european group
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "putorius and european"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))
rm(list = ls(pattern = "^DB_"))

## Nasua nasua ####
# only one sample location
species_temp <- "Nasua nasua"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Neovison vison ####
# 15 subspecies in msotw, with range descriptions from wiki
# I have multiple subspecies represented: vison, energumenos, evergladensis, and others 
# only one that I would be confident in excluding is nesolestes but Admiralty Island is not included in the range anyway
species_temp <- "Neovison vison"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Nyctereutes procyonoides ####
# Going to follow current IUCN treatment of japanese populations as subspecies 
# I don't have: koreensis, orestes, procyonoides, ussuriensis
"There are six recognized subspecies: albus, koreensis, orestes, procyonoides, ussuriensis and viverrinus. 
It has been suggested recently that Japanese Raccoon Dogs should be classified as a distinct species Nyctereutes 
viverrinus with two subspecies N. v. viverrinus and N. v. albus (Sang-In et al. 2015). The classification is based 
on chromosomal, molecular and morphological differences between Japanese and mainland populations (see also Kauhala 
and Saeki 2004a)." #IUCN

"The raccoon dogs from Hokkaido are sometimes recognized as a different subspecies from the mainland tanuki as 
Nyctereutes procyonoides albus (Hornaday, 1904) (or N. viverrinus albus if recognized as a distinct species). 
This taxon is synonymized with N. p. viverrinus in Mammal Species of the World,[13][16] but comparative morphometric 
analysis supports recognizing the Hokkaido population as a distinct subspecific unit.[13]" #wiki

species_temp <- "Nyctereutes procyonoides"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("island") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == species_temp & island == "Hokkaido" ~ "albus",
                              binomial == species_temp & str_detect(island, "Honshu|Kyushu|Sado|Shikoku") ~ "viverrinus",
                              binomial == species_temp & island == "Sakhalin" ~ "other",
                              binomial == species_temp & is.na(island) ~ "other",
                              TRUE ~ subgroup))

#map# tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp,]) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# GMPD subgroup assignment
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("albus", "viverrinus"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Odocoileus hemionus ####
# 10 subspecies according to msotw, also recognised by IUCN
# Following wiki description of ranges: 
# I have californicus, columbianus, crooki, hemionus, inyoensis
# Don't have cerrosensis, fuliginatus, peninsulae, sheldoni, sitkensis
# Don't need to exclude cerrosensis or sheldoni as Cedros Island and Tiburon island not included in polygon anyway

# the IUCN range polygon differs from the ranges reported on wiki
# if wiki is correct, the remaining subspecies not represented in gmpd are contiguous with those that are 
# if the IUCN range is correct, peninsulae subspecies is isolated

# Going to leave it be as gap in IUCN may be artefact and I have no other evidence of isolation from fuliginatus

species_temp <- "Odocoileus hemionus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Odocoileus virginianus ####
# following wiki images of subspecies distributions
# Most of North American range is continuous and generally well sampled, 
# and subspecies ranges do not match closely enough geographic features for me to be confident in splitting
# North american island subspecies: hiltonensis, mcilhennyi, rothschildi, taurinsulae, venatorius
# Some of these small islands are not represented in IUCN range polygons, but are inconsequential to analysis anyway

# Gap in wiki polygons is at Darien gap so I will split there
# https://commons.wikimedia.org/wiki/File:Odocoileus_virginianus_SA_map.svg

species_temp <- "Odocoileus virginianus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

buffer_temp <- countries50[countries50$continent == "North America",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with

clip_points_temp <- matrix(c(-96.3, 28.4,
                        -87.6, 13.3,
                        -79.7, 8.6,
                        -89.9, 29.3,
                        -76.3, 35.2,
                        -76.2, 36.1,
                        -76.4, 38.9,
                        -69.9, 41.5,
                        -70.1, 43.7,
                        -64.4, 44.2,
                        -66.4, 50.2),
                      ncol = 2, byrow = TRUE) # patching up gaps
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
buffer_temp <- gUnion(buffer_temp, terra::buffer(clip_points_temp, 50000))

range_temp_north <- raster::intersect(range_temp, buffer_temp)
range_temp_north@data$subgroup <- "northern"
range_temp_south <- range_temp - buffer_temp
range_temp_south@data$subgroup <- "southern"

range_temp_try <- raster::bind(range_temp_north, range_temp_south)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
# texanus sample overwritten
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "northern"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Oreamnos americanus ####
# no subspecies on iucn, wiki, or msotw
species_temp <- "Oreamnos americanus"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Otocolobus manul ####
# only one sample point :(
species_temp <- "Otocolobus manul"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Otocyon megalotis ####
# already got the subspecies <3
species_temp <- "Otocyon megalotis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("megalotis", "virgatus"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Ourebia ourebi ####
# only extant isolated subspecies is already labelled, so referring to main population as nominate subspecies:
"Numerous subspecies of the Oribi have been described but most of these reflect individual variation and have little 
or no validity. Haggard's Oribi (O. o. haggardi) of eastern coastal Kenya and adjacent Somalia is a geographically 
isolated subspecies which is well differentiated in size and colour from other Oribi. Another distinctive subspecies 
from East Africa, the Kenya Oribi (O. o. keniae) from the lower slopes of Mount Kenya, is now apparently extinct."

species_temp <- "Ourebia ourebi"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == species_temp & is.na(subgroup) ~ "oribi",
                              TRUE ~ subgroup)) 
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "oribi"

#map# tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp,]) + tm_polygons("subgroup") + 
  #map# tm_shape(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]) + tm_dots("subgroup") 

## Ovibos moschatus ####
# no subspecies on iucn, wiki, msotw
species_temp <- "Ovibos moschatus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Ovis ammon ####
# Sample locations don't correspond at all to Ovis ammon range
# see "Subspecies info" as well 
species_temp <- "Ovis ammon"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Ovis canadensis ####
# Not following the IUCN usage of 7 subspecies. following wiki instead as it's more recent
# I have: canadensis, sierrae, and nelsoni

"The 2016 genetics study suggested more modest divergence of this desert bighorn sheep into three lineages 
consistent with the earlier work of Cowan: Nelson's (O. c. nelsoni), Mexican (O. c. mexicana), and Peninsular 
(O. c. cremnobates). These three lineages occupy desert biomes that vary significantly in climate, suggesting 
exposure to different selection regimens."

# unable to split according to lineages or subspecies as their range descriptions are not specific enough

species_temp <- "Ovis canadensis"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Ovis dalli ####
# IUCN states there are two subspecies, but they are known to admix
# see Fannin sheep (O. d. fannini)
species_temp <- "Ovis dalli"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- NA

## Ozotoceros bezoarticus ####
# Not entirely sure which subspecies I have (only leucogaster, or bezoarticus as well) but can remove those I definitely don't have
# IUCN polygons do not match well with wiki description of ranges (wiki descriptions are more expansive)
"O. b. bezoarticus (Linnaeus, 1758), in the Cerrado ranging from eastern and central Brazil, south of the Amazon river 
between the plateau of Mato Grosso and the upper San Francisco river. While there is a population further north on the 
island of Marajo, at the mouth of the Amazon River (Rossetti and de Toledo 2006).

O. b. celer Cabrera, 1943, inhabiting the entire Argentinean pampas from the Atlantic coast to Sub Andean foothills and 
southward to the Rio Negro basin.

O. b. leucogaster (Goldfüss, 1817), located in the seasonally flooded grasslands of southwestern Brazil in southern 
Mato Grosso, southeastern Bolivia, Paraguay and in the Chaco savannas in northern Argentina (northern Santiago del 
Estero, Santa Fe, Formosa and Corrientes).

O. b. arerunguaensis González, Álvares-Valin and Maldonado, 2002. Type locality - Uruguay; northwestern Salto 
Department, Arerunguá, El Tapado, 31º41’51”S, 56º43’31”W.

O. b. uruguayensis González, Álvares-Valin and Maldonado, 2002. Type locality - Uruguay, grasslands of eastern Rocha 
Department, Sierra de Los Ajos, 33º50’01”S; 54º01’34”W."

species_temp <- "Ozotoceros bezoarticus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# celer
clip_points_temp <- matrix(c(-62.8, -38.8,
                        -66.3, -34.8,
                        -57.2, -36.7),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- raster::disaggregate(range_temp)
range_temp_celer <- range_temp_try[clip_points_temp,]
range_temp_celer@data$subgroup <- "celer"

# uruguayensis
clip_points_temp <- matrix(c(-54, -33.5),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_uruguayensis <- range_temp_try[clip_points_temp,]
range_temp_uruguayensis@data$subgroup <- "uruguayensis"

# arerunguaensis
clip_points_temp <- matrix(c(-56.7, -31.3),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_arerunguaensis <- range_temp_try[clip_points_temp,]
range_temp_arerunguaensis@data$subgroup <- "arerunguaensis"

# leucogaster
clip_points_temp <- matrix(c(-60.8, -28.2,
                        -56.2, -28,
                        -51, -27.8,
                        -57.4, -22.2,
                        -58.5, -13.4,
                        -62, -16.1,
                        -60.9, -14.5,
                        -63.2, -13,
                        -68.3, -12.6,
                        -58.5, -17),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_leucogaster <- range_temp_try[clip_points_temp,]
range_temp_leucogaster@data$subgroup <- "leucogaster"

# bezoarticus
clip_points_temp <- matrix(c(-51.8, -22.3,
                        -45.9, -21,
                        -52.2, -16.4,
                        -45.5, -12.9,
                        -52, -12.1, 
                        -49.6, -9.5,
                        -50.8, -6.5),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_bezoarticus <- range_temp_try[clip_points_temp,]
range_temp_bezoarticus@data$subgroup <- "bezoarticus"

range_temp_try <- raster::bind(range_temp_celer, range_temp_uruguayensis, range_temp_arerunguaensis, range_temp_leucogaster, range_temp_bezoarticus)

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("leucogaster", "bezoarticus"), buff = c(0.5, 1.5)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))

## Paguma larvata ####
# I only have samples in their introduced range
species_temp <- "Paguma larvata"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Panthera leo ####
# following the now outdated subspecies taxonomy of Asia (persica) and Africa (leo) as Asian population
# is so geographically isolated
species_temp <- "Panthera leo"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == species_temp & is.na(dist_comm) ~ "leo",
                              binomial == species_temp ~ "persica",
                              TRUE ~ subgroup))

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "leo"

#map# tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp,]) + tm_polygons("subgroup") + 
  #map# tm_shape(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]) + tm_dots() 

rm(list = ls(pattern = "^range_temp"))

## Panthera onca ####
# Cat group describes it as a monotypic species, 
# but there are sufficiently isolated subpopulations that I'm going to split
# Can't split South American subgroups because of lack of specificity over which Amazon tributaries
# Also does not appear to be as isolated as Northern and Central American subgroups
"The status of the subspecies is unclear. Although eight subspecies have been recognized (Seymour 1989), 
morphological and genetic analyses do not support the existence of discrete subspecies (Larson 1997, Eizirik 
et al. 2001, Ruiz-Garcia et al. 2006). While not elevating the regional differences to the subspecies level, 
Eizirk et al. (2001) found evidence for four incompletely isolated phylogeographic groups: Mexico and Guatemala,
southern Central America, northern South America, and South America south of the Amazon River. Similarly, 
Ruiz-Garcia et al. (2006) found that the Andes Mountains incompletely isolates Jaguar populations in Colombia."

species_temp <- "Panthera onca"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# Mexico and Guatemala
clip_poly_temp <- countries50[countries50$name %in% c("Mexico", "Guatemala", "Belize"), ]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.3)
range_temp_northern <- raster::intersect(range_temp, clip_poly_temp)
range_temp_northern@data$subgroup <- "northern"

# southern Central America
clip_poly_temp <- matrix(c(-71.3, 9.3,
                      -72.1, 8.9,
                      -72.1, 8.1,
                      -72.7, 7.5,
                      -72.7, 6,
                      -78.5, 0,
                      -80, -3,
                      -88, -3,
                      -88, 17,
                      -71.3, 17,
                      -71.3, 9.3),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_central <- raster::intersect(range_temp, clip_poly_temp)
range_temp_central@data$subgroup <- "central"

# Southern
range_temp_try <- terra::buffer(raster::bind(range_temp_central, range_temp_northern), 0.3)
range_temp_south <- range_temp - range_temp_try
range_temp_south@data$subgroup <- "southern"

range_temp_try <- raster::bind(range_temp_south, range_temp_central, range_temp_northern)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "southern"

rm(list = ls(pattern = "^range_temp"))

## Panthera pardus ####
# already has subspecies designations <3
species_temp <- "Panthera pardus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "pardus"

## Pecari tajacu ####
# splitting at Panama border
"DNA studies suggest that P. tajacu may consist of at least two major clades or lineages comprising specimens 
from North/Central and South America (Gongora et al. 2006, 2011) with structural chromosomal differences 
(Gongora et al. 2000, Adega et al. 2006). Additional studies are needed to clarify whether intra and/or 
inter-specific genetic and chromosomal variation is occurring within this species."

species_temp <- "Pecari tajacu"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

buffer_temp <- countries50[countries50$continent == "North America" & countries50$name != "Trinidad and Tobago",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with

clip_points_temp <- matrix(c(-96.3, 28.4,
                        -87.6, 13.3,
                        -79.7, 8.6),
                      ncol = 2, byrow = TRUE) # patching up gaps
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
buffer_temp <- gUnion(buffer_temp, terra::buffer(clip_points_temp, 10000))

range_temp_north <- raster::intersect(range_temp, buffer_temp)
range_temp_north@data$subgroup <- "northern"
range_temp_south <- range_temp - buffer_temp
range_temp_south@data$subgroup <- "southern"

range_temp_try <- raster::bind(range_temp_north, range_temp_south)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "northern"

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Pelea capreolus ####
# continuous range with samples across it, also no subspecies described by iucn, wiki, msotw
#map# species_temp <- "Pelea capreolus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Phacochoerus aethiopicus ####
# range only includes extant subspecies, delamerei
"the Cape Warthog, Phacochoerus aethiopicus aethiopicus, endemic to South Africa, became extinct in the 1870s"
"the Somali Warthog, P. a. delamerei, occurs in Kenya and the Horn of Africa."

species_temp <- "Phacochoerus aethiopicus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp, bbox = bbox(sp.gmpd.points)) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# Some suspicious GMPD samples
# odd ones are probably introduced for hunting

sp.gmpd.points <- sp.gmpd.points[terra::buffer(range_temp, 4),]
GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

"See Grubb and d’Huart (2010) for a detailed historic overview of the classification of Phacochoerus."
# https://doi.org/10.2982/028.099.0204

## Philantomba monticola ####
# lots of subspecies, most of which have contiguous ranges, broadly split into two groups which overlap geographically
# I seem to have bicolor which has a range description as follows:
# "The range extends from Zanzibar to the KwaZulu Natal region in South Africa."
# This description includes Northern and Southern polygons so will leave the range as is 
#map# species_temp <- "Philantomba monticola"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Procapra gutturosa ####
"Procapra gutturosa has not given rise to distinct geographic races. Specimens from the Mongolian Altai are
indistinguishable from those of eastern Mongolia and hence the subspecies P. g. altaica described from 
Bayan-Tsagan-Gobi is not accepted (Sokolov and Lushchekina 1997)." #IUCN

#map# species_temp <- "Procapra gutturosa"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

## Procyon lotor ####
# Many subspecies. Based on wiki I appear to have:
# lotor, elucus, fuscipes, hirtus, litoreus, marinus, megalodous, pacificus, psora, simus
# All of these are mainland subspecies except:
# litoreus: "Coastal strip and islands of Georgia."
# marinus: "Keys of the Ten Thousand Islands Group, and adjoining mainland of southwestern Florida from Cape Sable north through the Everglades to Lake Okeechobee."
## Both of these are hard to separate and not completely isolated from mainland subspecies, so will include with mainland

# Don't appear to have
# excelsus (contiguous with others), auspicatus (polygon not included), gloveralleni (extinct and not included),
# grinnelli (isolated), hernandezii (contiguous), incautus (polygon not included), inesperatus (not fully isolated), 
# insularis (isolated), maynardi (polygon not included), minor (polygon not included), pallidus (continguous), 
# pumilus (contiguous), vancouverensis (isolated but polygon not included)
## only separating fully isolated subspecies (grinnelli and insularis)

species_temp <- "Procyon lotor"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# grinnelli
clip_poly_temp <- matrix(c(-112.1, 27.9,
                      -107.7, 23,
                      -110.6, 21.5,
                      -115.3, 27.2,
                      -112.1, 27.9),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_grinnelli <- raster::intersect(range_temp, clip_poly_temp)
range_temp_grinnelli@data$subgroup <- "grinnelli"

# insularis
clip_poly_temp <- matrix(c(-106.3, 22.5,
                      -105.9, 21.5,
                      -106.3, 20.5,
                      -107.3, 21.5,
                      -106.3, 22.5),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_insularis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_insularis@data$subgroup <- "insularis"

# inesperatus
clip_poly_temp <- matrix(c(-80.46361,25.18777,
                      -80.40206,25.21110,
                      -80.36805,25.28172,
                      -80.26471,25.37177,
                      -80.21612,25.33873,
                      -80.35152,25.11391,
                      -80.58412,24.92214,
                      -80.62105,24.96879),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_inesperatus <- raster::intersect(range_temp, clip_poly_temp)
range_temp_inesperatus@data$subgroup <- "inesperatus"

# mainland subspecies
range_temp_mainland <- range_temp - range_temp_grinnelli - range_temp_insularis - range_temp_inesperatus
range_temp_mainland@data$subgroup <- "mainland"

range_temp_try <- raster::bind(range_temp_mainland, range_temp_grinnelli, range_temp_insularis, range_temp_inesperatus)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "mainland"

rm(list = ls(pattern = "^range_temp"))

## Procyon pygmaeus ####
# only one sample 
species_temp <- "Procyon pygmaeus"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Puma concolor ####
# cat group recognises two subspecies: 
# concolor "South America, possibly excluding W of Andes in north."
# couguar "North and Central America, possibly N South America W of Andes."
# Also going to remove Florida population as it's isolated and has experienced severe inbreeding depression

species_temp <- "Puma concolor"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# couguar
clip_poly_temp <- matrix(c(-77.2, 9.8,
                      -81, 6.1,
                      -108.1, 17.6,
                      -127.9, 39.7,
                      -132.4, 54,
                      -128.9, 57.7,
                      -123.5, 59.7,
                      -102.1, 60,
                      -86, 22.6,
                      -77.2, 9.8),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- raster::disaggregate(range_temp)
range_temp_couguar <- range_temp_try[clip_poly_temp,]
range_temp_couguar <- raster::aggregate(range_temp_couguar, by = names(range_temp_couguar))
range_temp_couguar@data$subgroup <- "couguar"

# concolor
florida_temp <- terra::buffer(states50[states50$name == "Florida",], 0.1)
range_temp_concolor <- range_temp - range_temp_couguar - florida_temp - clip_poly_temp
range_temp_concolor@data$subgroup <- "concolor"

range_temp_try <- raster::bind(range_temp_couguar, range_temp_concolor)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("concolor", "couguar"), buff = c(0, 2.5)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Rangifer tarandus ####
# 14 subspecies listed by msotw
# North America:
# I may have: caboti, caribou, granti, groenlandicus, pearyi, terraenovae
# I don't have: osborni (contiguous with others)

# Russia & Finland:
# I don't have: buskensis, fennicus, pearsoni, phylarchus, sibiricus, valentinae

# Europe:
# Have: platyrhynchus in Svalbard, tarandus in Norway

# Based on geographic overlap between subspecies, I'm grouping them as follows:
# all north american and Greenland subspecies = american
# Russia and Finland 
# Norway = tarandus
# Svalbard = platyrhynchus

# Excluding outlying points in North Norway because of overlap with semi-domestic herds (see Sami wiki page)
# Although caribou in North America are hunted, they remain undomesticated

species_temp <- "Rangifer tarandus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# american subspecies
buffer_temp <- countries50[countries50$continent == "North America",]
buffer_temp <- terra::buffer(buffer_temp, 0.3) 

range_temp_americas <- raster::intersect(range_temp, buffer_temp)
range_temp_americas@data$subgroup <- "american"

# platyrhynchus
norway_temp <- countries50[countries50$name == "Norway",]
norway_temp <- terra::buffer(norway_temp, 0.3) 

range_temp_tara_platy <- raster::intersect(range_temp, norway_temp)
range_temp_tara_platy@data$subgroup <- "platyrhynchus"

# tarandus
range_temp_tara_platy <- raster::disaggregate(range_temp_tara_platy)
range_temp_tara_platy@data[2, "subgroup"] <- "tarandus"

range_temp_otherssp <- range_temp - buffer_temp - norway_temp
range_temp_otherssp@data$subgroup <- "other"
range_temp_try <- raster::bind(range_temp_americas, range_temp_tara_platy, range_temp_otherssp)

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("american", "tarandus", "platyrhynchus"), buff = c(0, 1, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Raphicerus campestris ####
# already done <3
species_temp <- "Raphicerus campestris"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "campestris"

## Redunca arundinum ####
# iucn, wiki, msotw list no subspecies
#map# species_temp <- "Redunca arundinum"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Redunca fulvorufula ####
# already done <£
species_temp <- "Redunca fulvorufula"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "fulvorufula"

## Rupicapra pyrenaica ####
# partially done, and wiki states:
"R. p. pyrenaica (Pyrenean chamois): France and Spain
R. p. parva (Cantabrian chamois): Spain
R. p. ornata (Abruzzo chamois): Central and southern Italy"
# this agrees with IUCN geographic range description 

species_temp <- "Rupicapra pyrenaica"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

range_temp_ornata <- raster::intersect(range_temp, terra::buffer(countries50[countries50$name == "Italy",], 0))
range_temp_ornata@data$subgroup <- "ornata"

range_temp_parva <- raster::intersect(range_temp[is.na(range_temp$subspecies),], terra::buffer(countries50[countries50$name == "Spain",], 0))
range_temp_parva@data$subgroup <- "parva"

range_temp_try <- raster::bind(range_temp_ornata,
                         range_temp_parva,
                         range_temp[!is.na(range_temp$subspecies),])

range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("parva", "pyrenaica"), buff = c(0.5, 1)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))

## Rupicapra rupicapra ####
"The species R. rupicapra is categorized into seven subspecies:

R. r. asiatica (Anatolian chamois or Turkish chamois): Turkey
R. r. balcanica (Balkan chamois): Albania, Bosnia and Herzegovina, Bulgaria, Croatia, northern Greece 
(the Pindus Mountains),[4] North Macedonia, Serbia, Montenegro, and Slovenia (isolated populations)
R. r. carpatica (Carpathian chamois): Romania
R. r. cartusiana (Chartreuse chamois): France
R. r. caucasica (Caucasian chamois): Azerbaijan, Georgia, Russia
R. r. rupicapra (Alpine chamois): Austria, France, Germany, Italy, Switzerland, Slovenia,[5] Slovakia (Veľká Fatra, Slovak Paradise)
R. r. tatrica (Tatra chamois): Slovakia (Tatras and Low Tatras) and Poland (Tatras)"

# asiatica, caucasica, carpatica, and tatrica subspecies are isolated and not represented
# Don't have cartusiana based on range description here: http://www.wilddocu.de/chartreuse-chamois-rupicapra-rupicapra-cartusiana/
# Can't separate it anyway because of gbif data distribution
# balcanica is a bit hard to separate, but need to do it anyway

species_temp <- "Rupicapra rupicapra"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# tatrica
tatrica_temp <- terra::buffer(countries50[countries50$name == "Slovakia",], 0.1)
range_temp_tatrica <- raster::intersect(range_temp, tatrica_temp)
range_temp_tatrica@data$subgroup <- "tatrica"

# carpatica
carpatica_temp <- terra::buffer(countries50[countries50$name == "Romania",], 0)
range_temp_carpatica <- raster::intersect(range_temp, carpatica_temp)
range_temp_carpatica@data$subgroup <- "carpatica"

# asiatica
asiatica_temp <- terra::buffer(turkey_rangepol[2,], 1)
range_temp_asiatica <- raster::intersect(range_temp, asiatica_temp)
range_temp_asiatica@data$subgroup <- "asiatica"

# caucasica
caucasica_temp <- terra::buffer(countries50[countries50$name == "Georgia",], 1.5)
caucasica_temp <- caucasica_temp - asiatica_temp
range_temp_caucasica <- raster::intersect(range_temp, caucasica_temp)
range_temp_caucasica@data$subgroup <- "caucasica"

# balcanica
balcanica_temp <- countries50[countries50$name %in% c("Greece", "Bulgaria", "Macedonia", "Albania", 
                                                      "Montenegro", "Serbia", "Bosnia and Herz.", "Kosovo"),]
balcanica_temp <- terra::buffer(raster::aggregate(balcanica_temp), 0.25)

clip_poly_temp <- matrix(c(17, 44,
                      14.9, 44,
                      14.9, 45.3,
                      15.4, 45.5,
                      14.9, 45.8,
                      16.5, 46.5,
                      17, 44),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

balcanica_temp <- gUnion(balcanica_temp, clip_poly_temp)
range_temp_balcanica <- raster::intersect(range_temp, balcanica_temp) 
range_temp_balcanica@data$subgroup <- "balcanica"

# rupicapra and cartusiana
clip_poly_temp <- raster::bind(tatrica_temp, carpatica_temp, asiatica_temp, caucasica_temp, balcanica_temp)
range_temp_rupicapra <- range_temp - clip_poly_temp
range_temp_rupicapra@data$subgroup <- "rupicapra cartusiana"

range_temp_try <- raster::bind(range_temp_tatrica, range_temp_carpatica, range_temp_asiatica, range_temp_caucasica, range_temp_balcanica, range_temp_rupicapra)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# removing tatrica because GMPD subgroup says there's one rupi and one tatrica there 
# and they are v isolated
range_temp_try <- raster::disaggregate(range_temp_try)
sp.gmpd.points <- sp.gmpd.points[terra::buffer(range_temp_try[range_temp_try$subgroup == "rupicapra cartusiana",], 0.5),]
sp.gmpd.points@data$subgroup <- "rupicapra cartusiana"
GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Rusa unicolor ####
# only one sample
species_temp <- "Rusa unicolor"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Saiga tatarica ####
# only have one sample point
species_temp <- "Saiga tatarica"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Spilogale gracilis ####
# I only have two amphialus samples for which the location data is inaccurate. 
# Recorded as Santa Rosa Island and Santa Cruz Island but coordinates do not map to correct locations
# amphialus is also not recognised in range polygons 

species_temp <- "Spilogale gracilis"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Spilogale putorius ####
# only have one sample location
species_temp <- "Spilogale putorius"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Suricata suricatta ####
# only one sample
species_temp <- "Suricata suricatta"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Sus scrofa ####
# 16 subspecies, with four regional groupings:
# Indonesian: vittatus (not represented)
# Eastern: sibiricus, ussuricus, leucomystax, riukiuanus, taivanus, moupinensis (none represented)
# Indian: davidi, cristatus (none represented)
# Western: scrofa, meridionalis, algira, attila, lybicus, nigripes
# only scrofa represented, but not possible to distinguish range boundaries in Iran and Iraq
# also lots of overlap in range descriptions of attila and lybicus
# can at least split into regional, algira, meridionalis, and majori
species_temp <- "Sus scrofa"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# Indonesian
clip_poly_temp <- countries50[countries50$name %in% c("Malaysia", "Indonesia"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.12)
range_temp_vittatus <- raster::intersect(range_temp, clip_poly_temp)
range_temp_vittatus@data$subgroup <- "vittatus"

# Eastern
clip_poly_temp <- matrix(c(95.3, 48.6,
                      95.3, 55.3,
                      147, 55.3,
                      147, 7.1,
                      105.1, 7.1,
                      100.4, 13.3,
                      99.9, 20.9,
                      96.4, 24.3,
                      97.6, 28.5,
                      97.2, 30.1,
                      95.3, 48.6),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - countries50[countries50$name == "Myanmar",]

range_temp_eastern <- raster::intersect(range_temp, clip_poly_temp)
range_temp_eastern@data$subgroup <- "eastern"

# Indian
clip_poly_temp <- matrix(c(64.6, 25.1,
                      69.8, 31.9,
                      73.9, 35.2,
                      98.3, 30.9,
                      101.3, 21.6,
                      103.4, 10.9, 
                      104.9, 4.6,
                      64.6, 4.6,
                      64.6, 25.1),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - range_temp_eastern - range_temp_vittatus

range_temp_indian <- raster::intersect(range_temp, clip_poly_temp)
range_temp_indian@data$subgroup <- "indian"

# algira
clip_poly_temp <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.1)
range_temp_algira <- raster::intersect(range_temp, clip_poly_temp)
range_temp_algira@data$subgroup <- "algira"

# meridionalis
clip_poly_temp <- matrix(c(-8.9, 36,
                      -1.5, 38.8,
                      0.5, 36.4,
                      -8.9, 36),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_meridionalis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_meridionalis <- raster::bind(range_temp_meridionalis, range_temp[range_temp$island %in% c("Sardinia", "Corsica"),])
range_temp_meridionalis@data$subgroup <- "meridionalis"

# majori
clip_poly_temp <- matrix(c(8.51, 44.3,
                      8.52, 44.9,
                      13.2, 44.9,
                      19.6, 39,
                      16.4, 37.1,
                      8.51, 44.3),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_majori <- raster::intersect(range_temp, clip_poly_temp)
range_temp_majori@data$subgroup <- "majori"

# Western
range_temp_try <- raster::bind(range_temp_vittatus, range_temp_eastern, range_temp_indian, range_temp_algira, range_temp_majori, range_temp_meridionalis)
range_temp_western <- range_temp - range_temp_try
range_temp_western@data$subgroup <- "western"

# finishing off
range_temp_try <- raster::bind(range_temp_try, range_temp_western)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "western"

rm(list = ls(pattern = "^range_temp"))

## Sylvicapra grimmia ####
# can't really do much with the available information, and range is continuously connected
"Fourteen subspecies were recognized by Grubb and Groves (2001) and Wilson (2013). Distribution is continuous, 
there are many cases of intergradation and geographical boundaries between forms have not been delineated 
accurately"
#map# species_temp <- "Sylvicapra grimmia"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Syncerus caffer ####
# already done <3
species_temp <- "Syncerus caffer"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# GMPD subgroup assignment
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp, 
                           buff = data.frame(subgroup = c("caffer", "nanus"), buff = c(1.5, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Taxidea taxus ####
# wiki:
"Four subspecies have been recognized based on differences in skull size and pelage colour (Long 1972): 
Taxidea taxus berlandieri, in the southern United States; 
T. t. jacksoni, in the north-central United States and southern Ontario in Canada; 
T. t. taxus, in the Great Plains ecosystem from the United States into the prairie provinces of Canada; 
and T. t. jeffersonii, in the western United States and southern British Columbia."
# Likely to be taxus, but range is continuous and boundaries unclear

#map# species_temp <- "Taxidea taxus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Tragelaphus angasii ####
# no subspecies described by IUCN, wiki, or msotw
#map# species_temp <- "Tragelaphus angasii"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Tragelaphus eurycerus ####
# only one sample location
species_temp <- "Tragelaphus eurycerus"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Tragelaphus oryx ####
# "Three subspecies of Common Eland have been recognized, although their validity requires investigation" - IUCN
# oryx, livingstonii, and pattersonianus
# according to wiki I have all three and their ranges are contiguous without clear barriers
#map# species_temp <- "Tragelaphus oryx"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Tragelaphus scriptus ####
# very complex taxonomic history (see iucn taxonomic notes and wiki discussion)
# As a result (and because range is largely continuous) it's not possible to divide into subspecies based on info I have
# although this would be desirable considering there are two divergent lineages that may be distinct species
#map# species_temp <- "Tragelaphus scriptus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Tragelaphus spekii ####
# unable to split based on basins because of lack of specificity in descriptions
"The species might even be monotypic,[6] however, based on different drainage systems, three distinct subspecies are currently recognised:[14][15]

T. s. spekii (Speke, 1863): Nile sitatunga or East African sitatunga. Found in the Nile watershed.
T. s. gratus (Sclater, 1880): Congo sitatunga or forest sitatunga. Found in western and central Africa.
T. s. selousi (W. Rothschild, 1898): Southern sitatunga or Zambezi sitatunga. Found in southern Africa."

#map# species_temp <- "Tragelaphus spekii"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Tragelaphus strepsiceros ####
# seems as though I have strepsiceros, which is connected to chora by rift valley
# https://doi.org/10.1046/j.1365-294x.2001.01205.x
# Based on above paper's finding that their genotypes are quite divergent, even in geographically close samples
# I'm going to split the polygon in Kenya, though the exact line will be a bit arbitrary without knowledge of a specific barrier
"T. s. strepsiceros – southern parts of the range from southern Kenya to Namibia, Botswana, and South Africa
T. s. chora – northeastern Africa from northern Kenya through Ethiopia to eastern Sudan, Somalia, and Eritrea
T. s. cottoni – Chad and western Sudan"

species_temp <- "Tragelaphus strepsiceros"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

clip_poly_temp <- matrix(c(26, 8,
                      17, 8,
                      17, 15,
                      26, 15,
                      26, 8),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_cottoni <- raster::intersect(range_temp, clip_poly_temp)
range_temp_cottoni@data$subgroup <- "cottoni"

clip_poly_temp <- matrix(c(37, -1.2,
                      35, -0.9,
                      32, -0.9,
                      32, 21,
                      48, 21,
                      48, -1.2,
                      37, -1.2),
                    ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string=CRS(proj4string(IUCN_Native_Data)))

range_temp_chora <- raster::intersect(range_temp, clip_poly_temp)
range_temp_chora@data$subgroup <- "chora"

range_temp_strepsiceros <- range_temp - clip_poly_temp - range_temp_cottoni
range_temp_strepsiceros@data$subgroup <- "strepsiceros"

range_temp_try <- raster::bind(range_temp_cottoni, range_temp_chora, range_temp_strepsiceros)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "strepsiceros"

rm(list = ls(pattern = "^range_temp"))

## Urocyon cinereoargenteus ####
# Many many subspecies, but only venezuelae (Colombia and Venezuela) is geographically isolated
species_temp <- "Urocyon cinereoargenteus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

range_temp_try <- raster::disaggregate(range_temp)
range_temp_try@data[2, "subgroup"] <- "venezuelae"
range_temp_try@data[c(1,3), "subgroup"] <- "northern"

range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# extending range_temp_try because of overlapping buffers
clip_points_temp <- matrix(c(-112.5, 45.5),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- SpatialPolygonsDataFrame(terra::buffer(clip_points_temp, 100000), 
                                        data = range_temp_try@data[which(range_temp_try@data$subgroup == "northern"),])

range_temp_try <- raster::bind(range_temp_try, clip_points_temp)

sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try,
                           buff = data.frame(subgroup = c("northern", "venezuelae"), buff = c(0.1, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))

## Urocyon littoralis ####
# Only have one sample location for each island so not useful to split by subspecies
# and also not biologically informative otherwise
"Six distinct subspecies are recognized, one on each of the islands where they occur:

San Miguel Island Fox (Urocyon littoralis littoralis (Baird, 1858)), San Miguel Island,
Santa Rosa Island Fox (U. l. santarosae Grinnell & Linsdale, 1930), Santa Rosa Island,
Santa Cruz Island Fox (U. l. santacruzae Merriam, 1903), Santa Cruz Island,
Santa Catalina Island Fox (U. l. catalinae Merriam, 1903), Santa Caralina Island,
San Nicolas Island Fox (U. l. dickeyi Grinnell & Linsdale, 1930), San Nicolas Island, and
San Clemente Island Fox (U. l. clementae Merriam, 1903), San Clemente Island."

# I have Santa Catalina (catalinae), San Nicolas (dickeyi), and San Clemente (clementae)
species_temp <- "Urocyon littoralis"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Ursus americanus ####
# subspecific taxonomy is too complex and mainland populations are too admixed to separate reasonably
# coupled with the fact that I have samples across the full range
# can therefore only really exclude island subspecies populations
# I don't have Haida Gwaii (carlottae), Dall Island (pugnax), Vancouver Island (vancouveri), Kenai (perniger), Newfoundland (hamiltoni)
# afaik these are the only geographically isolated subpops
# https://doi.org/10.1093/molbev/msv114
species_temp <- "Ursus americanus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# Dall island and Haida Gwaii (pugnax and carlottae)
range_temp_haida <- matrix(c(-134, 54.5,
                      -132.7, 55,
                      -132, 54.4,
                      -131, 54.2,
                      -130.5, 51.5,
                      -133.5, 53,
                      -134, 54.5),
                    ncol = 2, byrow = TRUE)
range_temp_haida <- Polygon(range_temp_haida)
range_temp_haida <- SpatialPolygons(list(Polygons(list(range_temp_haida), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- range_temp - range_temp_haida
range_temp_haida <- raster::intersect(range_temp, range_temp_haida)
range_temp_haida@data$subgroup <- "carlottae pugnax"

# Kenai peninsula (perniger)
range_temp_kenai <- matrix(c(-150.5, 61.1,
                      -149.36, 60.91,
                      -149.156, 60.91,
                      -149.01, 60.83,
                      -148.7, 60.78,
                      -148.57, 60.81,
                      -148.46, 60.81,
                      -148.23, 60.77,
                      -147.4, 60.4,
                      -148, 59.9,
                      -151.5, 58.9,
                      -152.5, 59.4,
                      -151.6, 60.7, 
                      -150.5, 61.1),
                    ncol = 2, byrow = TRUE)
range_temp_kenai <- Polygon(range_temp_kenai)
range_temp_kenai <- SpatialPolygons(list(Polygons(list(range_temp_kenai), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- range_temp_try - range_temp_kenai
range_temp_kenai <- raster::intersect(range_temp, range_temp_kenai)
range_temp_kenai@data$subgroup <- "perniger"

# Vancouver Island (vancouveri)
clip_points_temp <- matrix(c(-125.8, 49.8,
                        -126.7, 49.7,
                        -126.1, 49.3,
                        -125.8, 49.2,
                        -123.5, 48.8,
                        -123.4, 48.9,
                        -126.9, 50.65),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- raster::disaggregate(range_temp_try)
range_temp_vancouveri <- range_temp_try[clip_points_temp,]
range_temp_vancouveri@data$subgroup <- "vancouveri"
range_temp_try <- range_temp_try - range_temp_vancouveri

# Newfoundland (hamiltoni)
clip_points_temp <- matrix(c(-56.3, 48.6),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_hamiltoni <- range_temp_try[clip_points_temp,]
range_temp_hamiltoni@data$subgroup <- "hamiltoni"
range_temp_try <- range_temp_try - range_temp_hamiltoni

# mainland
range_temp_try <- raster::aggregate(range_temp_try, by = names(range_temp_try))
range_temp_try@data$subgroup <- "mainland"

range_temp_try <- raster::bind(range_temp_try, range_temp_haida, range_temp_kenai, range_temp_vancouveri, range_temp_hamiltoni)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)
# overwriting floridanus samples
GMPD_Spatial@data[GMPD_Spatial@data$HostCorrectedName == species_temp, "subgroup"] <- "mainland"

rm(list = ls(pattern = "^range_temp"))

## Ursus arctos ####
# Brown bear taxonomy and subspecies classification has been described as "formidable and confusing,"
# North America:
"DNA analysis shows that, apart from recent human-caused population fragmentation,[39] brown bears 
in North America are generally part of a single interconnected population system, with the exception 
of the population (or subspecies) in the Kodiak Archipelago, which has probably been isolated since 
the end of the last Ice Age.[40][41] "

species_temp <- "Ursus arctos"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# North America
# middendorffi
range_temp_middendorffi <- range_temp[range_temp$island == "Kodiak" & !is.na(range_temp$island),]
range_temp_middendorffi@data$subgroup <- "middendorffi"

# ungavaensis 
clip_poly_temp <- states50[states50$name %in% c("Québec", "Newfoundland and Labrador"),]
clip_poly_temp <- terra::buffer(clip_poly_temp)
range_temp_ungavaensis <- raster::intersect(range_temp, clip_poly_temp)
range_temp_ungavaensis@data$subgroup <- "ungavaensis"

# sitkensis
range_temp_sitkensis <- range_temp[range_temp$island %in% c("Admiralty", "Chichago & Baranof"),]
range_temp_sitkensis@data$subgroup <- "sitkensis"

# horribilis
clip_poly_temp <- countries50[countries50$continent == "North America", ]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.6)

poly_temp <- matrix(c(-131.1, 75.8,
                      -90.4, 75.8,
                      -90.4, 64.9,
                      -131.1, 64.9,
                      -131.1, 75.8),
                    ncol = 2, byrow = TRUE)
poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- gUnion(clip_poly_temp, poly_temp)

range_temp_try <- raster::intersect(range_temp, clip_poly_temp)
range_temp_horribilis <- range_temp_try - range_temp_sitkensis - range_temp_ungavaensis - range_temp_middendorffi
range_temp_horribilis@data$subgroup <- "horribilis"

range_temp_america <- raster::bind(range_temp_middendorffi, range_temp_ungavaensis, range_temp_sitkensis, range_temp_horribilis)
#map# tm_shape(range_temp_america) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

# Eurasia and North Africa
# crowtheri
clip_poly_temp <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.1)
range_temp_crowtheri <- raster::intersect(range_temp, clip_poly_temp)
range_temp_crowtheri@data$subgroup <- "crowtheri"

# pyrenaicus
clip_poly_temp <- countries50[countries50$name == "Spain",]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.8) - range_temp_crowtheri
range_temp_pyrenaicus <- raster::intersect(range_temp, clip_poly_temp)
range_temp_pyrenaicus@data$subgroup <- "pyrenaicus"

# marsicanus
clip_poly_temp <- matrix(c(13, 44,
                      16, 42,
                      13, 40, 
                      9, 43,
                      13, 44),
                     ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
range_temp_marsicanus <- raster::intersect(range_temp, clip_poly_temp)
range_temp_marsicanus@data$subgroup <- "marsicanus"

# syriacus
clip_poly_temp <-  matrix(c(23.9, 39, 
                       35, 43,
                       38, 45.5,
                       62, 43.6,
                       64, 37,
                       77, 27,
                       25, 27, 
                       23.9, 39),
                     ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_try <- raster::disaggregate(range_temp)
range_temp_syriacus <- range_temp_try[clip_poly_temp,]
range_temp_syriacus@data$subgroup <- "syriacus"

# pruinosus, isabellinus, gobiensis
clip_poly_temp <-  matrix(c(67, 48,
                       91, 45, 
                       106, 46,
                       106, 27,
                       77, 27,
                       75.5, 32,
                       74.8, 33,
                       74, 34.8,
                       72, 35.5,
                       70, 37,
                       65, 37),
                     ncol = 2, byrow = TRUE)
clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
##map# tm_shape(clip_poly_temp) + tm_polygons(alpha = 0) +tm_shape(range_temp) + tm_polygons(alpha = 0)

range_temp_pruinosus <- range_temp_try[clip_poly_temp,]
range_temp_pruinosus@data$subgroup <- "pruinosus isabellinus gobiensis"

# arctos
# reusing yenisei/angara line from alces alces, but going the opposite way round
YA_coords <- YA_coords[nrow(YA_coords):1,]
YA_curr_coords <- rbind(YA_coords,
                         matrix(c(106.5, 50.3,
                                  105.6, 48.1,
                                  90.9, 44.5, 
                                  77.1, 47.5, 
                                  33.7, 45,
                                  
                                  29.1, 41.19,
                                  29.07, 41.16,
                                  29.08, 41.13,
                                  29.07, 41.11,
                                  29.07, 41.09,
                                  
                                  20.8, 37.6,
                                  13.2, 44.8,
                                  6.1, 43,
                                  3.8, 43,
                                  3.8, 71.8, 
                                  81.2, 71.8,
                                  YA_coords[1,]), 
                                ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(YA_curr_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_arctos <- raster::intersect(range_temp, clip_poly_temp)
range_temp_arctos@data$subgroup <- "arctos"

# collaris beringianus lasiotus
YA_curr_coords <- rbind(YA_coords,
                            matrix(c(106.5, 50.3,
                                     105.6, 48.1,
                                     112, 39.6,
                                     110, 24,
                                     -99.2, 14.4,
                                     -172.3, 52.5,
                                     -169.1, 65.9,
                                     -169.1, 76,
                                     83.5, 76,
                                     YA_coords[1,]), 
                                   ncol = 2, byrow = TRUE))

clip_poly_temp <- Polygon(YA_curr_coords)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_collaris <- range_temp - clip_poly_temp
range_temp_collaris@data$subgroup <- "collaris beringianus lasiotus"

# adding together
range_temp_try <- raster::bind(range_temp_america, range_temp_crowtheri, range_temp_pyrenaicus, range_temp_syriacus, 
                         range_temp_pruinosus, range_temp_marsicanus, range_temp_arctos, range_temp_collaris)

#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
buff_temp <- data.frame(subgroup = c("horribilis", "middendorffi", "marsicanus", "arctos"), 
                                 buff = c(0, 0, 0, 0.5))
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = buff_temp)

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))

## Ursus maritimus ####
# no subpsecies on IUCN or msotw, and wiki states: 
"When the polar bear was originally documented, two subspecies were identified [...] This 
distinction has since been invalidated."
#map# species_temp <- "Ursus maritimus"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Vulpes corsac ####
# only one sample
species_temp <- "Vulpes corsac"
#map# unique(GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp,]@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,]
GMPD_Spatial <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,]

## Vulpes lagopus ####
"Besides the nominate subspecies, the common Arctic fox, V. l. lagopus, four other subspecies of this fox have been described:

Bering Islands Arctic fox, V. l. beringensis
Greenland Arctic fox, V. l. foragoapusis
Iceland Arctic fox, V. l. fuliginosus
Pribilof Islands Arctic fox, V. l. pribilofensis" # (wiki)

# I don't have beringensis or pribilofensis
# I do have foragoapusis, fuliginosus, and lagopus

species_temp <- "Vulpes lagopus"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# foragoapusis (Iceland)
range_temp_foragoapusis <- range_temp[!is.na(range_temp$island) & range_temp@data$island == "Iceland",]
range_temp_foragoapusis@data$subgroup <- "foragoapusis"

# fuliginosus (Greenland)
range_temp_fuliginosus <- range_temp[!is.na(range_temp$island) & range_temp$island == "Greenland",]
range_temp_try <- raster::disaggregate(range_temp)

clip_points_temp <- matrix(c(-53.375, 69.829),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))

range_temp_fuliginosus <- raster::bind(range_temp_fuliginosus, range_temp_try[clip_points_temp,])
range_temp_fuliginosus@data$subgroup <- "fuliginosus"

# beringensis
# polygon is missing some islands in Bering Sea where there are gbif samples
bering_temp <- c("St. Matthew", "Hall", "St. Lawrence")
range_temp_beringensis <- range_temp[range_temp$island %in% bering_temp, ]
range_temp_beringensis@data$subgroup <- "beringensis"

a_temp <-  matrix(c(165.757016, 55.283572,
                    165.965757, 55.361706,
                    166.273374, 55.308592,
                    166.255908, 55.183336,
                    166.591977, 54.909539,
                    166.674375, 54.678355,
                    166.070127, 55.063977,
                    166.081113, 55.123701,
                    165.989103, 55.220562),
                  ncol = 2, byrow = TRUE)
a_temp <- Polygon(a_temp)

b_temp <-  matrix(c(167.446164, 54.866095,
                    167.731809, 54.761639,
                    167.935056, 54.652145,
                    168.102597, 54.497719,
                    167.748288, 54.647378,
                    167.451657, 54.818649),
                  ncol = 2, byrow = TRUE)
b_temp <- Polygon(b_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(a_temp, b_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
data_temp <- range_temp@data[nrow(range_temp),]
data_temp$island <- "Bering"
data_temp$subgroup <- "beringensis"

clip_poly_temp <- SpatialPolygonsDataFrame(Sr          = clip_poly_temp,
                                      data        = data_temp, 
                                      match.ID    = FALSE)

range_temp_beringensis <- raster::bind(range_temp_beringensis, clip_poly_temp)

# pribilofensis
pribilof_temp <- c("St. Paul", "St. George")
range_temp_pribilofensis <- range_temp[range_temp$island %in% pribilof_temp, ]
range_temp_pribilofensis@data$subgroup <- "pribilofensis"

# lagopus
range_temp_lagopus <- range_temp[!range_temp$island %in% c(bering_temp, pribilof_temp, "Greenland", "Iceland"),]
range_temp_lagopus <- range_temp_lagopus - range_temp_try[clip_points_temp,]
range_temp_lagopus@data$subgroup <- "lagopus"

range_temp_try <- raster::bind(range_temp_lagopus, range_temp_foragoapusis, range_temp_fuliginosus, range_temp_beringensis, range_temp_pribilofensis)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
buff_temp <- data.frame(subgroup = c("fuliginosus", "lagopus", "foragoapusis"), 
                                 buff = c(5, 1, 0))
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = buff_temp)

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

rm(list = ls(pattern = "^range_temp"))

## Vulpes macrotis ####
"most available data suggest that kit foxes in the San Joaquin Valley of California are likely 
to warrant a subspecific designation, V. m. mutica, due to geographical isolation, and that any 
other kit foxes may be included in a second subspecies, V. m. macrotis"
species_temp <- "Vulpes macrotis"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>% 
  mutate(subgroup = case_when(binomial == species_temp & is.na(subgroup) ~ "macrotis",
                              TRUE ~ subgroup))

# GMPD subgroup assignment
# including point near boundary as it lies just within the San Joaquin Valley which is the geographic
# boundary described by IUCN as separating the two subspecies
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]

range_temp_mutica <- range_temp[range_temp$subgroup == "mutica", ]
range_temp_mutica <- gBuffer(range_temp_mutica, byid = TRUE, width = 0.6)

range_temp_macrotis <- range_temp[range_temp$subgroup == "macrotis", ]
range_temp_macrotis <- gBuffer(range_temp_macrotis, byid = TRUE, width = 1.2)
range_temp_macrotis <- range_temp_macrotis - range_temp_mutica

range_temp_try <- raster::bind(range_temp_mutica, range_temp_macrotis)

sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try, 
                           buff = data.frame(subgroup = c("macrotis", "mutica"), buff = c(0, 0)))

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

## Vulpes velox ####
# no subspecies in IUCN, wiki, or msotw
#map# species_temp <- "Vulpes velox"
#map# range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
#map# sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

## Vulpes vulpes ####
# 45 subspecies with patchy range descriptions which do not cover full extent of iucn range polygon
# In Europe I have: crucifera, silacea, vulpes
# In Africa: barbara
# In Eurasia, moving West to East: 
# palaestina/arabica, kurdistanica, flavescens, unknown, japonica, schrencki

species_temp <- "Vulpes vulpes"
range_temp <- IUCN_Native_Data[IUCN_Native_Data$binomial == species_temp, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == species_temp, ]
#map# tm_shape(range_temp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots("subgroup") 

# barbara and atlantica
clip_poly_temp <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia", "Libya"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.1)
range_temp_barbara <- raster::intersect(range_temp, clip_poly_temp)
range_temp_barbara@data$subgroup <- "barbara atlantica"

# silacea
clip_poly_temp <- countries50[countries50$name %in% c( "Spain", "Portugal", "Andorra"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.15)
range_temp_silacea <- raster::intersect(range_temp, clip_poly_temp)
range_temp_silacea@data$subgroup <- "silacea"

# ichnusae
range_temp_ichnusae <- range_temp[range_temp$island %in% c("Sardinia", "Corsica"),]
range_temp_ichnusae@data$subgroup <- "ichnusae"

# niloticus
range_temp_niloticus <- raster::disaggregate(range_temp)[countries50[countries50$name == "Sudan",],]
range_temp_niloticus@data$subgroup <- "niloticus"

# crucifera
clip_poly_temp <- countries50[countries50$continent == "Europe" & !countries50$name %in% c("Russia", "Finland", "Norway", "Sweden", "Spain", "Portugal", "Andorra"),]
clip_poly_temp <- terra::buffer(clip_poly_temp, 0.2)
poly_temp <- terra::buffer(countries50[countries50$name %in% c("Russia", "Sweden"),], 0.07)
poly_temp <- gUnion(poly_temp, terra::buffer(countries50[countries50$name == "Turkey",], 0.1))
clip_poly_temp <- clip_poly_temp - poly_temp - range_temp_ichnusae - range_temp_silacea

poly_temp <-  matrix(c(20.7, 56, 
                       25, 54.2, 
                       18.5, 54,
                       20.7, 56),
                     ncol = 2, byrow = TRUE)
poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- gUnion(clip_poly_temp, poly_temp)

range_temp_crucifera <- raster::intersect(range_temp, clip_poly_temp)
range_temp_crucifera@data$subgroup <- "crucifera"
range_temp_crucifera@polygons[[2]]@ID <- "2"

# arabica?
# vulpes (split at the Urals, and removing caucasica)
clip_poly_temp <- matrix(c(67.4, 68.8,
                      66.1, 68,
                      65.7, 67.2,
                      63.5, 66.5,
                      62.8, 65.8, 
                      61.1, 65,
                      59.9, 65.1,
                      58.9, 59.4,
                      59.7, 55.3, 
                      57.8, 54.5, 
                      57.1, 52.5,
                      57.2, 50.6,
                      46.8, 49.3,
                      40.1, 49.3, 
                      3.4, 55,
                      3.4, 71.3, 
                      64.1, 71.3,
                      67.4, 68.8),
                    ncol = 2, byrow = TRUE)

clip_poly_temp <- Polygon(clip_poly_temp)
clip_poly_temp <- SpatialPolygons(list(Polygons(list(clip_poly_temp), ID = "a")), proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_poly_temp <- clip_poly_temp - terra::buffer(countries50[countries50$name == "Kazakhstan",], 0.15) - range_temp_crucifera

range_temp_vulpes <- raster::intersect(range_temp, clip_poly_temp)
range_temp_vulpes@data$subgroup <- "vulpes"

# japonica
range_temp_japonica <- range_temp[range_temp$island %in% c("Kyushu", "Honshu"),]
range_temp_japonica@data$subgroup <- "japonica"

# schrencki
range_temp_schrencki <- range_temp[range_temp$island %in% "Hokkaido",]
range_temp_schrencki@data$subgroup <- "schrencki"

# splendidissima
range_temp_splendidissima <- range_temp[range_temp$island %in% "Kuril Islands",]
range_temp_splendidissima@data$subgroup <- "splendidissima"

# other subspecies
range_temp_try <- raster::bind(range_temp_barbara, range_temp_silacea, range_temp_ichnusae, range_temp_niloticus,
                         range_temp_crucifera, range_temp_vulpes, range_temp_japonica, range_temp_schrencki, range_temp_splendidissima)
range_temp_otherssp <- range_temp - range_temp_try
range_temp_otherssp@data$subgroup <- "unknown"

range_temp_try <- raster::bind(range_temp_try, range_temp_otherssp)
#map# tm_shape(range_temp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != species_temp,],
                                 range_temp_try)

# GMPD subgroup assignment
# extending range_temp_try because of overlapping buffers
clip_points_temp <- matrix(c(11.27, 55.27,
                        22.57, 58.58),
                      ncol = 2, byrow = TRUE)
clip_points_temp <- SpatialPoints(clip_points_temp, proj4string = CRS(proj4string(IUCN_Native_Data)))
clip_points_temp <- SpatialPolygonsDataFrame(terra::buffer(clip_points_temp, 100000), 
                                        data = range_temp_crucifera@data[which(range_temp_crucifera@data$SHAPE_Area > 6000),],
                                        match.ID = FALSE)

range_temp_try <- raster::bind(range_temp_try, clip_points_temp)

buff_temp <- data.frame(subgroup = c("barbara atlantica", "silacea", "crucifera", "vulpes", "unknown", "japonica", "schrencki"), 
                                 buff = c(0, 0, 0, 0, 0, 0, 0.1))
sp.gmpd.points <- pip_test(species_temp, dat = GMPD_Spatial, range.polygon = range_temp_try,
                           buff = buff_temp)

GMPD_Spatial <- raster::bind(GMPD_Spatial[GMPD_Spatial$HostCorrectedName != species_temp,], sp.gmpd.points)

# Finishing up ####

GMPD_Spatial@data[is.na(GMPD_Spatial@data$subgroup), "subgroup"] <- "not used"
IUCN_Native_Data@data[is.na(IUCN_Native_Data@data$subgroup), "subgroup"] <- "not used"

GMPD_Subgroups <- as.data.frame(GMPD_Spatial)

# tidying
rm(list = ls(pattern = "^range_temp"))
rm(list = ls(pattern = "_temp$"))
rm(list = ls(pattern = "^YA_"))
rm(turkey_rangepol, sp.gmpd.points, list = ls(pattern = "^clip_"))
rm(drop_ID, species_temp, list = ls(pattern = "^split_"))
