# Phylogeny prep
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# libraries

# devtools::install_github("jinyizju/U.PhyloMaker")
library(U.PhyloMaker)
library(tidyverse)
library(DHARMa)
library(phyr)

megatree <- read.tree('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_megatree.tre')
sp.list  <- read.csv('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_sample_species_list.csv', sep=",")
gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_genus_list.csv', sep=",")

names(sp.list)

# Full phylo ####
sp.list <- GMPD_Analysis_Data %>% 
  mutate(species = case_when(HostCorrectedName == "Melogale subaurantiaca" ~ "Melogale moschata",
                             HostCorrectedName == "Cervus canadensis"      ~ "Cervus elaphus",
                             HostCorrectedName == "Lycalopex gymnocercus"  ~ "Pseudalopex gymnocercus",
                             TRUE ~ HostCorrectedName),
         genus = str_split_i(species, " ", 1)) %>% 
  select(species, genus) %>% 
  distinct()
# Down one because we merge Cervus canadensis and Cervus elaphus in this data source

result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result

# Create a distance matrix for phylo in residuals
phyloMat <- cophenetic.phylo(result) %>% as.data.frame.array() 
phyloMat <- phyloMat %>% 
  bind_rows(phyloMat[rownames(phyloMat) == "Cervus_elaphus",]) %>% 
  bind_cols(.$Cervus_elaphus) %>% 
  rename(Cervus_canadensis = "...109",
         # Cervus_elaphus = "Cervus_elaphus...66",
         Melogale_subaurantiaca = Melogale_moschata,
         Lycalopex_gymnocercus = Pseudalopex_gymnocercus) 

rownames(phyloMat)[rownames(phyloMat) == "Cervus_elaphus...109"] <- "Cervus_canadensis"
rownames(phyloMat)[rownames(phyloMat) == "Cervus_elaphus...66"] <- "Cervus_elaphus"
rownames(phyloMat)[rownames(phyloMat) == "Melogale_moschata"] <- "Melogale_subaurantiaca"
rownames(phyloMat)[rownames(phyloMat) == "Pseudalopex_gymnocercus"] <- "Lycalopex_gymnocercus"

Matrix_Data <- GMPD_IUCN_Species %>% 
  dplyr::select(HostCorrectedName) %>% 
  mutate(species = str_replace_all(HostCorrectedName, " ", "_")) %>% 
  dplyr::select(species) %>% 
  left_join(rownames_to_column(phyloMat), by = join_by(species == rowname)) 
  
distMat <- as.matrix(Matrix_Data[,Matrix_Data$species],
                     dimnames = Matrix_Data$species)

# model phylo ####
sp.list <- GMPD_IUCN_Species %>% 
  mutate(species = case_when(HostCorrectedName == "Melogale subaurantiaca" ~ "Melogale moschata",
                             HostCorrectedName == "Cervus canadensis"      ~ "Cervus elaphus",
                             HostCorrectedName == "Lycalopex gymnocercus"  ~ "Pseudalopex gymnocercus",
                             TRUE ~ HostCorrectedName),
         genus = str_split_i(species, " ", 1)) %>% 
  select(species, genus) %>% 
  distinct()
# Down one because we merge Cervus canadensis and Cervus elaphus in this data source

result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result

# Create a distance matrix for phylo in residuals
phyloMat <- cophenetic.phylo(result) %>% as.data.frame.array() 
phyloMat <- phyloMat %>% 
  bind_rows(phyloMat[rownames(phyloMat) == "Cervus_elaphus",]) %>% 
  bind_cols(.$Cervus_elaphus) %>% 
  rename(Cervus_canadensis = "...94",
         Melogale_subaurantiaca = Melogale_moschata,
         Lycalopex_gymnocercus = Pseudalopex_gymnocercus) 

rownames(phyloMat)[rownames(phyloMat) == "Cervus_elaphus...94"] <- "Cervus_canadensis"
rownames(phyloMat)[rownames(phyloMat) == "Cervus_elaphus...60"] <- "Cervus_elaphus"
rownames(phyloMat)[rownames(phyloMat) == "Melogale_moschata"] <- "Melogale_subaurantiaca"
rownames(phyloMat)[rownames(phyloMat) == "Pseudalopex_gymnocercus"] <- "Lycalopex_gymnocercus"

Matrix_Data <- GMPD_IUCN_Species %>% 
  dplyr::select(HostCorrectedName) %>% 
  mutate(species = str_replace_all(HostCorrectedName, " ", "_")) %>% 
  dplyr::select(species) %>% 
  left_join(rownames_to_column(phyloMat), by = join_by(species == rowname)) 

distMat <- as.matrix(Matrix_Data[,Matrix_Data$species],
                     dimnames = Matrix_Data$species)

sim_LatQMedQDist_Nested_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Nested_WIO$LatQMedQDist)
sim_LatQMedQDist_Nested_WIO2 <- recalculateResiduals(sim_LatQMedQDist_Nested_WIO, group = GMPD_IUCN_Species$HostCorrectedName)
testSpatialAutocorrelation(simulationOutput = sim_LatQMedQDist_Nested_WIO2, distMat = phyloMat)

# TODO: rerun model with altered phylo

# restricted model phylo ####
phylo_model_data <- GMPD_IUCN_Species %>% filter(HostCorrectedName != "Cervus canadensis")
phylo_model  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName:ParasiteCorrectedName), data = phylo_model_data, family = binomial, weights = HostsSampled), silent = TRUE)
phylo_model_basic  <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = phylo_model_data, family = binomial, weights = HostsSampled)

phylo_model_data <- phylo_model_data %>% 
  mutate(species = case_when(HostCorrectedName == "Melogale subaurantiaca" ~ "Melogale moschata",
                             HostCorrectedName == "Lycalopex gymnocercus"  ~ "Pseudalopex gymnocercus",
                             TRUE ~ HostCorrectedName)) 
sp.list <- phylo_model_data %>% 
  mutate(genus = str_split_i(species, " ", 1)) %>% 
  select(species, genus) %>% 
  distinct()

result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result

# Create a distance matrix for phylo in residuals
phyloMat <- cophenetic.phylo(result) %>% as.data.frame.array() 

rownames(phyloMat)[rownames(phyloMat) == "Melogale_moschata"] <- "Melogale_subaurantiaca"
rownames(phyloMat)[rownames(phyloMat) == "Pseudalopex_gymnocercus"] <- "Lycalopex_gymnocercus"

# Matrix_Data <- phylo_model_data %>% 
#   dplyr::select(HostCorrectedName) %>% 
#   mutate(species = str_replace_all(HostCorrectedName, " ", "_")) %>% 
#   dplyr::select(species) %>% 
#   left_join(rownames_to_column(phyloMat), by = join_by(species == rowname)) 
# 
# distMat <- as.matrix(Matrix_Data[, Matrix_Data$species],
#                      dimnames = Matrix_Data$species)

phylo_model_sim <- DHARMa::simulateResiduals(fittedModel = phylo_model)
phylo_model_sim2 <- recalculateResiduals(phylo_model_sim, group = phylo_model_data$HostCorrectedName)
testSpatialAutocorrelation(simulationOutput = phylo_model_sim2, distMat = phyloMat)
# looks alright?
phylo_model_basic_sim <- DHARMa::simulateResiduals(fittedModel = phylo_model_basic)
phylo_model_basic_sim2 <- recalculateResiduals(phylo_model_basic_sim, group = phylo_model_data$HostCorrectedName)
testSpatialAutocorrelation(simulationOutput = phylo_model_basic_sim2, distMat = phyloMat)

# Trying a phylogenetic model
phylo_model_data <- GMPD_IUCN_Species %>% filter(HostCorrectedName != "Cervus canadensis")
phylo_model_nested <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName:ParasiteCorrectedName), data = phylo_model_data, family = binomial, weights = HostsSampled), silent = TRUE)
phylo_model_basic  <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = phylo_model_data, family = binomial, weights = HostsSampled)
phylo_model_host   <- try(phyr::pglmm(formula = cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2) + 
                                        poly(CorrectedDistanceSquared, 2) + (1|species),
                                      data = phylo_model_data, family = "binomial",
                                      cov_ranef = list(species = result)), silent = FALSE)


