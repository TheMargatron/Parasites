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

# TODO: tidy up
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


