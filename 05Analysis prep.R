# Analysis and plotting of final dataset
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ###########################################################

library(here)               #
library(tidyverse)          #
library(corrplot)

source(here::here("Functions.R"))

GMPD_Climate_Data <- read.csv(here::here("Data/Data back ups/GMPD_Climate_Data_04.csv"), 
                              header = TRUE, 
                              stringsAsFactors = FALSE)

# temp until I rerun earlier scripts
GMPD_Parasite_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), 
                          header = TRUE, 
                          stringsAsFactors = FALSE) %>% 
  dplyr::select(ParasiteCorrectedName, ParType, ParPhylum, ParClass) %>% 
  dplyr::distinct() %>% 
  dplyr::right_join(GMPD_Climate_Data, by = join_by("ParasiteCorrectedName"))

# It was losing rows when both steps together
GMPD_Parasite_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), 
                               header = TRUE, 
                               stringsAsFactors = FALSE) %>% 
  dplyr::select(HostCorrectedName, Group) %>% 
  dplyr::distinct() %>% 
  dplyr::right_join(GMPD_Parasite_Data, by = join_by("HostCorrectedName")) %>% 
  dplyr::mutate(Group = case_when(!is.na(Group) ~ Group,
                                  str_detect(HostCorrectedName, "Tragelaphus oryx|Cervus canadensis|Equus quagga") ~ "ungulates",
                                  str_detect(HostCorrectedName, "Melogale subaurantiaca")                          ~ "carnivores"))


# Prep and data exploration ####################################################
## Hosts ####
# Check whether there are any host groups which should be removed 

pdf(here::here("Figures/Host Group.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = Group, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = Group, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = Group, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Host Order.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = HostOrder, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostOrder, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostOrder, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Host Family.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = HostFamily, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostFamily, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostFamily, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Host Species.pdf"), width = 10, height = 25)

basic_barplot(dat = GMPD_Climate_Data, xvar = HostCorrectedName, yvar = Prevalence) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
basic_barplot(dat = GMPD_Climate_Data, xvar = HostCorrectedName, yvar = Latitude) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
basic_barplot(dat = GMPD_Climate_Data, xvar = HostCorrectedName, yvar = AbsLatitude) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

dev.off()

## Group: uninformative considering we will include phylogeny
## Order: Doesn't give much more than Group. Arguably could remove Perissodactyla because of few rows, 
# but it's geographically restricted to regions where data is already lacking and will be accounted for in phylogeny
## Family: Too many factor levels to make much sense of it, and better to account for it in phylogeny

## Parasites ####
# Check whether there are any parasite groups which should be removed

pdf(here::here("Figures/Parasite Class.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = ParClass, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParClass, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParClass, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Parasite Type.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = ParType, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParType, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParType, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Parasite Phylum.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = ParPhylum, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParPhylum, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParPhylum, yvar = AbsLatitude)

dev.off()

## Class: too many to make sense of
## Type: Fungi and prions are a bit sketchy. Both have relatively few rows and are geographically restricted
# Fungi are all in Finland, UK, and Czechia/Slovakia. Prions are all in USA in one area
# removing these two groups from all analyses
## Phylum: Ascomycota and prions are removed above. 
# Sarcomastigophora is a bit more acceptable because it's within Protozoa, and
# because it has a broader geographic spread

GMPD_Analysis_Data <- GMPD_Parasite_Data %>% 
  dplyr::filter(!ParType %in% c("Prion", "Fungus") &
                  !is.na(ParType) &
                  HostsSampled != 0) # TODO: delete after rerunning?


# https://stats.stackexchange.com/questions/117497/comparing-between-random-effects-structures-in-a-linear-mixed-effects-model
# Can't compare AIC as I usually would for fixed effects so could use anova instead
# https://stats.stackexchange.com/questions/418593/understanding-anova-to-compare-mixed-model-with-a-gzlm

# Can't use anova because there is no difference in the number of parameters between models
# Could then return to comparing with AIC? The reason it's not used for random effects is 
# because it penalises according to the number of parameters and that's not the point for random
# But then Phylum will always come out on top just because it has more levels

# I'm going to go back to taking a more thoughtful approach

GMPD_IUCN_Species %>% 
  group_by(ParPhylum, ParType) %>% 
  summarise(n = n()) %>% 
  print() %>% 
  group_by(ParType) %>% 
  summarise(n = sum(n))

# The argument for Phylum is that it has more levels and so captures variation more accurately
# The flipside is true of Type which is coarser
# The argument against Phylum is that the levels are very uneven in their sampling
# That's where I think Type wins out, because the smaller samples coalesce without the
# groupings losing meaning as they are based on (hopefully) informative traits

  
## Correlations ####
GMPD_Analysis_Data %>% 
  dplyr::select(Prevalence, HostsSampled, Latitude, MedianProp, 
                CorrectedDensity, CorrectedDistance, CorrectedAngle,
                Axis1, Axis2, 
                wc2.1_10m_bio_5, wc2.1_10m_bio_6, wc2.1_10m_bio_12) %>% 
  cor() %>% 
  corrplot::corrplot()

# https://doi.org/10.1111/j.1600-0587.2012.07348.x
# None of the variables I plan on modelling 
# (Latitude, MedianProp, Axis1, Axis2, CorrectedDensity, CorrectedDistance, CorrectedAngle)
# have particularly high collinearity, though I will keep an eye on
# Latitude + Axis1, and CorrectedDistance + CorrectedDensity

# Analysis  ####################################################################
## Prepping IUCN Species data ####
GMPD_IUCN_Species <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "iucn",
                RangeTaxonLvl == "species") %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled        = base::scale(MedianProp),
                MedianPropSquared       = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled     = base::scale(MedianPropSquared),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                CorrectedDistanceLogged = base::scale(log(CorrectedDistance + 1)),
                CorrectedDistanceSqrt   = base::scale(sqrt(CorrectedDistance)),
                
                ParasiteDetected        = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected      = HostsSampled - ParasiteDetected,
                
                HostSubgroup            = paste(HostCorrectedName, subgroup, sep = "_"),
                HostParasite            = paste(HostCorrectedName, ParasiteCorrectedName, sep = "_"),
                
                RowID                   = row.names(.)) %>% 
  dplyr::arrange(desc(pmax(CorrectedDistance, CorrectedAngle))) %>% 
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistanceSquared = case_when(!TopHalf ~ CorrectedDistance * -1,
                                                     TRUE ~ CorrectedDistance),
                CorrectedDistanceSquScaled = base::scale(CorrectedDistanceSquared))

### Non-normal distributions ####
# TODO: put this in a sensible place
hist(GMPD_IUCN_Species$Latitude)
hist(abs(GMPD_IUCN_Species$Latitude))

hist(GMPD_IUCN_Species$MedianProp)
GMPD_IUCN_Species %>% 
  dplyr::mutate(MedianPropSquared = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                     TRUE       ~ MedianProp)) %>% 
  dplyr::pull(MedianPropSquared) %>% 
  hist()

hist(GMPD_IUCN_Species$Axis1)
hist(GMPD_IUCN_Species$Axis2)

hist(GMPD_IUCN_Species$CorrectedDensity)
hist(GMPD_IUCN_Species$CorrectedAngle)

hist(GMPD_IUCN_Species$CorrectedDistance)
hist(log(GMPD_IUCN_Species$CorrectedDistance+1))
hist(sqrt(GMPD_IUCN_Species$CorrectedDistance))
GMPD_IUCN_Species %>% 
  dplyr::arrange(desc(pmax(CorrectedDistance, CorrectedAngle))) %>% 
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistanceSquared = case_when(!TopHalf ~ CorrectedDistance * -1,
                                                     TRUE ~ CorrectedDistance)) %>% 
  dplyr::pull(CorrectedDistanceSquared) %>% 
  hist()
# I think the last is the most appropriate given that it's a true half-normal 
# distribution that was calculated from a 2d space
# Or at least, something approximating this approach
# Debatable whether I should use angle in sorting

## Prepping IUCN Subgroup data ####
GMPD_IUCN_Subgroup <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "iucn",
                RangeTaxonLvl == "subgroup") %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled        = base::scale(MedianProp),
                MedianPropSquared       = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled     = base::scale(MedianPropSquared),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                ParasiteDetected        = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected      = HostsSampled - ParasiteDetected,
                HostSubgroup            = paste(HostCorrectedName, subgroup))

## Prepping GBIF Species data ####
GMPD_GBIF_Species <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "gbif",
                RangeTaxonLvl == "species") %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled        = base::scale(MedianProp),
                MedianPropSquared       = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled     = base::scale(MedianPropSquared),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                ParasiteDetected        = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected      = HostsSampled - ParasiteDetected,
                HostSubgroup            = paste(HostCorrectedName, subgroup))

## Prepping GBIF Subgroup data ####
GMPD_GBIF_Subgroup <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "gbif",
                RangeTaxonLvl == "subgroup") %>%
  dplyr::mutate(LatitudeScaled            = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled          = base::scale(MedianProp),
                MedianPropSquared         = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled       = base::scale(MedianPropSquared),
                
                Axis1Scaled               = base::scale(Axis1),
                Axis2Scaled               = base::scale(Axis2),
                
                CorrectedDensityScaled    = base::scale(CorrectedDensity),
                CorrectedDistanceScaled   = base::scale(CorrectedDistance),
                CorrectedAngleScaled      = base::scale(CorrectedAngle),

                ParasiteDetected          = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected        = HostsSampled - ParasiteDetected,
                HostSubgroup              = paste(HostCorrectedName, subgroup))

## Frequentist ####
source(here::here("05.1Frequentist analysis.R"))

# TODO: go through asymmetric and reasoning
# TODO: add in model for : instead of * for where it doesn't converge
### Geographic ####
# Looking at fixed effects first, making sure to compare between unweighted, weighted, and random ID models

# We're relying mostly on random ID for final conclusions and generally disregarding 
# weighted models. When models are unweighted, the model is underdispersed (and 
# without weighting, incorrectly specified). With the model properly weighted, the 
# residuals are badly overdispersed and terms are spuriously accepted. To account 
# for overdispersion we add random ID intercepts for each row of the data. This 
# addresses the overdispersion well because it accounts for host traits, parasite 
# traits, and unmeasured environmental variables. It also permits us more flexibility 
# than modelling with a more interpretable random structure (e.g. random host slope) 
# because more complex random structures often fail to converge or have singular fit.

# TODO: example of model diagnostics?

#### Latitude ####
# starting with Latitude because that's what we're most interested in
# unweighted
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Lat, LM_IUCN_Species$Fixed$Null)
relative_likelihood(LM_IUCN_Species$Fixed$Lat, LM_IUCN_Species$Fixed$Null)

# weighted
summary(LM_IUCN_Species$Fixed_W$Lat)
AIC(LM_IUCN_Species$Fixed_W$Lat, LM_IUCN_Species$Fixed_W$Null)
relative_likelihood(LM_IUCN_Species$Fixed_W$Lat, LM_IUCN_Species$Fixed_W$Null)

# random ID
summary(LM_IUCN_Species$Fixed_ID$Lat)
AIC(LM_IUCN_Species$Fixed_ID$Lat, LM_IUCN_Species$FixedID$Null)
relative_likelihood(LM_IUCN_Species$Fixed_ID$Lat, LM_IUCN_Species$Fixed_ID$Null)

# Latitude is accepted

#### Latitude + MedianProp ####
# Next thing we're interested in is position within geographic range
# unweighted
summary(LM_IUCN_Species$Fixed$LatMed)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$Lat)

# weighted
summary(LM_IUCN_Species$Fixed_W$LatMed)
AIC(LM_IUCN_Species$Fixed_W$LatMed, LM_IUCN_Species$Fixed_W$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatMed, LM_IUCN_Species$Fixed_W$Lat)

# random ID
summary(LM_IUCN_Species$Fixed_ID$LatMed)
AIC(LM_IUCN_Species$Fixed_ID$LatMed, LM_IUCN_Species$Fixed_ID$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatMed, LM_IUCN_Species$Fixed_ID$Lat)

# MedianProp is accepted, but jumps from positive, negative, and back to positive between the models
# Because the random ID model is the most robust (and comparable to random host)
# That's the one we can rely on 
# Either way we want to include geographic range position

#### Latitude + Quadratic MedianProp ####
# Probably better to model it as a quadratic
summary(LM_IUCN_Species$Fixed$LatQMed)
AIC(LM_IUCN_Species$Fixed$LatQMed, LM_IUCN_Species$Fixed$LatMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMed, LM_IUCN_Species$Fixed$LatMed)

summary(LM_IUCN_Species$Fixed_W$LatQMed)
AIC(LM_IUCN_Species$Fixed_W$LatQMed, LM_IUCN_Species$Fixed_W$LatMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMed, LM_IUCN_Species$Fixed_W$LatMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMed)
AIC(LM_IUCN_Species$Fixed_ID$LatQMed, LM_IUCN_Species$Fixed_ID$LatMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMed, LM_IUCN_Species$Fixed_ID$LatMed)

# Definitely better as a quadratic
# Got the same thing of it switching signs, but same conclusion as previous

#### Latitude + Asymmetric Quadratic MedianProp ####
# Median prop might have an asymmetric pattern depending on which half of the range it's in
summary(LM_IUCN_Species$Fixed$LatQMedAsym)
AIC(LM_IUCN_Species$Fixed$LatQMedAsym, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsym, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsym)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsym, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsym, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsym) # fails to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsym, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsym, LM_IUCN_Species$Fixed_ID$LatQMed)

# It doesn't hold for the unweighted model but does for both weighted models

#### Latitude * Quadratic MedianProp ####
# Testing out the interaction term without the asymmetry first 
summary(LM_IUCN_Species$Fixed$LatQMedInt)
AIC(LM_IUCN_Species$Fixed$LatQMedInt, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedInt, LM_IUCN_Species$Fixed$LatQMed)
AIC(LM_IUCN_Species$Fixed$LatQMedInt, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedInt, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedInt)
AIC(LM_IUCN_Species$Fixed_W$LatQMedInt, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedInt, LM_IUCN_Species$Fixed_W$LatQMed)
AIC(LM_IUCN_Species$Fixed_W$LatQMedInt, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedInt, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedInt)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedInt, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedInt, LM_IUCN_Species$Fixed_ID$LatQMed)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedInt, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedInt, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# interaction isn't quite there in the unweighted model but better supported by the weighted ones
# It doesn't add as much as the asymmetry does to the model

#### Latitude * Asymmetric Quadratic MedianProp ####
# Testing out all terms at once
summary(LM_IUCN_Species$Fixed$LatQMedIntAsym)
AIC(LM_IUCN_Species$Fixed$LatQMedIntAsym, LM_IUCN_Species$Fixed$LatQMedInt)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedIntAsym, LM_IUCN_Species$Fixed$LatQMedInt)
# ^ comparing without asymmetry, v comparing without interaction
AIC(LM_IUCN_Species$Fixed$LatQMedIntAsym, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedIntAsym, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedIntAsym)
AIC(LM_IUCN_Species$Fixed_W$LatQMedIntAsym, LM_IUCN_Species$Fixed_W$LatQMedInt)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedIntAsym, LM_IUCN_Species$Fixed_W$LatQMedInt)
# ^ comparing without asymmetry, v comparing without interaction
AIC(LM_IUCN_Species$Fixed_W$LatQMedIntAsym, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedIntAsym, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedIntAsym)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedIntAsym, LM_IUCN_Species$Fixed_ID$LatQMedInt)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedIntAsym, LM_IUCN_Species$Fixed_ID$LatQMedInt)
# ^ comparing without asymmetry, v comparing without interaction
AIC(LM_IUCN_Species$Fixed_ID$LatQMedIntAsym, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedIntAsym, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# Asymmetry isn't quite there in the unweighted model, is included by the weighted, and random ID doesn't converge
# Same goes for the interaction.
# Will leave the interaction out for now as it added the least to the model

### Climatic ####
#### Latitude + Quadratic MedianProp + Axis1 ####
# First looking at pca axes as conceptual twin to Latitude
# without the asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAxis1)
AIC(LM_IUCN_Species$Fixed$LatQMedAxis1, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAxis1, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedAxis1)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAxis1, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAxis1, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAxis1)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAxis1, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAxis1, LM_IUCN_Species$Fixed_ID$LatQMed)

# with the asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAsymAxis1)
AIC(LM_IUCN_Species$Fixed$LatQMedAsymAxis1, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsymAxis1, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsymAxis1)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsymAxis1, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsymAxis1, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis1) # failed to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis1, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis1, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# Definitely not 

#### Latitude + Quadratic MedianProp + Axis2 ####
# without asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAxis2)
AIC(LM_IUCN_Species$Fixed$LatQMedAxis2, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAxis2, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedAxis2)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAxis2, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAxis2, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAxis2)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAxis2, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAxis2, LM_IUCN_Species$Fixed_ID$LatQMed)

# with asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAsymAxis2)
AIC(LM_IUCN_Species$Fixed$LatQMedAsymAxis2, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsymAxis2, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsymAxis2)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsymAxis2, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsymAxis2, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis2) # failed to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis2, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis2, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# Definitely not in unweighted or random ID, so not overall

#### Latitude + Quadratic MedianProp + Distance ####
# Distance from niche centre as a conceptual twin to MedianProp 
# without asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedDist)
AIC(LM_IUCN_Species$Fixed$LatQMedDist, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedDist, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedDist)
AIC(LM_IUCN_Species$Fixed_W$LatQMedDist, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedDist, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedDist)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedDist, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedDist, LM_IUCN_Species$Fixed_ID$LatQMed)

# with asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedDist)
AIC(LM_IUCN_Species$Fixed$LatQMedDist, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedDist, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedDist)
AIC(LM_IUCN_Species$Fixed_W$LatQMedDist, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedDist, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedDist) # Failed to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedDist, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedDist, LM_IUCN_Species$Fixed_ID$LatQMed)

# Nope, based on unweighted and random ID
# But that's quite nice really because I want Dist to be symmetric and therefore flat
# It should be modelled as a quadratic

#### Latitude + Quadratic MedianProp + Quadratic Distance ####
# Distance from niche centre as a conceptual twin to MedianProp 
# without asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedQDist)
AIC(LM_IUCN_Species$Fixed$LatQMedQDist, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedQDist, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedQDist)
AIC(LM_IUCN_Species$Fixed_W$LatQMedQDist, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedQDist, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedQDist)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedQDist, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedQDist, LM_IUCN_Species$Fixed_ID$LatQMed)

# with asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAsymQDist)
AIC(LM_IUCN_Species$Fixed$LatQMedAsymQDist, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsymQDist, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsymQDist)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsymQDist, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsymQDist, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# Only just improves the weighted and random ID models (strongly in weighted)
# Tentatively accept

#### Latitude + Quadratic MedianProp + Density ####
# Density  
# without asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedDens)
AIC(LM_IUCN_Species$Fixed$LatQMedDens, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedDens, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedDens)
AIC(LM_IUCN_Species$Fixed_W$LatQMedDens, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedDens, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedDens)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedDens, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedDens, LM_IUCN_Species$Fixed_ID$LatQMed)

# with asymmetry  
summary(LM_IUCN_Species$Fixed$LatQMedAsymDens)
AIC(LM_IUCN_Species$Fixed$LatQMedAsymDens, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsymDens, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsymDens)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsymDens, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsymDens, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymDens) # failed to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsymDens, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsymDens, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# Nope, based on unweighted and random ID

#### Latitude + Quadratic MedianProp + Quadratic Density ####
# Testing density as a quadratic
# without asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedQDens)
AIC(LM_IUCN_Species$Fixed$LatQMedQDens, LM_IUCN_Species$Fixed$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedQDens, LM_IUCN_Species$Fixed$LatQMed)

summary(LM_IUCN_Species$Fixed_W$LatQMedQDens)
AIC(LM_IUCN_Species$Fixed_W$LatQMedQDens, LM_IUCN_Species$Fixed_W$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedQDens, LM_IUCN_Species$Fixed_W$LatQMed)

summary(LM_IUCN_Species$Fixed_ID$LatQMedQDens)
AIC(LM_IUCN_Species$Fixed_ID$LatQMedQDens, LM_IUCN_Species$Fixed_ID$LatQMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedQDens, LM_IUCN_Species$Fixed_ID$LatQMed)

# with asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAsymQDens)
AIC(LM_IUCN_Species$Fixed$LatQMedAsymQDens, LM_IUCN_Species$Fixed$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsymQDens, LM_IUCN_Species$Fixed$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsymQDens)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsymQDens, LM_IUCN_Species$Fixed_W$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsymQDens, LM_IUCN_Species$Fixed_W$LatQMedAsym)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDens) 
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDens, LM_IUCN_Species$Fixed_ID$LatQMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDens, LM_IUCN_Species$Fixed_ID$LatQMedAsym)

# Nope again

#### Latitude * Quadratic Distance + Quadratic MedianProp ####
# Distance from niche centre as a conceptual twin to MedianProp 
# without asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedQIDist)
AIC(LM_IUCN_Species$Fixed$LatQMedQIDist, LM_IUCN_Species$Fixed$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedQIDist, LM_IUCN_Species$Fixed$LatQMedQDist)

summary(LM_IUCN_Species$Fixed_W$LatQMedQIDist)
AIC(LM_IUCN_Species$Fixed_W$LatQMedQIDist, LM_IUCN_Species$Fixed_W$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedQIDist, LM_IUCN_Species$Fixed_W$LatQMedQDist)

summary(LM_IUCN_Species$Fixed_ID$LatQMedQIDist) # failed to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedQIDist, LM_IUCN_Species$Fixed_ID$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedQIDist, LM_IUCN_Species$Fixed_ID$LatQMedQDist)

# with asymmetry
summary(LM_IUCN_Species$Fixed$LatQMedAsymQIDist)
AIC(LM_IUCN_Species$Fixed$LatQMedAsymQIDist, LM_IUCN_Species$Fixed$LatQMedAsymQDist)
relative_likelihood(LM_IUCN_Species$Fixed$LatQMedAsymQIDist, LM_IUCN_Species$Fixed$LatQMedAsymQDist)

summary(LM_IUCN_Species$Fixed_W$LatQMedAsymQIDist)
AIC(LM_IUCN_Species$Fixed_W$LatQMedAsymQIDist, LM_IUCN_Species$Fixed_W$LatQMedAsymQDist)
relative_likelihood(LM_IUCN_Species$Fixed_W$LatQMedAsymQIDist, LM_IUCN_Species$Fixed_W$LatQMedAsymQDist)

summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymQIDist) # failed to converge
AIC(LM_IUCN_Species$Fixed_ID$LatQMedAsymQIDist, LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist)
relative_likelihood(LM_IUCN_Species$Fixed_ID$LatQMedAsymQIDist, LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist)

# Only just improves the weighted and random ID models (strongly in weighted)
# Tentatively accept


### Comparisons ####
#### Axis1 against Latitude ####
summary(LM_IUCN_Species$Fixed$Axis1)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Axis1, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$Axis1, LM_IUCN_Species$Fixed$Lat)

summary(LM_IUCN_Species$Fixed_W$Axis1)
summary(LM_IUCN_Species$Fixed_W$Lat)
AIC(LM_IUCN_Species$Fixed_W$Axis1, LM_IUCN_Species$Fixed_W$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_W$Axis1, LM_IUCN_Species$Fixed_W$Lat)

summary(LM_IUCN_Species$Fixed_ID$Axis1)
summary(LM_IUCN_Species$Fixed_ID$Lat)
AIC(LM_IUCN_Species$Fixed_ID$Axis1, LM_IUCN_Species$Fixed_ID$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_ID$Axis1, LM_IUCN_Species$Fixed_ID$Lat)

# Latitude outperforms Axis1

#### Axis2 against Latitude ####
summary(LM_IUCN_Species$Fixed$Axis2)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Axis2, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$Axis2, LM_IUCN_Species$Fixed$Lat)

summary(LM_IUCN_Species$Fixed_W$Axis2)
summary(LM_IUCN_Species$Fixed_W$Lat)
AIC(LM_IUCN_Species$Fixed_W$Axis2, LM_IUCN_Species$Fixed_W$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_W$Axis2, LM_IUCN_Species$Fixed_W$Lat)

summary(LM_IUCN_Species$Fixed_ID$Axis2)
summary(LM_IUCN_Species$Fixed_ID$Lat)
AIC(LM_IUCN_Species$Fixed_ID$Axis2, LM_IUCN_Species$Fixed_ID$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_ID$Axis2, LM_IUCN_Species$Fixed_ID$Lat)

# Latitude outperforms Axis2

#### Axes against Latitude ####
summary(LM_IUCN_Species$Fixed$Axes)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Axes, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$Axes, LM_IUCN_Species$Fixed$Lat)

summary(LM_IUCN_Species$Fixed_W$Axes)
summary(LM_IUCN_Species$Fixed_W$Lat)
AIC(LM_IUCN_Species$Fixed_W$Axes, LM_IUCN_Species$Fixed_W$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_W$Axes, LM_IUCN_Species$Fixed_W$Lat)

summary(LM_IUCN_Species$Fixed_ID$Axes)
summary(LM_IUCN_Species$Fixed_ID$Lat)
AIC(LM_IUCN_Species$Fixed_ID$Axes, LM_IUCN_Species$Fixed_ID$Lat)
relative_likelihood(LM_IUCN_Species$Fixed_ID$Axes, LM_IUCN_Species$Fixed_ID$Lat)

# Latitude outperforms Axes

#### Distance against MedianProp ####
summary(LM_IUCN_Species$Fixed$Dist)
summary(LM_IUCN_Species$Fixed$Med)
AIC(LM_IUCN_Species$Fixed$Dist, LM_IUCN_Species$Fixed$Med)
relative_likelihood(LM_IUCN_Species$Fixed$Dist, LM_IUCN_Species$Fixed$Med)

summary(LM_IUCN_Species$Fixed_W$Dist)
summary(LM_IUCN_Species$Fixed_W$Med)
AIC(LM_IUCN_Species$Fixed_W$Dist, LM_IUCN_Species$Fixed_W$Med)
relative_likelihood(LM_IUCN_Species$Fixed_W$Dist, LM_IUCN_Species$Fixed_W$Med)

summary(LM_IUCN_Species$Fixed_ID$Dist)
summary(LM_IUCN_Species$Fixed_ID$Med)
AIC(LM_IUCN_Species$Fixed_ID$Dist, LM_IUCN_Species$Fixed_ID$Med)
relative_likelihood(LM_IUCN_Species$Fixed_ID$Dist, LM_IUCN_Species$Fixed_ID$Med)

# MedianProp is better

#### Quadratic Distance against Quadratic MedianProp ####
summary(LM_IUCN_Species$Fixed$QDist)
summary(LM_IUCN_Species$Fixed$QMed)
AIC(LM_IUCN_Species$Fixed$QDist, LM_IUCN_Species$Fixed$QMed)
relative_likelihood(LM_IUCN_Species$Fixed$QDist, LM_IUCN_Species$Fixed$QMed)

summary(LM_IUCN_Species$Fixed_W$QDist)
summary(LM_IUCN_Species$Fixed_W$QMed)
AIC(LM_IUCN_Species$Fixed_W$QDist, LM_IUCN_Species$Fixed_W$QMed)
relative_likelihood(LM_IUCN_Species$Fixed_W$QDist, LM_IUCN_Species$Fixed_W$QMed)

summary(LM_IUCN_Species$Fixed_ID$QDist)
summary(LM_IUCN_Species$Fixed_ID$QMed)
AIC(LM_IUCN_Species$Fixed_ID$QDist, LM_IUCN_Species$Fixed_ID$QMed)
relative_likelihood(LM_IUCN_Species$Fixed_ID$QDist, LM_IUCN_Species$Fixed_ID$QMed)

# MedianProp is still better

# Don't really have anything to compare Density with yet

# TODO: write a summary for frequentist analysis so far
# Look at differences in these effects when I add Host
# Same again for Parasite Type
# Same again for both Type and Host, but say where it ran out of beans

### Fixed effects Conclusions ####
# Final fixed effects model
summary(LM_IUCN_Species$Fixed$LatQMedAsymQDist)
summary(LM_IUCN_Species$Fixed_W$LatQMedAsymQDist)
summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist)

# diagnostics
sim_LatQMedAsymQDist_Fixed <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed$LatQMedAsymQDist)
sim_LatQMedAsymQDist_Fixed_W <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed_W$LatQMedAsymQDist)
sim_LatQMedAsymQDist_Fixed_ID <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist)

plot(sim_LatQMedAsymQDist_Fixed)
plot(sim_LatQMedAsymQDist_Fixed_W)
plot(sim_LatQMedAsymQDist_Fixed_ID)
# pretty rubbish without the random ID
# will have to see what the best random structure is

### Random Effects ####
# When looking at the random effects, although the maximal fixed effects model was:
summary(LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist)
# this repeatedly failed to converge across different random effect structures so we fall back on:
summary(LM_IUCN_Species$Fixed_ID$LatQMedQDist)

# The diagnostics are comparable to the more maximal model
sim_LatQMedQDist_Fixed <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed$LatQMedQDist)
sim_LatQMedQDist_Fixed_W <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed_W$LatQMedQDist)
sim_LatQMedQDist_Fixed_ID <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed_ID$LatQMedQDist)

plot(sim_LatQMedQDist_Fixed)
plot(sim_LatQMedQDist_Fixed_W)
plot(sim_LatQMedQDist_Fixed_ID)


#### Host focus ####
##### Host ####
summary(LM_IUCN_Species$Host_W$LatQMedQDist)
# Can't look at host slopes

# weighted
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Fixed_W$LatQMedQDist)
AIC(LM_IUCN_Species$Host_WIO$LatQMedQDist, LM_IUCN_Species$Fixed_W$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Host_WIO$LatQMedQDist, LM_IUCN_Species$Fixed_W$LatQMedQDist)
# Host is definitely better than nothing for both the weighted and unweighted models

sim_LatQMedQDist_Host_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Host_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Fixed_ID)
plot(sim_LatQMedQDist_Host_WIO)
# Pretty wiggly with host

##### Group/Host ####
summary(LM_IUCN_Species$HostGroup_W$LatQMedQDist)
# Can't look at slopes

# weighted
summary(LM_IUCN_Species$HostGroup_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$HostGroup_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$HostGroup_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
# Identical to host because there's no crossing

sim_LatQMedQDist_HostGroup_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$HostGroup_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_HostGroup_WIO)
# Identical to host

##### Group as fixed ####
# weighted
summary(LM_IUCN_Species$Fixed_Grp_ID$LatQMedQDist)
summary(LM_IUCN_Species$Fixed_ID$LatQMedQDist)

AIC(LM_IUCN_Species$Fixed_Grp_ID$LatQMedQDist, LM_IUCN_Species$Fixed_ID$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Fixed_Grp_ID$LatQMedQDist, LM_IUCN_Species$Fixed_ID$LatQMedQDist)
# Group makes it a little better

sim_LatQMedQDist_Fixed_Grp_ID <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed_Grp_ID$LatQMedQDist)

plot(sim_LatQMedQDist_Fixed_Grp_ID)

##### Group as fixed, Host random ####
summary(LM_IUCN_Species$Group_W$LatQMedQDist)
# Can't look at slopes

# weighted
summary(LM_IUCN_Species$Group_WIO$LatQMedQDist)
summary(LM_IUCN_Species$HostGroup_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Group_WIO$LatQMedQDist, LM_IUCN_Species$HostGroup_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Group_WIO$LatQMedQDist, LM_IUCN_Species$HostGroup_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$Group_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Group_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
# Group does marginally better as a fixed effect, but not enough that I would include it

sim_LatQMedQDist_Group_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Group_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_HostGroup_WIO)
plot(sim_LatQMedQDist_Group_WIO)
# Pretty wiggly still, but the shape has changed 

#### Parasite focus ####
##### Parasite species ####
# weighted
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Parasite_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Parasite_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
# Parasite does way better than host

sim_LatQMedQDist_Parasite_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Parasite_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_Parasite_WIO)
# Much better with parasite
# The slope goes away

##### Parasite Type ####
# weighted
summary(LM_IUCN_Species$Type_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Type_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Type_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Type_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Type_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
# Better than host but not as good as parasite species

sim_LatQMedQDist_Type_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Type_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_Type_WIO)
# Gets rid of the slope!

##### ParType/Parasite ####
# weighted
summary(LM_IUCN_Species$ParType_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Type_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$ParType_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$ParType_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$ParType_WIO$LatQMedQDist, LM_IUCN_Species$Type_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$ParType_WIO$LatQMedQDist, LM_IUCN_Species$Type_WIO$LatQMedQDist)
# Just about better than parasite alone, and definitely better than type alone

sim_LatQMedQDist_ParType_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$ParType_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_ParType_WIO)
plot(sim_LatQMedQDist_Parasite_WIO)
plot(sim_LatQMedQDist_Type_WIO)
# Looks nice

##### ParType fixed ####
# weighted
summary(LM_IUCN_Species$Fixed_Typ_ID$LatQMedQDist)
# failure to converge
summary(LM_IUCN_Species$Fixed_ID$LatQMedQDist)

AIC(LM_IUCN_Species$Fixed_Typ_ID$LatQMedQDist, LM_IUCN_Species$Fixed_ID$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Fixed_Typ_ID$LatQMedQDist, LM_IUCN_Species$Fixed_ID$LatQMedQDist)
# Similar story to group, it's good as a fixed effect

sim_LatQMedQDist_Fixed_Typ_ID <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed_Typ_ID$LatQMedQDist)

plot(sim_LatQMedQDist_Fixed_Typ_ID)

##### ParType fixed, Parasite random ####
# weighted
summary(LM_IUCN_Species$TypePar_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Type_WIO$LatQMedQDist)
summary(LM_IUCN_Species$ParType_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$TypePar_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$TypePar_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$TypePar_WIO$LatQMedQDist, LM_IUCN_Species$Type_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$TypePar_WIO$LatQMedQDist, LM_IUCN_Species$Type_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$TypePar_WIO$LatQMedQDist, LM_IUCN_Species$ParType_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$TypePar_WIO$LatQMedQDist, LM_IUCN_Species$ParType_WIO$LatQMedQDist)
# Similar story to group, it's good as a fixed effect but not much better than as a random

sim_LatQMedQDist_TypePar_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$TypePar_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_TypePar_WIO)


#### Host and Parasite ####
##### Host + ParType ####
# weighted
summary(LM_IUCN_Species$Both_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$Both_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Both_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Both_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Both_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Both_WIO$LatQMedQDist, LM_IUCN_Species$Type_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Both_WIO$LatQMedQDist, LM_IUCN_Species$Type_WIO$LatQMedQDist)
# Better than host, better than type, but still not as good as parasite species

sim_LatQMedQDist_Both_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Both_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_Type_WIO)
plot(sim_LatQMedQDist_Both_WIO)
# Doesn't get rid of the slope this time!

##### Host:Parasite Type ####
summary(LM_IUCN_Species$Nested_W$LatQMedQDist)
# Can't look at slopes

# weighted
summary(LM_IUCN_Species$Nested_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$Combo_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Nested_WIO$LatQMedQDist, LM_IUCN_Species$Combo_WIO$LatQMedQDist)
# Better than host, better than type, better than parasite species

sim_LatQMedQDist_Nested_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Nested_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_Nested_WIO)
# Looks great!

##### Host + Parasite ####
# weighted
summary(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist, LM_IUCN_Species$ParType_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist, LM_IUCN_Species$ParType_WIO$LatQMedQDist)
# Better than host, better than type, better than parasite species
# But not better than both species combined

sim_LatQMedQDist_BothSpecies_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_BothSpecies_WIO)
# Looks great!

##### Host Parasite combo ####
summary(LM_IUCN_Species$Combo_W$LatQMedQDist)
# Can't look at host slopes

# weighted
summary(LM_IUCN_Species$Combo_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
summary(LM_IUCN_Species$Host_WIO$LatQMedQDist)

AIC(LM_IUCN_Species$Combo_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Combo_WIO$LatQMedQDist, LM_IUCN_Species$Host_WIO$LatQMedQDist)
AIC(LM_IUCN_Species$Combo_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
relative_likelihood(LM_IUCN_Species$Combo_WIO$LatQMedQDist, LM_IUCN_Species$Parasite_WIO$LatQMedQDist)
# Better than host or parasite!

sim_LatQMedQDist_Combo_WIO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Combo_WIO$LatQMedQDist)

plot(sim_LatQMedQDist_Host_WIO)
plot(sim_LatQMedQDist_Combo_WIO)
# Gets rid of the slope!
# TODO: rewrite notes a little

#### Host Phylo ####


## Bayesian ####
# source(here::here("05.3Bayesian analysis.R"))
### Conclusions ####

# Plots ########################################################################

