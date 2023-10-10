# Analyse distance metrics and climate density
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(MCMCglmm)
library(lme4)
library(here)
library(postMCMCglmm)
library(beepr)
# source(here::here("Functions.R"))

GMPD_Kernel_Data <- readRDS(here::here("Data/Data back ups/GMPD_Kernel_Data_04.rds"))

Parasite_KD <- GMPD_Kernel_Data$samples.out %>% 
  dplyr::bind_rows()

Parasite_KD <- Model_Data <- Parasite_KD %>%
  dplyr::filter(RestrAll,
                RangeMethod   == "iucn",
                RangeTaxonLvl == "species") %>%
  mutate(LatitudeScaled            = base::scale(abs(Latitude)),
         EquatorwardsPropScaled    = base::scale(EquatorwardsProp),
         
         MedianPropScaled          = base::scale(MedianProp),
         MedianPropSquared         = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled       = base::scale(MedianPropSquared),
         
         CorrectedDistanceScaled   = base::scale(CorrectedDistance),
         UncorrectedDistanceScaled = base::scale(UncorrectedDistance),
         
         Axis1Scaled               = base::scale(Axis1),
         Axis2Scaled               = base::scale(Axis2),
         
         ParasiteDetected          = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected        = HostsSampled - ParasiteDetected,
         HostSubgroup              = paste(HostCorrectedName, subgroup))

test_model_0 <- glm(formula = Prevalence ~ Axis1Scaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_0) 

test_model_1 <- glm(formula = Prevalence ~ Axis2Scaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_1) 

test_model_2 <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_2) 

test_model_3 <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + LatitudeScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_3) 

test_model_4 <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_4) 
# Keep axis1

test_model_5 <- glm(formula = Prevalence ~ Axis2Scaled + LatitudeScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_5) 
# Drop axis2

test_model_6 <- glm(formula = Prevalence ~ LatitudeScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_6) 
# 


test_model_7 <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_7) 

test_model_7a <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2),
                     data = Parasite_KD, family = binomial)
summary(test_model_7a) 


test_model_8 <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + distanceScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_8) 

test_model_8a <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + poly(distanceScaled, 2),
                     data = Parasite_KD, family = binomial)
summary(test_model_8a) 


test_model_9 <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + MedianPropScaled,
                    data = Parasite_KD, family = binomial)
summary(test_model_9) 

test_model_9a <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + poly(MedianPropSquared, 2),
                     data = Parasite_KD, family = binomial)
summary(test_model_9a) 


test_model_10 <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + MedianPropScaled + distanceScaled,
                     data = Parasite_KD, family = binomial)
summary(test_model_10) 

test_model_10a <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + poly(MedianPropSquared, 2) + distanceScaled,
                      data = Parasite_KD, family = binomial)
summary(test_model_10a) # best

test_model_10b <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + MedianPropScaled + poly(distanceScaled, 2),
                      data = Parasite_KD, family = binomial)
summary(test_model_10b) # worse

test_model_10c <- glm(formula = Prevalence ~ Axis1Scaled + LatitudeScaled + poly(MedianPropSquared, 2) + poly(distanceScaled, 2),
                      data = Parasite_KD, family = binomial)
summary(test_model_10c) # worse 


test_model_11 <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + distanceScaled,
                     data = Parasite_KD, family = binomial)
summary(test_model_11) # worse

test_model_11a <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + distanceScaled,
                      data = Parasite_KD, family = binomial)
summary(test_model_11a) # worse


test_model_12 <- glm(formula = Prevalence ~ Axis1Scaled + MedianPropScaled + distanceScaled,
                     data = Parasite_KD, family = binomial)
summary(test_model_12) # worse

test_model_12a <- glm(formula = Prevalence ~ Axis1Scaled + poly(MedianPropSquared, 2) + distanceScaled,
                      data = Parasite_KD, family = binomial)
summary(test_model_12a) # worse

test_model_13 <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + MedianPropScaled + distanceScaled,
                     data = Parasite_KD, family = binomial)
summary(test_model_13) # worse

test_model_13a <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + poly(MedianPropSquared, 2) + distanceScaled,
                      data = Parasite_KD, family = binomial)
summary(test_model_13a) # worse

test_model_14 <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled + Axis2Scaled + MedianPropScaled + distanceScaled,
                     data = Parasite_KD, family = binomial)
summary(test_model_14) # worse

test_model_14a <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled + Axis2Scaled + poly(MedianPropSquared, 2) + distanceScaled,
                      data = Parasite_KD, family = binomial)
summary(test_model_14a) # worse


saveRDS(ks.out, file = "KS_data")
saveRDS(ks.out2, file = "KS_data2")

ks.out <- ks.out2


