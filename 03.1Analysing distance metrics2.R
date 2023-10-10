# Analysing distance metrics 
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(MCMCglmm)
library(lme4)
library(here)
library(postMCMCglmm)
library(beepr)
source(here::here("Functions.R"))

GMPD_Distances_Data <- read.csv(here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), 
                                header = TRUE, 
                                stringsAsFactors = FALSE)

# Model naming convention:
# LM/MC     = Linear model (LM) or MCMCglmm (MC)
# IUCN/GBIF = method used for restricting data and calculating position within range
# SP/SG     = position within range calculated from position within species (SP) or subgroup (SG) range
# F/R/G     = Fixed effects only (F) or random as well (R)

# Number order:
# 1 = Latitude 
# 2 = Latitude + MedianProp
# 3 = Latitude * MedianProp
# 4 = Latitude + poly(MedianPropSquared)
# 5 = Latitude + poly(MedianPropSquared):AboveMedn

# Letter order for random effects:
# no letter = estimated full vcv matrix for host
# a = No covariance estimated, only variance in host
# b = Misc alterations to host
# c = estimated full vcv matrix for host:subgroup
# d = No covariance estimated, only variance in host:subgroup
# e = Misc alterations to host:subgroup

# IUCN species #################################################################

## Prepping data ###############################################################
GMPD_IUCN_Species <- GMPD_Distances_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "iucn",
                RangeTaxonLvl == "species") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared),
         ParasiteDetected       = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected     = HostsSampled - ParasiteDetected,
         HostSubgroup           = paste(HostCorrectedName, subgroup))

## GLMER #######################################################################
LM_IUCN_SP <- list()

### Fixed effects only #########################################################

# 1 ____________________________________________________________________________
LM_IUCN_SP_F_1 <- glm(formula = Prevalence ~ LatitudeScaled,
                      data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_F_1) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_F_1)

# Accept Latitude

# 2 ____________________________________________________________________________
LM_IUCN_SP_F_2 <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled,
                      data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_F_2) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_F_2)

# Accept MedianProp

# 3 ____________________________________________________________________________
LM_IUCN_SP_F_3 <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropScaled,
                      data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_F_3) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_F_3)

# AIC is reduced but interaction term is non-significant
# Reject interaction between Latitude and MedianProp

# 4 ____________________________________________________________________________
LM_IUCN_SP_F_4 <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2),
                      data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_F_4) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_F_4)

# Accept quadratic based on AIC 

# 5 ____________________________________________________________________________
LM_IUCN_SP_F_5 <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                      data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_F_5) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_F_5)

# Accept AboveMedn

# 6
LM_IUCN_SP_F_6 <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropSquared:AboveMedn,
                      data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_F_6) 


### Random effects #############################################################

# 1 ____________________________________________________________________________
# ----- Host
LM_IUCN_SP_R_1 <- glmer(formula = Prevalence ~ LatitudeScaled +
                          (1 + LatitudeScaled|HostCorrectedName),
                        data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_1) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_1)

# ----- Subgroup
LM_IUCN_SP_R_1c <- glmer(formula = Prevalence ~ LatitudeScaled +
                           (1 + LatitudeScaled|HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_1c) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_1c)

# Accept both host and subgroup

# 2 ____________________________________________________________________________
# ----- Host
LM_IUCN_SP_R_2 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                          (1 + LatitudeScaled + MedianPropScaled|HostCorrectedName),
                        data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_2)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_2)

# singular fit

LM_IUCN_SP_R_2a <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                           (1 + LatitudeScaled + MedianPropScaled||HostCorrectedName),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_2a) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_2a)

# Ideally would model full vcv but this will do for now

# ----- Subgroup
LM_IUCN_SP_R_2c <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                           (1 + LatitudeScaled + MedianPropScaled|HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_2c)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_2c)

# failed to converge

LM_IUCN_SP_R_2d <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                           (1 + LatitudeScaled + MedianPropScaled||HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_2d) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_2d)

# AIC is reduced by adding both host and subgroup
# Accept host and subgroup as random effects

# 3 ____________________________________________________________________________
# Skipping 3 because interaction term was not included in fixed effects

# 4 ____________________________________________________________________________
# ----- Host
LM_IUCN_SP_R_4 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                          (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                        data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4)

# failed to converge

LM_IUCN_SP_R_4a <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2)||HostCorrectedName),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4a)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4a)

# failed to converge

LM_IUCN_SP_R_4b <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 1)|HostCorrectedName),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4b)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4b)

# singular fit

LM_IUCN_SP_R_4bb <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1)||HostCorrectedName),
                          data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4bb) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4bb)

# 4 would be nice, but 4bb is the best we can do
# based on AIC it's better to ignore quadratic term for now

# ----- Subgroup
LM_IUCN_SP_R_4c <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4c)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4c)

# failed to converge

LM_IUCN_SP_R_4d <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2)||HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4d)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4d)

# failed to converge

LM_IUCN_SP_R_4e <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 1)|HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4e)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4e)

# failed to converge

LM_IUCN_SP_R_4ee <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1)||HostCorrectedName:HostSubgroup),
                          data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_4ee) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_4ee)

# again, based on AIC it's better to ignore quadratic term for now
# But accept host and subgroup

# 5 ____________________________________________________________________________
# ----- Host
LM_IUCN_SP_R_5 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                          (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn|HostCorrectedName),
                        data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5)

# failed to converge

LM_IUCN_SP_R_5a <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn||HostCorrectedName),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5a)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5a)

# failed to converge

LM_IUCN_SP_R_5b <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn|HostCorrectedName),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5b)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5b)

# singular fit

LM_IUCN_SP_R_5bb <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn||HostCorrectedName),
                          data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5bb) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5bb)

# 5 would be nice, but 5bb is the best we can do
# based on AIC it's better to ignore quadratic term for now

# ----- Subgroup
LM_IUCN_SP_R_5c <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn|HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5c)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5c)

# failed to converge

LM_IUCN_SP_R_5d <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn||HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5d)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5d)

# failed to converge

LM_IUCN_SP_R_5e <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn|HostCorrectedName:HostSubgroup),
                         data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5e)
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5e)

# failed to converge

LM_IUCN_SP_R_5ee <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn||HostCorrectedName:HostSubgroup),
                          data = GMPD_IUCN_Species, family = binomial)
summary(LM_IUCN_SP_R_5ee) 
LM_IUCN_SP <- append(LM_IUCN_SP, LM_IUCN_SP_R_5ee)

# again, based on AIC it's better to ignore quadratic term for now
# But accept host and subgroup

# from linear models we would accept Latitude and MedianProp, with host and subgroup as nested random effects
# But have to leave out the quadratic term in order to include random effects

write_rds(LM_IUCN_SP, here::here("Data/Model back ups/LM_IUCN_SP.rds"))

## MCMC ########################################################################
MC_IUCN_SP <- list()

### Fixed effects only #########################################################

# 1 ____________________________________________________________________________
MC_IUCN_SP_F_1 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_F_1)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_F_1)

# Accept latitude

# 2 ____________________________________________________________________________
MC_IUCN_SP_F_2 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_F_2)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_F_2)

# Accept MedianProp

# 3 ____________________________________________________________________________
MC_IUCN_SP_F_3 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled * MedianPropScaled,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_F_3)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_F_3)

# Could accept interaction between Latitude and MedianProp based on "significance"
# But the DIC is lower without the interaction term which agrees with LM approach
# Reject interaction between Latitude and MedianProp

# 4 ____________________________________________________________________________
MC_IUCN_SP_F_4 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_F_4)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_F_4)

# Accept quadratic based on AIC

# 5 ____________________________________________________________________________
MC_IUCN_SP_F_5 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_F_5) 
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_F_5)

# Accept AboveMedn

# Fixed effect results came out roughly the same for MC as for LM

### Random effects #############################################################

# 1 ____________________________________________________________________________
# ----- Host
MC_IUCN_SP_R_1 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                               random = ~ us(1 + LatitudeScaled):HostCorrectedName,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_R_1)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_1)

# Accept host

# ----- Subgroup
MC_IUCN_SP_R_1c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                                random = ~ us(1 + LatitudeScaled):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_1c)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_1c)

# DIC increases when I add subgroup in
# Accept host as a random effect but not subgroup?

# 2 ____________________________________________________________________________
# ----- Host
MC_IUCN_SP_R_2 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                           random = ~ us(1 + LatitudeScaled + MedianPropScaled):HostCorrectedName,
                           data    = GMPD_IUCN_Species,
                           family  = "multinomial2",
                           nitt    = 130000,
                           thin    = 100,
                           burnin  = 30000),
                   silent = TRUE)
summary(MC_IUCN_SP_R_2)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_2)

MC_IUCN_SP_R_2a <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                           random = ~ idh(1 + LatitudeScaled + MedianPropScaled):HostCorrectedName,
                           data    = GMPD_IUCN_Species,
                           family  = "multinomial2",
                           nitt    = 130000,
                           thin    = 100,
                           burnin  = 30000),
                    silent = TRUE)
summary(MC_IUCN_SP_R_2a)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_2a)

# ----- Subgroup
MC_IUCN_SP_R_2c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                           random = ~ us(1 + LatitudeScaled + MedianPropScaled):(HostCorrectedName + HostSubgroup),
                           data    = GMPD_IUCN_Species,
                           family  = "multinomial2",
                           nitt    = 130000,
                           thin    = 100,
                           burnin  = 30000),
                    silent = TRUE)
summary(MC_IUCN_SP_R_2c)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_2c)

# Often get "ill-conditioned G/R structure"
# DIC of host vs host and subgroup were identical when it worked so hard to say which is best random structure

MC_IUCN_SP_R_2d <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                random = ~ idh(1 + LatitudeScaled + MedianPropScaled):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_2d)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_2d)

# but definitely include host

# 4 ____________________________________________________________________________
# ----- Host
MC_IUCN_SP_R_4 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                               random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2)):HostCorrectedName,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_R_4)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_4)

# ill-conditioned G/R structure

MC_IUCN_SP_R_4a <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2)):HostCorrectedName,
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_4a)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_4a)

# Effective sample size is small for random effects so may need to discard
# But DIC is better with random host
# And is better than 2

# ----- Subgroup
MC_IUCN_SP_R_4c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2)):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_4c)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_4c)

MC_IUCN_SP_R_4d <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2)):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_4d)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_4d)

# DIC is increased with subgroup

# 5 ____________________________________________________________________________
# ----- Host
MC_IUCN_SP_R_5 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                               random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):HostCorrectedName,
                               data    = GMPD_IUCN_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_IUCN_SP_R_5)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_5)

# ill

MC_IUCN_SP_R_5a <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):HostCorrectedName,
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_5a)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_5a)

# Effective sample size is small for random effects so may need to discard
# But DIC is better with random host
# And is better than 2

# ----- Subgroup
MC_IUCN_SP_R_5c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_5c)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_5c)

# ill

MC_IUCN_SP_R_5d <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_IUCN_SP_R_5d)
MC_IUCN_SP <- append(MC_IUCN_SP, MC_IUCN_SP_R_5d)

# DIC is increased with subgroup




# Overall from MCMC, accept Latitude, MedianProp and quadratic term with host as a random effect
# Debatable whether to include subgroup

write_rds(MC_IUCN_SP, here::here("Data/Model back ups/MC_IUCN_SP.rds"))


# GBIF species #################################################################

## Prepping data ###############################################################
GMPD_GBIF_Species <- GMPD_Distances_Data %>%
  dplyr::filter(CleanAll, 
                RangeMethod   == "gbif",
                RangeTaxonLvl == "species") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared),
         ParasiteDetected       = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected     = HostsSampled - ParasiteDetected,
         HostSubgroup           = paste(HostCorrectedName, subgroup))

## GLMER #######################################################################
LM_GBIF_SP <- list()

### Fixed effects only #########################################################

# 1 ____________________________________________________________________________
LM_GBIF_SP_F_1 <- glm(formula = Prevalence ~ LatitudeScaled,
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_F_1) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_F_1)

# Accept Latitude

# 2 ____________________________________________________________________________
LM_GBIF_SP_F_2 <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled,
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_F_2) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_F_2)

# Accept MedianProp

# 3 ____________________________________________________________________________
LM_GBIF_SP_F_3 <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropScaled,
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_F_3) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_F_3)

# AIC is increased and interaction term is non-significant
# Reject interaction between Latitude and MedianProp

# 4 ____________________________________________________________________________
LM_GBIF_SP_F_4 <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2),
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_F_4) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_F_4)

# Reject AboveMedn in favour of quadratic based on AIC

# 5 ____________________________________________________________________________
LM_GBIF_SP_F_5 <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_F_5) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_F_5)

# Accept quadratic over AboveMedn based on AIC 



### Random effects #############################################################

# 1 ____________________________________________________________________________
# ----- Host
LM_GBIF_SP_R_1 <- glmer(formula = Prevalence ~ LatitudeScaled +
                          (1 + LatitudeScaled|HostCorrectedName),
                        data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_1) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_1)

# ----- Subgroup
LM_GBIF_SP_R_1c <- glmer(formula = Prevalence ~ LatitudeScaled +
                           (1 + LatitudeScaled|HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_1c) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_1c)

# Accept both host and subgroup

# 2 ____________________________________________________________________________
# ----- Host
LM_GBIF_SP_R_2 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                        (1 + LatitudeScaled + MedianPropScaled|HostCorrectedName),
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_2)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_2)

# singular fit

LM_GBIF_SP_R_2a <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                           (1 + LatitudeScaled + MedianPropScaled||HostCorrectedName),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_2a) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_2a)

# Ideally would model full vcv but this will do for now

# ----- Subgroup
LM_GBIF_SP_R_2c <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                           (1 + LatitudeScaled + MedianPropScaled|HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_2c)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_2c)

# failed to converge

LM_GBIF_SP_R_2d <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                           (1 + LatitudeScaled + MedianPropScaled||HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_2d) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_2d)

# AIC is reduced by adding both host and subgroup
# Accept host and subgroup as random effects

# 3 ____________________________________________________________________________
# Skipping 3 because it was rejected in the fixed effects

# 4 ____________________________________________________________________________
# ----- Host
LM_GBIF_SP_R_4 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                        (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4)

# failed to converge

LM_GBIF_SP_R_4a <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2)||HostCorrectedName),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4a)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4a)

# failed to converge

LM_GBIF_SP_R_4b <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                          (1 + LatitudeScaled + poly(MedianPropSquared, 1)|HostCorrectedName),
                        data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4b)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4b)

# singular fit

LM_GBIF_SP_R_4bb <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1)||HostCorrectedName),
                          data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4bb) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4bb)

# 4 would be nice, but 4bb is the best we can do
# based on AIC it's better to ignore quadratic term for now

# ----- Subgroup
LM_GBIF_SP_R_4c <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName:HostSubgroup),
                          data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4c)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4c)

# failed to converge

LM_GBIF_SP_R_4d <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2)||HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4d)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4d)

# failed to converge

LM_GBIF_SP_R_4e <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 1)|HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4e)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4e)

# failed to converge

LM_GBIF_SP_R_4ee <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1)||HostCorrectedName:HostSubgroup),
                          data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_4ee) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_4ee)

# again, based on AIC it's better to ignore quadratic term for now
# But accept host and subgroup

# 5 ____________________________________________________________________________
# ----- Host
LM_GBIF_SP_R_5 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                        (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn|HostCorrectedName),
                      data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5)

# failed to converge

LM_GBIF_SP_R_5a <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn||HostCorrectedName),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5a)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5a)

# failed to converge

LM_GBIF_SP_R_5b <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                          (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn|HostCorrectedName),
                        data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5b)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5b)

# singular fit

LM_GBIF_SP_R_5bb <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn||HostCorrectedName),
                          data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5bb) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5bb)

# 5 would be nice, but 5bb is the best we can do
# based on AIC it's better to ignore quadratic term for now

# ----- Subgroup
LM_GBIF_SP_R_5c <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn|HostCorrectedName:HostSubgroup),
                          data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5c)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5c)

# failed to converge

LM_GBIF_SP_R_5d <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn||HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5d)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5d)

# failed to converge

LM_GBIF_SP_R_5e <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                           (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn|HostCorrectedName:HostSubgroup),
                         data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5e)
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5e)

# failed to converge

LM_GBIF_SP_R_5ee <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn +
                            (1 + LatitudeScaled + poly(MedianPropSquared, 1):AboveMedn||HostCorrectedName:HostSubgroup),
                          data = GMPD_GBIF_Species, family = binomial)
summary(LM_GBIF_SP_R_5ee) 
LM_GBIF_SP <- append(LM_GBIF_SP, LM_GBIF_SP_R_5ee)

# again, based on AIC it's better to ignore quadratic term for now
# But accept host and subgroup

# from linear models we would accept Latitude and MedianProp, with host and subgroup as nested random effects
# But have to leave out the quadratic term in order to include random effects
write_rds(LM_GBIF_SP, here::here("Data/Model back ups/LM_GBIF_SP.rds"))

## MCMC ########################################################################
MC_GBIF_SP <- list()

### Fixed effects only #########################################################

# 1 ____________________________________________________________________________
MC_GBIF_SP_F_1 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_F_1)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_F_1)

# Accept latitude

# 2 ____________________________________________________________________________
MC_GBIF_SP_F_2 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_F_2)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_F_2)

# Accept MedianProp

# 3 ____________________________________________________________________________
MC_GBIF_SP_F_3 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled * MedianPropScaled,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_F_3)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_F_3)

# Could accept interaction between Latitude and MedianProp based on "significance"
# But the DIC is lower without the interaction term which agrees with LM approach
# Reject interaction between Latitude and MedianProp

# 4 ____________________________________________________________________________
MC_GBIF_SP_F_4 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_F_4)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_F_4)

# Reject AboveMedn in favour of quadratic based on AIC

# 5 ____________________________________________________________________________
MC_GBIF_SP_F_5 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_F_5)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_F_5)

# Accept quadratic over AboveMedn based on AIC

# Fixed effect results came out roughly the same for MC as for LM

### Random effects #############################################################

# 1 ____________________________________________________________________________
# ----- Host
MC_GBIF_SP_R_1 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                               random = ~ us(1 + LatitudeScaled):HostCorrectedName,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_R_1)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_1)

# Accept host

# ----- Subgroup
MC_GBIF_SP_R_1c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                                random = ~ us(1 + LatitudeScaled):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_1c)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_1c)

# DIC increases when I add subgroup in
# Accept host as a random effect but not subgroup?

# 2 ____________________________________________________________________________
# ----- Host
MC_GBIF_SP_R_2 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                               random = ~ us(1 + LatitudeScaled + MedianPropScaled):HostCorrectedName,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_R_2)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_2)

MC_GBIF_SP_R_2a <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                random = ~ idh(1 + LatitudeScaled + MedianPropScaled):HostCorrectedName,
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_2a)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_2a)

# ----- Subgroup
MC_GBIF_SP_R_2c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                random = ~ us(1 + LatitudeScaled + MedianPropScaled):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_2c)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_2c)

# ill

MC_GBIF_SP_R_2d <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                random = ~ idh(1 + LatitudeScaled + MedianPropScaled):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_2d)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_2d)

# DIC of host vs host and subgroup are identical so hard to say which is best random structure
# but definitely include host

# 3 ____________________________________________________________________________
# Skipping because interaction was rejected in fixed effects

# 4 ____________________________________________________________________________
# ----- Host
MC_GBIF_SP_R_4 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                               random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2)):HostCorrectedName,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_R_4)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_4)

# ill

MC_GBIF_SP_R_4a <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2)):HostCorrectedName,
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_4a)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_4a)

# Effective sample size is small for random effects so may need to discard
# But DIC is better with random host
# And is better than 2

# ----- Subgroup
MC_GBIF_SP_R_4c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2)):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_4c)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_4c)

# ill

MC_GBIF_SP_R_4d <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2)):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_4d)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_4d)

# 5 ____________________________________________________________________________
# ----- Host
MC_GBIF_SP_R_5 <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                               random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):HostCorrectedName,
                               data    = GMPD_GBIF_Species,
                               family  = "multinomial2",
                               nitt    = 130000,
                               thin    = 100,
                               burnin  = 30000),
                      silent = TRUE)
summary(MC_GBIF_SP_R_5)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_5)

MC_GBIF_SP_R_5a <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):HostCorrectedName,
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_5a)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_5a)

# Effective sample size is small for random effects so may need to discard
# But DIC is better with random host
# And is better than 2

# ----- Subgroup
MC_GBIF_SP_R_5c <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                random = ~ us(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_5c)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_5c)

MC_GBIF_SP_R_5d <- try(MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                random = ~ idh(1 + LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn):(HostCorrectedName + HostSubgroup),
                                data    = GMPD_GBIF_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000),
                       silent = TRUE)
summary(MC_GBIF_SP_R_5d)
MC_GBIF_SP <- append(MC_GBIF_SP, MC_GBIF_SP_R_5d)

# DIC is increased with subgroup

# Overall from MCMC, accept Latitude, MedianProp and quadratic term with host as a random effect
# Debatable whether to include subgroup

write_rds(MC_GBIF_SP, here::here("Data/Model back ups/MC_GBIF_SP.rds"))
