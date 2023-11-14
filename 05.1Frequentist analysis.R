# Frequentist stats
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# libraries and data ###########################################################
library(lme4)
library(here)
library(DHARMa)

# TODO: refactor with modelling function?

# Model naming convention ######################################################
## Model main lists:
# LM = Frequentist rather than Bayesian
# IUCN/GBIF = method used for restricting data and calculating position within range
# Species/Subgroup = position within range calculated from position within species (SP) or subgroup (SG) range
LM_IUCN_Species  <- list()
LM_IUCN_Subgroup <- list()
LM_GBIF_Species  <- list()
LM_GBIF_Subgroup <- list()

# TODO: update
## Model sub lists
# Fixed     = Fixed effects only
# Host      = Host random effects
# ParType   = ParType random effects
# Both      = Host and ParType random effects
# Phylo     = Host with phylogeny
# _W suffix = weighted by HostsSampled

LM_IUCN_Species$Fixed    <- list()
LM_IUCN_Species$Fixed_W  <- list()
LM_IUCN_Species$Fixed_ID <- list()

# Host focus
LM_IUCN_Species$Host          <- list()
LM_IUCN_Species$Host_IO       <- list()
LM_IUCN_Species$Host_W        <- list()
LM_IUCN_Species$Host_WIO      <- list()

LM_IUCN_Species$HostGroup     <- list()
LM_IUCN_Species$HostGroup_IO  <- list()
LM_IUCN_Species$HostGroup_W   <- list()
LM_IUCN_Species$HostGroup_WIO <- list()

# LM_IUCN_Species$Group         <- list()
# LM_IUCN_Species$Group_IO      <- list()
# LM_IUCN_Species$Group_W       <- list()
# LM_IUCN_Species$Group_WIO     <- list()

# Parasite focus
LM_IUCN_Species$Parasite     <- list()
LM_IUCN_Species$Parasite_IO  <- list()
LM_IUCN_Species$Parasite_W   <- list()
LM_IUCN_Species$Parasite_WIO <- list()

LM_IUCN_Species$ParType      <- list()
LM_IUCN_Species$ParType_IO   <- list()
LM_IUCN_Species$ParType_W    <- list()
LM_IUCN_Species$ParType_WIO  <- list()

# LM_IUCN_Species$ParFixed     <- list()
# LM_IUCN_Species$ParFixed_IO  <- list()
# LM_IUCN_Species$ParFixed_W   <- list()
# LM_IUCN_Species$ParFixed_WIO <- list()

# Both
LM_IUCN_Species$Both       <- list()
LM_IUCN_Species$Both_IO    <- list()
LM_IUCN_Species$Both_W     <- list()
LM_IUCN_Species$Both_WIO   <- list()

LM_IUCN_Species$Nested     <- list()
LM_IUCN_Species$Nested_IO  <- list()
LM_IUCN_Species$Nested_W   <- list()
LM_IUCN_Species$Nested_WIO <- list()

LM_IUCN_Species$Combo      <- list()
LM_IUCN_Species$Combo_IO   <- list()
LM_IUCN_Species$Combo_W    <- list()
LM_IUCN_Species$Combo_WIO  <- list()

#___________________________________
LM_IUCN_Subgroup$Fixed     <- list()
LM_IUCN_Subgroup$Host      <- list()
LM_IUCN_Subgroup$ParType   <- list()
LM_IUCN_Subgroup$Both      <- list()
LM_IUCN_Subgroup$Phylo     <- list()
  
LM_GBIF_Species$Fixed      <- list()
LM_GBIF_Species$Host       <- list()
LM_GBIF_Species$ParType    <- list()
LM_GBIF_Species$Both       <- list()
LM_GBIF_Species$Phylo      <- list()
  
LM_GBIF_Subgroup$Fixed     <- list()
LM_GBIF_Subgroup$Host      <- list()
LM_GBIF_Subgroup$ParType   <- list()
LM_GBIF_Subgroup$Both      <- list()
LM_GBIF_Subgroup$Phylo     <- list()

# IUCN with Species ############################################################
#///////////////////////////////////////////////////////////////////////////////
## Geographic niche ############################################################
#///////////////////////////////////////////////////////////////////////////////
### Null #######################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Null       <- glm(formula = Prevalence ~ 1, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Null     <- glm(formula = Prevalence ~ 1, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Fixed_Grp$LatQMed    <- glm(formula = Prevalence ~ Group, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_Grp_W$LatQMed  <- glm(formula = Prevalence ~ Group, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_Grp_ID$LatQMed <- glm(formula = Prevalence ~ Group, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Fixed_Typ$LatQMed    <- glm(formula = Prevalence ~ ParType, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_Typ_W$LatQMed  <- glm(formula = Prevalence ~ ParType, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_Typ_ID$LatQMed <- glm(formula = Prevalence ~ ParType, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

# Host focus
#_______________________________________________________________________________
LM_IUCN_Species$Host$Null        <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$Host_IO$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Host_W$Null      <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$Host_WIO$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$HostGroup$Null       <- lme4::glmer(formula = Prevalence ~ 1 + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$HostGroup_IO$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$HostGroup_W$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$HostGroup_WIO$Null   <- lme4::glmer(formula = Prevalence ~ 1 + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Group$Null       <- lme4::glmer(formula = Prevalence ~ Group + (1 + Group|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Group_IO$Null    <- lme4::glmer(formula = Prevalence ~ Group + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Group_W$Null     <- lme4::glmer(formula = Prevalence ~ Group + (1 + Group|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Group_WIO$Null   <- lme4::glmer(formula = Prevalence ~ Group + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

# Parasite focus
#_______________________________________________________________________________
LM_IUCN_Species$Parasite$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$Parasite_IO$Null  <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Parasite_W$Null   <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$Parasite_WIO$Null <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Type$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$Type_IO$Null  <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Type_W$Null   <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$Type_WIO$Null <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$ParType$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$ParType_IO$Null  <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$ParType_W$Null   <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$ParType_WIO$Null <- lme4::glmer(formula = Prevalence ~ 1 + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$TypePar$Null     <- lme4::glmer(formula = Prevalence ~ ParType + (1 + ParType|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$TypePar_IO$Null  <- lme4::glmer(formula = Prevalence ~ ParType + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$TypePar_W$Null   <- lme4::glmer(formula = Prevalence ~ ParType + (1 + ParType|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$TypePar_WIO$Null <- lme4::glmer(formula = Prevalence ~ ParType + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

# Both
#_______________________________________________________________________________
LM_IUCN_Species$Both$Null        <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$Both_IO$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Both_W$Null      <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$Both_WIO$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$BothSpecies$Null        <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$BothSpecies_IO$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$BothSpecies_W$Null      <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$BothSpecies_WIO$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Nested$Null      <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$Nested_IO$Null   <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Nested_W$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$Nested_WIO$Null  <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Combo$Null       <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
# LM_IUCN_Species$Combo_IO$Null    <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Combo_W$Null     <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
# LM_IUCN_Species$Combo_WIO$Null   <- lme4::glmer(formula = Prevalence ~ 1 + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude ###################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Lat       <- glm(formula = Prevalence ~ LatitudeScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Lat     <- glm(formula = Prevalence ~ LatitudeScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Lat    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### MedianProp #################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Med       <- glm(formula = Prevalence ~ MedianPropSquScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Med     <- glm(formula = Prevalence ~ MedianPropSquScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Med    <- lme4::glmer(formula = Prevalence ~ MedianPropSquScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Quadratic MedianProp #######################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QMed       <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QMed     <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QMed    <- lme4::glmer(formula = Prevalence ~ poly(MedianPropSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Asymmetric Quadratic MedianProp ############################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QMedAsym       <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QMedAsym     <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QMedAsym    <- lme4::glmer(formula = Prevalence ~ poly(MedianPropSquared, 2) * AboveMedn + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude + MedianProp ######################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatMed       <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropSquScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatMed     <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropSquScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatMed    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropSquScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude + Asymmetric MedianProp ###########################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatMedAsym       <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropSquScaled * AboveMedn, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatMedAsym     <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropSquScaled * AboveMedn, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatMedAsym    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropSquScaled * AboveMedn + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude + Quadratic MedianProp ###########################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMed       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMed     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMed    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Fixed_Grp$LatQMed    <- glm(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_Grp_W$LatQMed  <- glm(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_Grp_ID$LatQMed <- glm(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Fixed_Typ$LatQMed    <- glm(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_Typ_W$LatQMed  <- glm(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_Typ_ID$LatQMed <- glm(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

# Host focus
#_______________________________________________________________________________
LM_IUCN_Species$Host$LatQMed        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Host_IO$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Host_W$LatQMed      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Host_WIO$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$HostGroup$LatQMed       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$HostGroup_IO$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$HostGroup_W$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$HostGroup_WIO$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Group$LatQMed       <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + (1 + Group + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Group_IO$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Group_W$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + (1 + Group + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Group_WIO$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

# Parasite focus
#_______________________________________________________________________________
LM_IUCN_Species$Parasite$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Parasite_IO$LatQMed  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Parasite_W$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Parasite_WIO$LatQMed <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Type$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Type_IO$LatQMed  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Type_W$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Type_WIO$LatQMed <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$ParType$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$ParType_IO$LatQMed  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$ParType_W$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$ParType_WIO$LatQMed <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$TypePar$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + (1 + ParType + LatitudeScaled + poly(MedianPropSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$TypePar_IO$LatQMed  <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$TypePar_W$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + (1 + ParType + LatitudeScaled + poly(MedianPropSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$TypePar_WIO$LatQMed <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

# Both
#_______________________________________________________________________________
LM_IUCN_Species$Both$LatQMed        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Both_IO$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Both_W$LatQMed      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Both_WIO$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$BothSpecies$LatQMed        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$BothSpecies_IO$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$BothSpecies_W$LatQMed      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$BothSpecies_WIO$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Nested$LatQMed      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Nested_IO$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Nested_W$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Nested_WIO$LatQMed  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Combo$LatQMed       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Combo_IO$LatQMed    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Combo_W$LatQMed     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Combo_WIO$LatQMed   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

#///////////////////////////////////////////////////////////////////////////////
### Latitude + Asymmetric Quadratic MedianProp #################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsym       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsym     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsym    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude * MedianProp ######################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatMedInt       <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropSquScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatMedInt     <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropSquScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatMedInt    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled * MedianPropSquScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude * Asymmetric MedianProp ######################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatMedIntAsym       <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropSquScaled * AboveMedn, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatMedIntAsym     <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropSquScaled * AboveMedn, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatMedIntAsym    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled * MedianPropSquScaled * AboveMedn + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude * Quadratic MedianProp ######################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedInt       <- glm(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedInt     <- glm(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedInt    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Latitude * Asymmetric Quadratic MedianProp ######################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedIntAsym       <- glm(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedIntAsym     <- glm(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedIntAsym    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2) * AboveMedn + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
## Climatic niche ##############################################################
#///////////////////////////////////////////////////////////////////////////////
### PCA axes ###################################################################
#///////////////////////////////////////////////////////////////////////////////
#### Axis1 #####################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Axis1       <- glm(formula = Prevalence ~ Axis1Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Axis1     <- glm(formula = Prevalence ~ Axis1Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Axis1    <- lme4::glmer(formula = Prevalence ~ Axis1Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Axis2 #####################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Axis2       <- glm(formula = Prevalence ~ Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Axis2     <- glm(formula = Prevalence ~ Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Axis2    <- lme4::glmer(formula = Prevalence ~ Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Axes ######################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Axes       <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Axes     <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Axes    <- lme4::glmer(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Axis1 ##########################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatAxis1       <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled , data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatAxis1     <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled , data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatAxis1    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + Axis1Scaled  + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Axis2Scaled ####################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatAxis2       <- glm(formula = Prevalence ~ LatitudeScaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatAxis2     <- glm(formula = Prevalence ~ LatitudeScaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatAxis2    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Axes ###########################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatAxes       <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatAxes     <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatAxes    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + Axis1Scaled + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Quadratic MedianProp + Axis1 ##############################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QMedAxis1       <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QMedAxis1     <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QMedAxis1    <- lme4::glmer(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Quadratic MedianProp + Axis2Scaled ########################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QMedAxis2       <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QMedAxis2     <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QMedAxis2    <- lme4::glmer(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Quadratic MedianProp + Axes ###############################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QMedAxes       <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QMedAxes     <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QMedAxes    <- lme4::glmer(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Axis1 ###################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAxis1       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAxis1     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAxis1    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Axis1 ###################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymAxis1       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis1Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymAxis1     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis1Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis1    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis1Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Axis2 #############################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAxis2       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAxis2     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAxis2    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Axis2 #############################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymAxis2       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymAxis2     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymAxis2    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Axes ####################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAxes       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAxes     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAxes    <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Axes ####################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymAxes     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymAxes   <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis1Scaled + Axis2Scaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymAxes  <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + Axis1Scaled + Axis2Scaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
### Kernel Density #############################################################
#///////////////////////////////////////////////////////////////////////////////
#### Squared Distance ##########################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Dist                 <- glm(formula = Prevalence ~ CorrectedDistanceSquScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Dist               <- glm(formula = Prevalence ~ CorrectedDistanceSquScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Dist              <- lme4::glmer(formula = Prevalence ~ CorrectedDistanceSquScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Quadratic Squared Distance ################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QDist                <- glm(formula = Prevalence ~ poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QDist              <- glm(formula = Prevalence ~ poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QDist             <- lme4::glmer(formula = Prevalence ~ poly(CorrectedDistanceSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Density ###################################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$Dens                 <- glm(formula = Prevalence ~ CorrectedDensityScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$Dens               <- glm(formula = Prevalence ~ CorrectedDensityScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$Dens              <- lme4::glmer(formula = Prevalence ~ CorrectedDensityScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Quadratic Density #########################################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$QDens                <- glm(formula = Prevalence ~ poly(CorrectedDensity, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$QDens              <- glm(formula = Prevalence ~ poly(CorrectedDensity, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$QDens             <- lme4::glmer(formula = Prevalence ~ poly(CorrectedDensity, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Distance ################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedDist          <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceSquScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedDist        <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceSquScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedDist       <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceSquScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Distance ################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymDist      <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + CorrectedDistanceSquScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymDist    <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + CorrectedDistanceSquScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymDist   <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + CorrectedDistanceSquScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Density #################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedDens          <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDensityScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedDens        <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDensityScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedDens       <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDensityScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Density #################################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymDens      <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + CorrectedDensityScaled, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymDens    <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + CorrectedDensityScaled, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymDens   <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + CorrectedDensityScaled + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Quadratic Distance ######################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedQDist        <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedQDist      <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedQDist     <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Fixed_Grp$LatQMedQDist    <- glm(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_Grp_W$LatQMedQDist  <- glm(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_Grp_ID$LatQMedQDist <- lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

LM_IUCN_Species$Fixed_Typ$LatQMedQDist    <- glm(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_Typ_W$LatQMedQDist  <- glm(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_Typ_ID$LatQMedQDist <- lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

# Host focus
#_______________________________________________________________________________
LM_IUCN_Species$Host$LatQMedQDist          <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Host_IO$LatQMedQDist       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Host_W$LatQMedQDist        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Host_WIO$LatQMedQDist      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$HostGroup$LatQMedQDist     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$HostGroup_IO$LatQMedQDist  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$HostGroup_W$LatQMedQDist   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$HostGroup_WIO$LatQMedQDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|Group/HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Group$LatQMedQDist         <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Group_IO$LatQMedQDist      <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Group_W$LatQMedQDist       <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Group_WIO$LatQMedQDist     <- try(lme4::glmer(formula = Prevalence ~ Group + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

# Parasite focus
#_______________________________________________________________________________
LM_IUCN_Species$Parasite$LatQMedQDist     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Parasite_IO$LatQMedQDist  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Parasite_W$LatQMedQDist   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Parasite_WIO$LatQMedQDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Type$LatQMedQDist         <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Type_IO$LatQMedQDist      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Type_W$LatQMedQDist       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Type_WIO$LatQMedQDist     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$ParType$LatQMedQDist      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$ParType_IO$LatQMedQDist   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$ParType_W$LatQMedQDist    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$ParType_WIO$LatQMedQDist  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParType/ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$TypePar$LatQMedQDist      <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$TypePar_IO$LatQMedQDist   <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$TypePar_W$LatQMedQDist    <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$TypePar_WIO$LatQMedQDist  <- try(lme4::glmer(formula = Prevalence ~ ParType + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

# Both
#_______________________________________________________________________________
LM_IUCN_Species$Both$LatQMedQDist            <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Both_IO$LatQMedQDist         <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Both_W$LatQMedQDist          <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Both_WIO$LatQMedQDist        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName) + (1|ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$BothSpecies$LatQMedQDist     <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$BothSpecies_IO$LatQMedQDist  <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$BothSpecies_W$LatQMedQDist   <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$BothSpecies_WIO$LatQMedQDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName) + (1|ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Nested$LatQMedQDist          <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Nested_IO$LatQMedQDist       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Nested_W$LatQMedQDist        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Nested_WIO$LatQMedQDist      <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName:ParType), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

LM_IUCN_Species$Combo$LatQMedQDist           <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Combo_IO$LatQMedQDist        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial), silent = TRUE)
LM_IUCN_Species$Combo_W$LatQMedQDist         <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2)|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)
LM_IUCN_Species$Combo_WIO$LatQMedQDist       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceSquared, 2) + (1|HostCorrectedName:ParasiteCorrectedName), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled), silent = TRUE)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Quadratic Distance ######################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymQDist     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymQDist   <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + poly(CorrectedDistanceSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymQDist  <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + poly(CorrectedDistanceSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Quadratic MedianProp + Quadratic Density #######################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedQDens         <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDensity, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedQDens       <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDensity, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedQDens      <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDensity, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude + Asymmetric Quadratic MedianProp + Quadratic Density #######################
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymQDens     <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + poly(CorrectedDensity, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymQDens   <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + poly(CorrectedDensity, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymQDens  <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) * AboveMedn + poly(CorrectedDensity, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude * Quadratic Distance + Quadratic MedianProp ###########
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedQIDist        <- glm(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistanceSquared, 2) + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedQIDist      <- glm(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistanceSquared, 2) + poly(MedianPropSquared, 2), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedQIDist     <- lme4::glmer(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistanceSquared, 2) + poly(MedianPropSquared, 2) + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)

#///////////////////////////////////////////////////////////////////////////////
#### Latitude * Quadratic Distance + Asymmetric Quadratic MedianProp ###########
#///////////////////////////////////////////////////////////////////////////////
LM_IUCN_Species$Fixed$LatQMedAsymQIDist    <- glm(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistanceSquared, 2) + poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial)
LM_IUCN_Species$Fixed_W$LatQMedAsymQIDist  <- glm(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistanceSquared, 2) + poly(MedianPropSquared, 2) * AboveMedn, data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)
LM_IUCN_Species$Fixed_ID$LatQMedAsymQIDist <- lme4::glmer(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistanceSquared, 2) + poly(MedianPropSquared, 2) * AboveMedn + (1|RowID), data = GMPD_IUCN_Species, family = binomial, weights = HostsSampled)


# saveRDS(LM_IUCN_Species, here::here("Data/Model back ups/LM_IUCN_Species.rds"))


