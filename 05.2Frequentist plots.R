# Plotting frequentist analysis
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(ggplot2) 
library(here)
extrafont::loadfonts(device = "win")

# Fixed effects ####
summary(LM_IUCN_Species$Fixed$VeryFull)

## Latitude ####
LM_Fixed_Effects_Data <- GMPD_IUCN_Species %>% 
  dplyr::select(Latitude, MedianProp, MedianPropSquared, CorrectedDistance, Prevalence) %>% 
  dplyr::bind_cols(predict(poly(.$MedianPropSquared, 2), .$MedianPropSquared)) %>% 
  dplyr::rename(MedianProp1 = "1", MedianProp2 = "2") %>% 
  dplyr::bind_cols(predict(poly(.$CorrectedDistance, 2), .$CorrectedDistance)) %>% 
  dplyr::rename(CorrectedDistance1 = "1", CorrectedDistance2 = "2")

### Bigline ####

# medIntercept    <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["(Intercept)", "Estimate"]
# medLatitude     <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["LatitudeScaled", "Estimate"]
# 
# medMedian1      <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["poly(MedianPropSquared, 2)1", "Estimate"]
# medMedian2      <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["poly(MedianPropSquared, 2)2", "Estimate"]
# medMedianInt1   <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]
# medMedianInt2   <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]
# 
# medDistance1    <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["poly(CorrectedDistance, 2)1", "Estimate"]
# medDistance2    <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["poly(CorrectedDistance, 2)2", "Estimate"]
# medDistanceInt1 <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]
# medDistanceInt2 <- summary(LM_IUCN_Species$Fixed$VeryFull)$Coefficients["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"]

LM_Fixed_Effects_Coeff <- summary(LM_IUCN_Species$Fixed$VeryFull)$coefficients %>% 
  as.data.frame() %>% 
  mutate(CI = 2 * `Std. Error`,
         Upper = Estimate + CI,
         Lower = Estimate - CI)

meanLat   <- mean(abs(LM_Fixed_Effects_Data$Latitude))
midLat    <- 40
lowLat    <- 0
highLat   <- 70

meanMed   <- mean(LM_Fixed_Effects_Data$MedianPropSquared)
meanMed1  <- mean(LM_Fixed_Effects_Data$MedianProp1)
meanMed2  <- mean(LM_Fixed_Effects_Data$MedianProp2)

meanDist  <- mean(LM_Fixed_Effects_Data$CorrectedDistance)
meanDist1 <- mean(LM_Fixed_Effects_Data$CorrectedDistance1)
meanDist2 <- mean(LM_Fixed_Effects_Data$CorrectedDistance2)

Prevalence = Intercept + 
  Latitude * LatSlope +
  
  Median1 * Med1Slope + 
  Median2 * Med2Slope +
  Distance1 * Dist1Slope +
  Distance2 * Dist2Slope +
  
  Latitude * Median1 * MedInteraction1 +
  Latitude * Median2 * MedInteraction2 +
  Latitude * Distance1 * DistInteraction1 +
  Latitude * Distance2 * DistInteraction2

# https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
LM_Fixed_Effects_Lines <- data.frame(LatitudeRange    = seq(min(abs(LM_Fixed_Effects_Data$Latitude)), 
                                                            max(abs(LM_Fixed_Effects_Data$Latitude)), 
                                                            length.out = 2000),
                                     
                                     MedianPropRange  = seq(min(LM_Fixed_Effects_Data$MedianPropSquared), 
                                                            max(LM_Fixed_Effects_Data$MedianPropSquared), 
                                                            length.out = 2000),
                                     MedianPropRange1 = seq(min(LM_Fixed_Effects_Data$MedianProp1), 
                                                            max(LM_Fixed_Effects_Data$MedianProp1), 
                                                            length.out = 2000),
                                     MedianPropRange2 = seq(min(LM_Fixed_Effects_Data$MedianProp2), 
                                                            max(LM_Fixed_Effects_Data$MedianProp2), 
                                                            length.out = 2000),
                                     
                                     DistanceRange    = seq(min(LM_Fixed_Effects_Data$CorrectedDistance), 
                                                            max(LM_Fixed_Effects_Data$CorrectedDistance), 
                                                            length.out = 2000),
                                     DistanceRange1   = seq(min(LM_Fixed_Effects_Data$CorrectedDistance1), 
                                                            max(LM_Fixed_Effects_Data$CorrectedDistance1), 
                                                            length.out = 2000),
                                     DistanceRange2   = seq(min(LM_Fixed_Effects_Data$CorrectedDistance2), 
                                                            max(LM_Fixed_Effects_Data$CorrectedDistance2), 
                                                            length.out = 2000)) %>% 
  
  dplyr::mutate(LatitudeRangeScaled   = (LatitudeRange   - meanLat)/(sd(abs(LM_Fixed_Effects_Data$Latitude) - meanLat)),
                MedianPropRangeScaled = (MedianPropRange - meanMed)/(sd(LM_Fixed_Effects_Data$MedianPropSquared - meanMed)),
                DistanceRangeScaled   = (DistanceRange   - meanDist)/(sd(LM_Fixed_Effects_Data$CorrectedDistance - meanDist))) %>%
  
  dplyr::mutate(logitLat = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                              (LatitudeRangeScaled * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                              (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                              (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                              
                              (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                              (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                              
                              (LatitudeRangeScaled * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                              (LatitudeRangeScaled * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                              
                              (LatitudeRangeScaled * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                              (LatitudeRangeScaled * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitLatUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                   (LatitudeRangeScaled * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                   (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                   (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                   
                                   (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                   (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                   
                                   (LatitudeRangeScaled * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                   (LatitudeRangeScaled * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                   
                                   (LatitudeRangeScaled * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                   (LatitudeRangeScaled * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitLatLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                   (LatitudeRangeScaled * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                   (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                   (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                   
                                   (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                   (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                   
                                   (LatitudeRangeScaled * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                   (LatitudeRangeScaled * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                   
                                   (LatitudeRangeScaled * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                   (LatitudeRangeScaled * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________________________
                
                logitMed = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                              (meanLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                              (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                              (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                              
                              (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                              (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                              
                              (meanLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                              (meanLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                              
                              (meanLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                              (meanLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitMedUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                   (meanLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                   (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                   (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                   
                                   (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                   (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                   
                                   (meanLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                   (meanLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                   
                                   (meanLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                   (meanLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitMedLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                   (meanLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                   (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                   (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                   
                                   (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                   (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                   
                                   (meanLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                   (meanLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                   
                                   (meanLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                   (meanLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________________________
                
                logitMedLow = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                                 (lowLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                 (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                 (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                 
                                 (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                                 (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                                 
                                 (lowLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                                 (lowLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                                 
                                 (lowLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                                 (lowLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),

                logitMedLowUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                      (lowLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                      (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                      (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                      
                                      (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                      (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                      
                                      (lowLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                      (lowLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                      
                                      (lowLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                      (lowLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitMedLowLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                      (lowLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                      (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                      (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                      
                                      (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                      (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                      
                                      (lowLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                      (lowLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                      
                                      (lowLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                      (lowLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________
                
                logitMedMid = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                                 (midLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                 (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                 (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                 
                                 (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                                 (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                                 
                                 (midLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                                 (midLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                                 
                                 (midLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                                 (midLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitMedMidUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                      (midLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                      (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                      (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                      
                                      (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                      (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                      
                                      (midLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                      (midLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                      
                                      (midLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                      (midLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitMedMidLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                      (midLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                      (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                      (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                      
                                      (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                      (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                      
                                      (midLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                      (midLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                      
                                      (midLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                      (midLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________________________
                
                logitMedHigh = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                                  (highLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                  (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                  (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                  
                                  (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                                  (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                                  
                                  (highLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                                  (highLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                                  
                                  (highLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                                  (highLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),

                logitMedHighUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                       (highLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                       (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                       (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                       
                                       (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                       (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                       
                                       (highLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                       (highLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                       
                                       (highLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                       (highLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitMedHighLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                       (highLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                       (MedianPropRange1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                       (MedianPropRange2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                       
                                       (meanDist1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                       (meanDist2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                       
                                       (highLat * MedianPropRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                       (highLat * MedianPropRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                       
                                       (highLat * meanDist1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                       (highLat * meanDist2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________
                                
                logitDist = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                               (meanLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                               (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                               (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                               
                               (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                               (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                               
                               (meanLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                               (meanLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                               
                               (meanLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                               (meanLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitDistUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                    (meanLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                    (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                    (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                    
                                    (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                    (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                    
                                    (meanLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                    (meanLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                    
                                    (meanLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                    (meanLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitDistLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                    (meanLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                    (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                    (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                    
                                    (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                    (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                    
                                    (meanLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                    (meanLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                    
                                    (meanLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                    (meanLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________________________                
                
                logitDistLow = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                                  (lowLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                  (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                  (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                  
                                  (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                                  (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                                  
                                  (lowLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                                  (lowLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                                  
                                  (lowLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                                  (lowLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitDistLowUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                       (lowLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                       (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                       (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                       
                                       (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                       (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                       
                                       (lowLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                       (lowLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                       
                                       (lowLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                       (lowLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitDistLowLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                       (lowLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                       (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                       (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                       
                                       (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                       (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                       
                                       (lowLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                       (lowLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                       
                                       (lowLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                       (lowLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________________________
                
                logitDistMid = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                                  (midLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                  (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                  (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                  
                                  (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                                  (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                                  
                                  (midLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                                  (midLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                                  
                                  (midLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                                  (midLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitDistMidUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                       (midLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                       (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                       (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                       
                                       (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                       (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                       
                                       (midLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                       (midLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                       
                                       (midLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                       (midLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitDistMidLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                       (midLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                       (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                       (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                       
                                       (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                       (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                       
                                       (midLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                       (midLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                       
                                       (midLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                       (midLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"])),
                
                #_______________________________________________________________________________
                
                logitDistHigh = (LM_Fixed_Effects_Coeff["(Intercept)", "Estimate"] +
                                   (highLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                   (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                   (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                   
                                   (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Estimate"]) +
                                   (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Estimate"]) +
                                   
                                   (highLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Estimate"]) +
                                   (highLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Estimate"]) +
                                   
                                   (highLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Estimate"]) +
                                   (highLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Estimate"])),
                
                logitDistHighUpper = (LM_Fixed_Effects_Coeff["(Intercept)", "Upper"] +
                                        (highLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                        (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                        (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                        
                                        (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Upper"]) +
                                        (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Upper"]) +
                                        
                                        (highLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Upper"]) +
                                        (highLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Upper"]) +
                                        
                                        (highLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Upper"]) +
                                        (highLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Upper"])),
                
                logitDistHighLower = (LM_Fixed_Effects_Coeff["(Intercept)", "Lower"] +
                                        (highLat * LM_Fixed_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                        (meanMed1 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                        (meanMed2 * LM_Fixed_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                        
                                        (DistanceRange1 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)1", "Lower"]) +
                                        (DistanceRange2 * LM_Fixed_Effects_Coeff["poly(CorrectedDistance, 2)2", "Lower"]) +
                                        
                                        (highLat * meanMed1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)1", "Lower"]) +
                                        (highLat * meanMed2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(MedianPropSquared, 2)2", "Lower"]) +
                                        
                                        (highLat * DistanceRange1 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)1", "Lower"]) +
                                        (highLat * DistanceRange2 * LM_Fixed_Effects_Coeff["LatitudeScaled:poly(CorrectedDistance, 2)2", "Lower"]))
  ) %>% 
  dplyr::mutate(across(starts_with("logit"), 
                       function(x) {1/(1+exp(-x))}, 
                       .names = "{col}_line")) %>% 
  dplyr::rename_with(~str_replace(., pattern = "logit", replacement = ""), 
                     ends_with("_line"))

# geom_rug(aes(y = visited, colour = lvisited), data = wasp)
### Latitude ####
ggplot(LM_Fixed_Effects_Lines, aes(LatitudeRange, Lat_line)) +
  geom_ribbon(aes(ymin = LatLower_line, ymax = LatUpper_line),
              alpha = 0.1) +
  geom_line(col = "#DE8F6E", linewidth = 2, lineend = "round") +
  
  geom_rug(data = LM_Fixed_Effects_Data, aes(Latitude, Prevalence), sides = "b") +
  
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        panel.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
  ) +
  
  labs(x = "Latitude", y = "Parasitism rate") +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), limits = c(0,90)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) 


### MedianProp ####
ggplot(LM_Fixed_Effects_Lines, aes(MedianPropRange, MedMid_line)) +
  geom_ribbon(aes(ymin = MedMidLower_line, ymax = MedMidUpper_line),
              alpha = 0.1) +
  geom_ribbon(aes(ymin = MedLowLower_line, ymax = MedLowUpper_line),
              alpha = 0.1) +
  geom_ribbon(aes(ymin = MedHighLower_line, ymax = MedHighUpper_line),
              alpha = 0.1) +
  
  geom_rug(data = LM_Fixed_Effects_Data, aes(MedianPropSquared, Prevalence), sides = "b") +
  
  geom_line(col = "#DE8F6E", linewidth = 2, lineend = "round") +
  geom_line(aes(y = MedLow_line), 
            col = "#88AB75", 
            linewidth = 2, 
            lineend = "round",
            linetype = "dashed") +
  geom_line(aes(y = MedHigh_line), 
            col = "#2D93AD", 
            linewidth = 2, 
            lineend = "round",
            linetype = "dashed") +
  
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        panel.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
  ) +
  
labs(x = "Range position", y = "Parasitism rate") +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1,1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) 

### CorrectedDistance ####
ggplot(LM_Fixed_Effects_Lines, aes(DistanceRange, Dist_line)) +
  geom_ribbon(aes(ymin = DistMidLower_line, ymax = DistMidUpper_line),
              alpha = 0.1) +
  geom_ribbon(aes(ymin = DistLowLower_line, ymax = DistLowUpper_line),
              alpha = 0.1) +
  geom_ribbon(aes(ymin = DistHighLower_line, ymax = DistHighUpper_line),
              alpha = 0.1) +
  
  geom_rug(data = LM_Fixed_Effects_Data, aes(CorrectedDistance, Prevalence), sides = "b") +
  
  geom_line(col = "#DE8F6E", linewidth = 2, lineend = "round") +
  geom_line(aes(y = DistLow_line), 
            col = "#88AB75", 
            linewidth = 2, 
            lineend = "round",
            linetype = "dashed") +
  geom_line(aes(y = DistHigh_line), col = "#2D93AD", linewidth = 2, lineend = "round",
            linetype = "dashed") +
  
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        panel.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
  ) +
  
labs(x = "Niche position", y = "Parasitism rate") +
  scale_x_continuous(limits = c(0, 36)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) 

  

