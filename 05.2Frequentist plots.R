# Plotting frequentist analysis
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(ggplot2) 
library(here)
extrafont::loadfonts(device = "win")

# LatQMedQDist + Host:Parasite ####
summary(LM_IUCN_Species$Combo_WIO$LatQMedQDist)

### Bigline ####

LM_Final_Effects_Coeff <- summary(LM_IUCN_Species$Combo_WIO$LatQMedQDist)$coefficients %>% 
  as.data.frame() %>% 
  mutate(CI = 2 * `Std. Error`,
         Upper = Estimate + CI,
         Lower = Estimate - CI)

meanLat   <- mean(abs(GMPD_IUCN_Species$Latitude))
midLat    <- 40
lowLat    <- 0
highLat   <- 80

meanLat.sc <- (meanLat - meanLat)/(sd(abs(GMPD_IUCN_Species$Latitude) - meanLat))
midLat.sc  <- (midLat - meanLat)/(sd(abs(GMPD_IUCN_Species$Latitude) - meanLat))
lowLat.sc  <- (lowLat - meanLat)/(sd(abs(GMPD_IUCN_Species$Latitude) - meanLat))
highLat.sc <- (highLat - meanLat)/(sd(abs(GMPD_IUCN_Species$Latitude) - meanLat))

zeroMed1 <- predict(poly(GMPD_IUCN_Species$MedianPropSquared, 2), 0)[[1]]
zeroMed2 <- predict(poly(GMPD_IUCN_Species$MedianPropSquared, 2), 0)[[2]]

zeroDist1 <- predict(poly(GMPD_IUCN_Species$CorrectedDistanceSquared, 2), 0)[[1]]
zeroDist2 <- predict(poly(GMPD_IUCN_Species$CorrectedDistanceSquared, 2), 0)[[2]]


# Model equation:
# Prevalence = Intercept + 
#   Latitude * LatSlope +
#   
#   Median1 * Med1Slope + 
#   Median2 * Med2Slope +
#   
#   Distance1 * Dist1Slope +
#   Distance2 * Dist2Slope 

# https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
LM_Final_Effects_Lines <- data.frame(LatitudeRange    = seq(0, max(abs(GMPD_IUCN_Species$Latitude)), 
                                                            length.out = 2000),
                                     
                                     MedianPropRange  = seq(min(GMPD_IUCN_Species$MedianPropSquared), 
                                                            max(GMPD_IUCN_Species$MedianPropSquared), 
                                                            length.out = 2000),
                                     
                                     DistanceRange    = seq(min(GMPD_IUCN_Species$CorrectedDistanceSquared), 
                                                            max(GMPD_IUCN_Species$CorrectedDistanceSquared), 
                                                            length.out = 2000)) %>% 
  dplyr::bind_cols(predict(poly(GMPD_IUCN_Species$MedianPropSquared, 2), 
                           seq(min(GMPD_IUCN_Species$MedianPropSquared), 
                               max(GMPD_IUCN_Species$MedianPropSquared),
                               length.out = 2000))) %>% 
  dplyr::rename(MedianPropRange1 = "1", MedianPropRange2 = "2") %>% 
  dplyr::bind_cols(predict(poly(GMPD_IUCN_Species$CorrectedDistanceSquared, 2), 
                           seq(min(GMPD_IUCN_Species$CorrectedDistanceSquared), 
                               max(GMPD_IUCN_Species$CorrectedDistanceSquared),
                               length.out = 2000))) %>% 
  dplyr::rename(DistanceRange1 = "1", DistanceRange2 = "2") %>%   
  
  
  dplyr::mutate(LatitudeRangeScaled   = (LatitudeRange   - meanLat)/(sd(abs(GMPD_IUCN_Species$Latitude) - meanLat))) %>%
  
  dplyr::mutate(logitLat      = (LM_Final_Effects_Coeff["(Intercept)", "Estimate"] +
                                   (LatitudeRangeScaled * LM_Final_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                   
                                   (zeroMed1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                   (zeroMed2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                   
                                   (zeroDist1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Estimate"]) + 
                                   (zeroDist2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Estimate"])), 
                
                logitLatUpper = (LM_Final_Effects_Coeff["(Intercept)", "Upper"] +
                                   (LatitudeRangeScaled * LM_Final_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                   
                                   (zeroMed1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                   (zeroMed2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                   
                                   (zeroDist1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Upper"]) + 
                                   (zeroDist2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Upper"])),
                
                logitLatLower = (LM_Final_Effects_Coeff["(Intercept)", "Lower"] +
                                   (LatitudeRangeScaled * LM_Final_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                   
                                   (zeroMed1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                   (zeroMed2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                   
                                   (zeroDist1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Lower"]) + 
                                   (zeroDist2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Lower"])),
                
                #_______________________________________________________________________________
                
                logitMed      = (LM_Final_Effects_Coeff["(Intercept)", "Estimate"] +
                                   (midLat.sc * LM_Final_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                                   
                                   (MedianPropRange1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                                   (MedianPropRange2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                                   
                                   (zeroDist1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Estimate"]) +
                                   (zeroDist2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Estimate"])),
                
                logitMedUpper = (LM_Final_Effects_Coeff["(Intercept)", "Upper"] +
                                   (midLat.sc * LM_Final_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                   (MedianPropRange1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                   (MedianPropRange2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                   
                                   (zeroDist1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Upper"]) +
                                   (zeroDist2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Upper"])),
                
                logitMedLower = (LM_Final_Effects_Coeff["(Intercept)", "Lower"] +
                                   (midLat.sc * LM_Final_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                   (MedianPropRange1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                   (MedianPropRange2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                   
                                   (zeroDist1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Lower"]) +
                                   (zeroDist2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Lower"])),
                
                #_______________________________________________________________
                
                logitDist = (LM_Final_Effects_Coeff["(Intercept)", "Estimate"] +
                               (midLat.sc * LM_Final_Effects_Coeff["LatitudeScaled", "Estimate"]) +
                               (zeroMed1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Estimate"]) +
                               (zeroMed2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Estimate"]) +
                               
                               (DistanceRange1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Estimate"]) +
                               (DistanceRange2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Estimate"])),
                
                logitDistUpper = (LM_Final_Effects_Coeff["(Intercept)", "Upper"] +
                                    (midLat.sc * LM_Final_Effects_Coeff["LatitudeScaled", "Upper"]) +
                                    (zeroMed1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Upper"]) +
                                    (zeroMed2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Upper"]) +
                                    
                                    (DistanceRange1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Upper"]) +
                                    (DistanceRange2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Upper"])),
                
                logitDistLower = (LM_Final_Effects_Coeff["(Intercept)", "Lower"] +
                                    (midLat.sc * LM_Final_Effects_Coeff["LatitudeScaled", "Lower"]) +
                                    (zeroMed1 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)1", "Lower"]) +
                                    (zeroMed2 * LM_Final_Effects_Coeff["poly(MedianPropSquared, 2)2", "Lower"]) +
                                    
                                    (DistanceRange1 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)1", "Lower"]) +
                                    (DistanceRange2 * LM_Final_Effects_Coeff["poly(CorrectedDistanceSquared, 2)2", "Lower"]))
  ) %>% 
  dplyr::mutate(across(starts_with("logit"), 
                       function(x) {1/(1+exp(-x))}, 
                       .names = "{col}_line")) %>% 
  dplyr::rename_with(~str_replace(., pattern = "logit", replacement = ""), 
                     ends_with("_line"))

### Latitude ####
ggplot(LM_Final_Effects_Lines, aes(LatitudeRange, Lat_line)) +
  geom_ribbon(aes(ymin = LatLower_line, ymax = LatUpper_line),
              alpha = 0.1) +
  geom_line(col = "#DE8F6E", linewidth = 2, lineend = "round") +
  
  geom_rug(data = GMPD_IUCN_Species, aes(Latitude, Prevalence), sides = "b") +
  
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
ggplot(LM_Final_Effects_Lines, aes(MedianPropRange, Med_line)) +
  geom_ribbon(aes(ymin = MedLower_line, ymax = MedUpper_line),
              alpha = 0.1) +

  geom_rug(data = GMPD_IUCN_Species, aes(MedianPropSquared, Prevalence), sides = "b") +
  
  geom_line(col = "#DE8F6E", linewidth = 2, lineend = "round") +

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

### CorrectedDistanceSquared ####
ggplot(LM_Final_Effects_Lines, aes(DistanceRange, Dist_line)) +
  geom_ribbon(aes(ymin = DistLower_line, ymax = DistUpper_line),
              alpha = 0.1) +
  
  geom_rug(data = GMPD_IUCN_Species, aes(CorrectedDistanceSquared, Prevalence), sides = "b") +
  
  geom_line(col = "#DE8F6E", linewidth = 2, lineend = "round") +
  
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



