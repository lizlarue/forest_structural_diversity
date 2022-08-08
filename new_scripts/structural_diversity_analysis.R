# This code analyzes forest structural diversity in relation to degree of invasion.

#load libraries
#x <- c("tidyverse", "sf","scales", "ggplot2", "doBy", "lme4", "lmerTest", "stargazer","ggpubr")
x <- c("tidyverse", "scales", "ggplot2", "doBy", "lme4", "lmerTest", "stargazer", "ggpubr")
lapply(x, library, character.only = TRUE, verbose = FALSE)
theme_set(theme_classic())

#load data
strdiv <- as.data.frame(read_csv("new_scripts/output/structural_metrics_by_plot.csv"))

head(strdiv)

unique(strdiv$Dominant_NLCD_Classes)

#replace NAs in exotic_cover and exotic_SR with 0s
strdiv$exotic_cover[is.na(strdiv$exotic_cover)] <- 0
strdiv$exotic_richness[is.na(strdiv$exotic_richness)] <- 0


#remove outliers
strdiv$max.canopy.ht.aop <- ifelse(strdiv$max.canopy.ht.aop > 200, NA, strdiv$max.canopy.ht.aop) 
strdiv$mean.max.canopy.ht.aop <- ifelse(strdiv$mean.max.canopy.ht.aop > 200, NA, strdiv$mean.max.canopy.ht.aop)


#reclass NLCD
strdiv <- strdiv %>%
  mutate(new_class = ifelse(Dominant_NLCD_Classes == "Evergreen Forest, Grassland/Herbaceous, Shrub/Scrub" | Dominant_NLCD_Classes == "Evergreen Forest, Shrub/Scrub, Woody Wetlands" | Dominant_NLCD_Classes == "Evergreen Forest, Woody Wetlands" | Dominant_NLCD_Classes == "Dwarf Scrub, Evergreen Forest, Shrub/Scrub" | Dominant_NLCD_Classes == "Evergreen Forest, Shrub/Scrub" | Dominant_NLCD_Classes == "Evergreen Forest, Grassland/Herbaceous" | Dominant_NLCD_Classes == "Emergent Herbaceous Wetlands, Evergreen Forest, Woody Wetlands" | Dominant_NLCD_Classes == "Evergreen Forest", 'evergreen',
                    ifelse(Dominant_NLCD_Classes == "Deciduous Forest, Pasture/Hay" | Dominant_NLCD_Classes == "Deciduous Forest, Grassland/Herbaceous"  | Dominant_NLCD_Classes == "Deciduous Forest, Woody Wetlands" | Dominant_NLCD_Classes == "Deciduous Forest" |  Dominant_NLCD_Classes == "Cultivated Crops, Deciduous Forest" , 'deciduous', 'mixed forest')))

########################################
#linear mixed models with predictors including structral diversity metrics and forest type (fixed effects) and site as random effect
#exotic plant cover
final_mods <- list()

#notes: vert.sd.aop = external heterogeneity, sd.sd.aop = internal heterogeneity
#all predictors
final_mods$cover1 <- lmer(exotic_cover ~ new_class * max.canopy.ht.aop 
                         * mean.max.canopy.ht.aop *deepgap.fraction.aop * vert.sd.aop * sd.sd.aop 
                         + (1|siteID), data=strdiv, REML=FALSE)

#subset1
final_mods$cover2 <- lmer(exotic_cover ~ new_class * max.canopy.ht.aop 
                          * mean.max.canopy.ht.aop + (1|siteID), data=strdiv, REML=FALSE)

#subset2
final_mods$cover3 <- lmer(exotic_cover ~ new_class * deepgap.fraction.aop 
                          + (1|siteID), data=strdiv, REML=FALSE)


#subset3
final_mods$cover4 <- lmer(exotic_cover ~ new_class * vert.sd.aop * sd.sd.aop 
                          + (1|siteID), data=strdiv, REML=FALSE)


summary(final_mods$cover1)
summary(final_mods$cover2)
summary(final_mods$cover3)
summary(final_mods$cover4)


write.table(summary(final_mods$cover1))
#In addition: Warning message:
#  Some predictor variables are on very different scales: consider rescaling 





#exotic plant species richness 




####################################

#regression with exotic cover or SR (response) as a function of structural diversity (predictors)
lmcover1 = lm(exotic_cover~ max.canopy.ht.aop, data = strdiv)
summary (lmcover1)
#significantly related

lmcover2 = lm(exotic_cover~ mean.max.canopy.ht.aop, data = strdiv)
summary (lmcover2)
#significantly related

lmcover3 = lm(exotic_cover~ deepgap.fraction.aop, data = strdiv)
summary (lmcover3)
#not sig related

lmcover4 = lm(exotic_cover~ vert.sd.aop, data = strdiv)
summary (lmcover4)
#not sig related

lmcover5 = lm(exotic_cover~ sd.sd.aop, data = strdiv)
summary (lmcover5)
#not sig related


#regression with exotic SR 
lmcover6 = lm(exotic_richness~ max.canopy.ht.aop, data = strdiv)
summary (lmcover6)
#not sig related

lmcover7 = lm(exotic_richness~ mean.max.canopy.ht.aop, data = strdiv)
summary (lmcover7)
#sig related

lmcover8 = lm(exotic_richness~ deepgap.fraction.aop, data = strdiv)
summary (lmcover8)
#not sig related

lmcover9 = lm(exotic_richness~ vert.sd.aop, data = strdiv)
summary (lmcover9)
#highly significant

lmcover10 = lm(exotic_richness~ sd.sd.aop, data = strdiv)
summary (lmcover10)
#highly significant



#####
#plotting
ggplot(strdiv, aes(x=max.canopy.ht.aop, y=exotic_cover)) + 
  geom_point()

ggplot(strdiv, aes(x=mean.max.canopy.ht.aop, y=exotic_cover)) + 
  geom_point()

ggplot(strdiv, aes(x=deepgap.fraction.aop, y=exotic_cover)) + 
  geom_point()

ggplot(strdiv, aes(x=vert.sd.aop, y=exotic_cover)) + 
  geom_point()

ggplot(strdiv, aes(x=sd.sd.aop, y=exotic_cover)) + 
  geom_point()








ggplot(strdiv, aes(x=max.canopy.ht.aop, y=exotic_richness)) + 
  geom_point()

ggplot(strdiv, aes(x=mean.max.canopy.ht.aop, y=exotic_richness)) + 
  geom_point()

ggplot(strdiv, aes(x=deepgap.fraction.aop, y=exotic_richness)) + 
  geom_point()

ggplot(strdiv, aes(x=vert.sd.aop, y=exotic_richness)) + 
  geom_point()

ggplot(strdiv, aes(x=sd.sd.aop, y=exotic_richness)) + 
  geom_point()

