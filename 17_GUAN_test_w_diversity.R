#structural diversity tutorial
#https://www.neonscience.org/structural-diversity-discrete-return

library(lidR)
library(gstat)
library(neondiversity)

############### Set working directory ######
#set the working of the downloaded data
wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

############ Read in LiDAR data ###########
#ABBY <- readLAS(paste0(wd,"NEON_D16_ABBY_DP1_550000_5064000_classified_point_cloud_colorized.laz"))
ABBY <- readLAS(paste0(wd,"NEON_D16_ABBY_DP1_550000_5064000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 202 -drop_z_above 774")

summary(ABBY)
#plot(ABBY)


#set center of plot based on extent
x <- ((max(ABBY$X) - min(ABBY$X))/2)+ min(ABBY$X)
y <- ((max(ABBY$Y) - min(ABBY$Y))/2)+ min(ABBY$Y)

data.200m <- lasclipRectangle(ABBY, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

data.200m <- lasnormalize(data.200m, dtm)

data.40m <- lasclipRectangle(data.200m, 
                             xleft = (x - 20), ybottom = (y - 20),
                             xright = (x + 20), ytop = (y + 20))
data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
plot(data.40m)
#looks good

#Zip up all the code we previously used and write function to 
#run all 13 metrics in a single function. 
structural_diversity_metrics <- function(data.40m) {
  chm <- grid_canopy(data.40m, res = 1, dsmtin()) 
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  rumple <- rumple_index(chm) 
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) 
  deepgap.fraction <- deepgaps/cells 
  cover.fraction <- 1 - deepgap.fraction 
  vert.sd <- cloud_metrics(data.40m, sd(Z, na.rm = TRUE)) 
  sd.1m2 <- grid_metrics(data.40m, sd(Z), 1) 
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  Zs <- data.40m@data$Z
  Zs <- Zs[!is.na(Zs)]
  entro <- entropy(Zs, by = 1) 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
  GFP.AOP <- mean(gap_frac$gf) 
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  out.plot <- data.frame(
    matrix(c(x, y, mean.max.canopy.ht,max.canopy.ht, 
             rumple,deepgaps, deepgap.fraction, 
             cover.fraction, top.rugosity, vert.sd, 
             sd.sd, entro, GFP.AOP, VAI.AOP,VCI.AOP),
           ncol = 15)) 
  colnames(out.plot) <- 
    c("easting", "northing", "mean.max.canopy.ht.aop",
      "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
      "deepgap.fraction.aop", "cover.fraction.aop",
      "top.rugosity.aop","vert.sd.aop","sd.sd.aop", 
      "entropy.aop", "GFP.AOP.aop",
      "VAI.AOP.aop", "VCI.AOP.aop") 
  print(out.plot)
}

ABBY_structural_diversity <- structural_diversity_metrics(data.40m)



#####################################
#diversity data
devtools::insABBY_github("admahood/neondiversity")

# load packages
library(neonUtilities)
library(geoNEON)
library(dplyr, quietly=T)
library(downloader)
library(ggplot2)
library(tidyr)
library(doBy)
library(sf)
library(sp)
library(devtools)
library(neondiversity)

#
cover <- loadByProduct (dpID = "DP1.10058.001", site = 'ABBY')

coverDiv <- cover[[2]]

unique(coverDiv$divDataType)

cover2 <- coverDiv %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2$scientificName))

summary(cover2$nativeStatusCode)

#subset of invasive only
inv <- cover2 %>%
  filter(nativeStatusCode=="I")

exotic_SR <-length(unique(inv$scientificName))


ABBY_table <- cbind(ABBY_structural_diversity, all_SR, exotic_SR)

ABBY_table <- ABBY_table %>%
  mutate(Site.ID = "ABBY")

ABBY_table <- ABBY_table %>%
  select(-easting, -northing)

ABBY_table <- ABBY_table %>%
  left_join(veg_types)


combo15 <- rbind(combo14, ABBY_table)
combo15

write.table(combo15, file = "prelim_results.csv", sep = ",", row.names = FALSE)

#combo15 <- combo15 %>%
#mutate(cover = if_else(Dominant.NLCD.Classes == "Shrub/Scrub", "shrub", "forest"))

library(ggplot2)
#library(ggpmisc)
#my.formula <- y ~ x

#external heterogeneity
ggplot(combo15, aes(x = top.rugosity.aop, y = exotic_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")
#expect to see negative relationship here

fit <- lm(exotic_SR ~ top.rugosity.aop, data = combo15)
summary(fit)
#r2 = 0.094, p=0.285

#internal heterogeneity
ggplot(combo15, aes(x = sd.sd.aop, y = exotic_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_SR ~ sd.sd.aop, data = combo15)
summary(fit)
#r2 = 0.013, p=0.69

#mean canopy height
ggplot(combo15, aes(x = mean.max.canopy.ht.aop, y = exotic_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_SR ~ mean.max.canopy.ht.aop, data = combo15)
summary(fit)
#r2 = 0.0006, p=0.93

#gap fraction
ggplot(combo15, aes(x = deepgap.fraction.aop, y = exotic_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_SR ~ deepgap.fraction.aop, data = combo15)
summary(fit)
#r2 = 0.013, p=0.908

#max canopy height
ggplot(combo15, aes(x = max.canopy.ht.aop, y = exotic_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

fit <- lm(exotic_SR ~ max.canopy.ht.aop, data = combo15)
summary(fit)
#r2 = 0.056, p=0.417

#ratio of outer canopy surface area to ground surface area 
ggplot(combo15, aes(x = rumple.aop, y = exotic_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

fit <- lm(exotic_SR ~ rumple.aop, data = combo15)
summary(fit)
#r2 = 0.00058, p=0.93



###
###these are not that useful
ggplot(combo15, aes(x = mean.max.canopy.ht.aop, y = exotic_SR/all_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo15, aes(x = max.canopy.ht.aop, y = exotic_SR/all_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo15, aes(x = rumple.aop, y = exotic_SR/all_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo15, aes(x = deepgap.fraction.aop, y = exotic_SR/all_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo15, aes(x = top.rugosity.aop, y = exotic_SR/all_SR, color = cover, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)