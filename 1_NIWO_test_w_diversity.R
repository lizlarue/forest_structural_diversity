#structural diversity tutorial
#https://www.neonscience.org/structural-diversity-discrete-return

library(lidR)
library(gstat)
library(neondiversity)

wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

#NIWO <- readLAS(paste0(wd,"NEON_D13_NIWO_DP1_454000_4425000_classified_point_cloud_colorized.laz"))
NIWO <- readLAS(paste0(wd,"NEON_D13_NIWO_DP1_454000_4425000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 2460 -drop_z_above 2947")
summary(NIWO)

#Let's correct for elevation and measure structural diversity for NIWO
x <- ((max(NIWO$X) - min(NIWO$X))/2)+ min(NIWO$X)
y <- ((max(NIWO$Y) - min(NIWO$Y))/2)+ min(NIWO$Y)

#x <- 454500 
#y <- 4425608

data.200m <- lasclipRectangle(NIWO, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

data.200m <- lasnormalize(data.200m, dtm)

data.40m <- lasclipRectangle(data.200m, 
                             xleft = (x - 20), ybottom = (y - 20),
                             xright = (x + 20), ytop = (y + 20))
data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
plot(data.40m)


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

NIWO_structural_diversity <- structural_diversity_metrics(data.40m)





#####################################
#diversity data
devtools::install_github("admahood/neondiversity")

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

coverN <- loadByProduct (dpID = "DP1.10058.001", site = 'NIWO')

coverDivN <- coverN[[2]]

unique(coverDivN$divDataType)

cover2N <- coverDivN %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2N$scientificName))

summary(cover2N$nativeStatusCode)
#no 'I's...what are NIs?

#subset of invasive only
inv <- cover2N %>%
  filter(nativeStatusCode=="I")

exotic_SR <-length(unique(inv$scientificName))


NIWO_table <- cbind(NIWO_structural_diversity, all_SR, exotic_SR)

NIWO_table <- NIWO_table %>%
  mutate(Site.ID = "NIWO")

NIWO_table <- NIWO_table %>%
  select(-easting, -northing)

veg_types <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/field-sites.csv') %>%
  select(Site.ID, Dominant.NLCD.Classes)

NIWO_table <- NIWO_table %>%
  left_join(veg_types)

NIWO_table



#####
#for all sites
#all_sites_table <- all_sites_table %>%
  #select(-easting, -northing)
  

#carSpeeds <- read.csv(file = 'data/car-speeds.csv')
#veg_types <- read.csv(file = 'field-sites.csv') %>%
  #select(Site.ID, Dominant.NLCD.Classes)

#add veg class to table 
#all_sites_table <- all_sites_table %>%
  #left_join(veg_types)


######################
#soil chem for NIWO
NIWO_soil_chem_0 <- read.csv(file = 'NEON.D13.NIWO.DP1.00097.001.mgc_perbiogeosample.2014-08.basic.20191014T212752Z.csv')
  

NIWO_soil_chem <- read.csv(file = 'NEON.D13.NIWO.DP1.00097.001.mgc_perbiogeosample.2014-08.basic.20191014T212752Z.csv') %>%
  select(siteID, carbonTot, horizonName) %>%
  rename(Site.ID = siteID) %>%
  mutate(horizon = ifelse(horizonName == "A1" | horizonName == "A2", "A", ifelse(horizonName == "Bw1" | horizonName == "Bw2", "B", "C")))


