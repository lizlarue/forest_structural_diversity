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
#ONAQ <- readLAS(paste0(wd,"NEON_D15_ONAQ_DP1_367000_4454000_classified_point_cloud_colorized.laz"))
ONAQ <- readLAS(paste0(wd,"NEON_D15_ONAQ_DP1_367000_4454000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 1960 -drop_z_above 2650")

summary(ONAQ)
#plot(ONAQ)

#set center of plot based on extent
x <- ((max(ONAQ$X) - min(ONAQ$X))/2)+ min(ONAQ$X)
y <- ((max(ONAQ$Y) - min(ONAQ$Y))/2)+ min(ONAQ$Y)

#x <- 367500 #easting 
#y <- 4454500 #northing

data.200m <- lasclipRectangle(ONAQ, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

data.200m <- lasnormalize(data.200m, dtm)

data.40m <- lasclipRectangle(data.200m, 
                             xleft = (x - 20), ybottom = (y - 20),
                             xright = (x + 20), ytop = (y + 20))
data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
plot(data.40m)
#this does not look right!!!

############# Structural diversity metrics  ######

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

ONAQ_structural_diversity <- structural_diversity_metrics(data.40m)





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

#
coverO <- loadByProduct (dpID = "DP1.10058.001", site = 'ONAQ', check.size = FALSE)

coverDivO <- coverO[[2]]

unique(coverDivO$divDataType)

cover2O <- coverDivO %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2O$scientificName))

summary(cover2O$nativeStatusCode)

#subset of invasive only
inv <- cover2O %>%
  filter(nativeStatusCode=="I")

#total SR of exotics across all plots
exotic_SR <-length(unique(inv$scientificName))

#mean plot percent cover of exotics
exotic_cover <- inv %>%
  group_by(plotID) %>%
  summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
  summarize(exotic_cov = mean(sumz))


#exotic_cover <- sum(inv$percentCover, na.rm = TRUE)


ONAQ_table <- cbind(ONAQ_structural_diversity, all_SR, exotic_SR, exotic_cover)

ONAQ_table <- ONAQ_table %>%
  mutate(Site.ID = "ONAQ")


ONAQ_table <- ONAQ_table %>%
  select(-easting, -northing)

ONAQ_table <- ONAQ_table %>%
  left_join(veg_types)


combo4 <- rbind(combo3, ONAQ_table)
combo4


######################
#soil chem for ONAQ
ONAQ_soil_chem_0 <- read.csv(file = 'NEON.D15.ONAQ.DP1.00096.001.mgp_perbiogeosample.2014-06.basic.20201006T215549Z.csv')


ONAQ_soil_chem <- read.csv(file = 'NEON.D15.ONAQ.DP1.00096.001.mgp_perbiogeosample.2014-06.basic.20201006T215549Z.csv') %>%
  select(siteID, carbonTot, horizonName) %>%
  rename(Site.ID = siteID) %>%
  mutate(horizon = ifelse(horizonName == "A1" | horizonName == "A2", "A", ifelse(horizonName == "Bw" | horizonName == "Bk"| horizonName == "2Bkm", "B", "C")))

