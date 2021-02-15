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
#HARV <- readLAS(paste0(wd,"NEON_D01_HARV_DP1_725000_4705000_classified_point_cloud_colorized.laz"))
HARV <- readLAS(paste0(wd,"NEON_D01_HARV_DP1_725000_4705000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below -163 -drop_z_above 660")
#should I filter this to zero?

summary(HARV)
#plot(HARV)


#set center of plot based on extent
x <- ((max(HARV$X) - min(HARV$X))/2)+ min(HARV$X)
y <- ((max(HARV$Y) - min(HARV$Y))/2)+ min(HARV$Y)

data.200m <- lasclipRectangle(HARV, 
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

HARV_structural_diversity <- structural_diversity_metrics(data.40m)



#####################################
#diversity data
devtools::insHARV_github("admahood/neondiversity")

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


coverH <- loadByProduct (dpID = "DP1.10058.001", site = 'HARV', check.size = FALSE)

coverDivH <- coverH[[2]]

unique(coverDivH$divDataType)

cover2H <- coverDivH %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2H$scientificName))

summary(cover2H$nativeStatusCode)

#subset of invasive only
inv <- cover2H %>%
  filter(nativeStatusCode=="I")

#total SR of exotics across all plots
exotic_SR <-length(unique(inv$scientificName))

#mean plot percent cover of exotics
exotic_cover <- inv %>%
  group_by(plotID) %>%
  summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
  summarize(exotic_cov = mean(sumz))


HARV_table <- cbind(HARV_structural_diversity, all_SR, exotic_SR, exotic_cover)

HARV_table <- HARV_table %>%
  mutate(Site.ID = "HARV")

HARV_table <- HARV_table %>%
  select(-easting, -northing)

HARV_table <- HARV_table %>%
  left_join(veg_types)


combo10 <- rbind(combo9, HARV_table)
combo10




####################################
write.table(combo10, file = "prelim_results.csv", sep = ",", row.names = FALSE)

library(ggplot2)
ggplot(combo10, aes(x = mean.max.canopy.ht.aop, y = exotic_SR))+
  geom_point()

ggplot(combo10, aes(x = max.canopy.ht.aop, y = exotic_SR))+
  geom_point()

ggplot(combo10, aes(x = rumple.aop, y = exotic_SR))+
  geom_point()
