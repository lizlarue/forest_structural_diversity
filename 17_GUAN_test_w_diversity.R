#structural diversity tutorial
#https://www.neonscience.org/structural-diversity-discrete-return

library(lidR)
library(gstat)
library(neondiversity)

############### Set working directory ######
#set the working of the downloaded data
wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

############ No LiDAR data for GUAN###########
#GUAN <- readLAS(paste0(wd,"NEON_D16_GUAN_DP1_550000_5064000_classified_point_cloud_colorized.laz"))
#GUAN <- readLAS(paste0(wd,"NEON_D16_GUAN_DP1_550000_5064000_classified_point_cloud_colorized.laz"),
                #filter = "-drop_z_below 202 -drop_z_above 774")

#summary(GUAN)
#plot(GUAN)


#set center of plot based on extent
#x <- ((max(GUAN$X) - min(GUAN$X))/2)+ min(GUAN$X)
#y <- ((max(GUAN$Y) - min(GUAN$Y))/2)+ min(GUAN$Y)

#data.200m <- lasclipRectangle(GUAN, 
                              #xleft = (x - 100), ybottom = (y - 100),
                              #xright = (x + 100), ytop = (y + 100))

#dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

#data.200m <- lasnormalize(data.200m, dtm)

#data.40m <- lasclipRectangle(data.200m, 
                             #xleft = (x - 20), ybottom = (y - 20),
                             #xright = (x + 20), ytop = (y + 20))
#data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
#plot(data.40m)
#looks good

#Zip up all the code we previously used and write function to 
#run all 13 metrics in a single function. 
#structural_diversity_metrics <- function(data.40m) {
  #chm <- grid_canopy(data.40m, res = 1, dsmtin()) 
  #mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  #max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  #rumple <- rumple_index(chm) 
  #top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  #cells <- length(chm@data@values) 
  #chm.0 <- chm
  #chm.0[is.na(chm.0)] <- 0 
  #zeros <- which(chm.0@data@values == 0) 
  #deepgaps <- length(zeros) 
  #deepgap.fraction <- deepgaps/cells 
  #cover.fraction <- 1 - deepgap.fraction 
  #vert.sd <- cloud_metrics(data.40m, sd(Z, na.rm = TRUE)) 
  #sd.1m2 <- grid_metrics(data.40m, sd(Z), 1) 
  #sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  #Zs <- data.40m@data$Z
  #Zs <- Zs[!is.na(Zs)]
  #entro <- entropy(Zs, by = 1) 
  #gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
  #GFP.AOP <- mean(gap_frac$gf) 
  #LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  #VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  #VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  #out.plot <- data.frame(
    #matrix(c(x, y, mean.max.canopy.ht,max.canopy.ht, 
             #rumple,deepgaps, deepgap.fraction, 
             #cover.fraction, top.rugosity, vert.sd, 
             #sd.sd, entro, GFP.AOP, VAI.AOP,VCI.AOP),
          #ncol = 15)) 
  #colnames(out.plot) <- 
    #c("easting", "northing", "mean.max.canopy.ht.aop",
     # "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
      #"deepgap.fraction.aop", "cover.fraction.aop",
      #"top.rugosity.aop","vert.sd.aop","sd.sd.aop", 
      #"entropy.aop", "GFP.AOP.aop",
      #"VAI.AOP.aop", "VCI.AOP.aop") 
  #print(out.plot)
#}

#GUAN_structural_diversity <- structural_diversity_metrics(data.40m)


#put in NAs since there are no lidar data
head(combo15)
GUAN_structural_diversity <- data.frame("easting" = NA, "northing" = NA, "mean.max.canopy.ht.aop" = NA, "max.canopy.ht.aop" = NA, "rumple.aop" = NA, "deepgaps.aop" = NA, "deepgap.fraction.aop" = NA, "cover.fraction.aop" = NA, "top.rugosity.aop" = NA, "vert.sd.aop" = NA, "sd.sd.aop" = NA, "entropy.aop" = NA, "GFP.AOP.aop" = NA, "VAI.AOP.aop" = NA, "VCI.AOP.aop" = NA)
#GUAN_structural_diversity$easting <-NA
#GUAN_structural_diversity$northing <- NA
#GUAN_structural_diversity$mean.max.canopy.ht.aop <- NA
#GUAN_structural_diversity$max.canopy.ht.aop <- NA
#GUAN_structural_diversity$rumple.aop <- NA
#GUAN_structural_diversity$deepgaps.aop <- NA
#GUAN_structural_diversity$deepgap.fraction.aop <- NA
#GUAN_structural_diversity$cover.fraction.aop <- NA
#GUAN_structural_diversity$top.rugosity.aop <- NA
#GUAN_structural_diversity$vert.sd.aop <- NA
#GUAN_structural_diversity$sd.sd.aop <- NA
#GUAN_structural_diversity$entropy.aop <- NA
#GUAN_structural_diversity$GFP.AOP.aop <- NA
#GUAN_structural_diversity$VAI.AOP.aop <- NA
#GUAN_structural_diversity$VCI.AOP.aop <- NA




#####################################
#diversity data
devtools::insGUAN_github("admahood/neondiversity")

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
cover <- loadByProduct (dpID = "DP1.10058.001", site = 'GUAN', check.size = FALSE)

coverDiv <- cover[[2]]

unique(coverDiv$divDataType)

cover2 <- coverDiv %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2$scientificName))

summary(cover2$nativeStatusCode)

#subset of invasive only
inv <- cover2 %>%
  filter(nativeStatusCode=="I")

#total SR of exotics across all plots
exotic_SR <-length(unique(inv$scientificName))

#mean plot percent cover of exotics
exotic_cover <- inv %>%
  group_by(plotID) %>%
  summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
  summarize(exotic_cov = mean(sumz))



GUAN_table <- cbind(GUAN_structural_diversity, all_SR, exotic_SR, exotic_cover)

GUAN_table <- GUAN_table %>%
  mutate(Site.ID = "GUAN")

GUAN_table <- GUAN_table %>%
  select(-easting, -northing)

GUAN_table <- GUAN_table %>%
  left_join(veg_types)



#############################################
#calculate spectral reflectance as CV
#as defined here: https://www.mdpi.com/2072-4292/8/3/214/htm
f <- paste0(wd,"NEON_D04_GUAN_DP3_730000_1990000_reflectance.h5")


###
#for each of the 426 bands, I need to calculate the mean reflectance and the SD reflectance across all pixels 

myNoDataValue <- as.numeric(reflInfo$Data_Ignore_Value)

dat <- data.frame()

for (i in 1:426){
  #extract one band
  b <- h5read(f,"/GUAN/Reflectance/Reflectance_Data",index=list(i,1:nCols,1:nRows)) 
  
  # set all values equal to -9999 to NA
  b[b == myNoDataValue] <- NA
  
  #calculate mean and sd
  meanref <- mean(b, na.rm = TRUE)
  SDref <- sd(b, na.rm = TRUE)
  
  rowz <- cbind(i, meanref, SDref)
  
  dat <- rbind(dat, rowz)
}


dat$calc <- dat$SDref/dat$meanref

CV <- sum(dat$calc)/426


GUAN_table$specCV <- CV


combo16 <- rbind(combo15, GUAN_table)
combo16







write.table(combo16, file = "prelim_results.csv", sep = ",", row.names = FALSE)

#remove SRER (shrub) and CUPE (no cover data)
combo16sub <-subset(combo16, Site.ID!="CUPE" & Site.ID!="SRER")
#combo16sub <-subset(combo15, Site.ID!="CUPE" & Site.ID!="SRER" & Site.ID!="SCBI")

combo16sub$native_SR <-combo16sub$all_SR - combo16sub$exotic_SR


library(ggplot2)


#linear regression with structural data and exotic SR
#external heterogeneity
ggplot(combo16sub, aes(x = top.rugosity.aop, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")
#expect to see negative relationship here

fit <- lm(exotic_SR ~ top.rugosity.aop, data = combo16sub)
summary(fit)
#r2 = 0.115, p=0.23

#internal heterogeneity
ggplot(combo16sub, aes(x = sd.sd.aop, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_SR ~ sd.sd.aop, data = combo16sub)
summary(fit)
#r2 = 0.011, p=0.71

#mean canopy height
ggplot(combo16sub, aes(x = mean.max.canopy.ht.aop, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_SR ~ mean.max.canopy.ht.aop, data = combo16sub)
summary(fit)
#r2 = 0.023, p=0.88

#gap fraction
ggplot(combo16sub, aes(x = deepgap.fraction.aop, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_SR ~ deepgap.fraction.aop, data = combo16sub)
summary(fit)
#r2 = 2.396e-06, p=0.99

#max canopy height
ggplot(combo16sub, aes(x = max.canopy.ht.aop, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

fit <- lm(exotic_SR ~ max.canopy.ht.aop, data = combo16sub)
summary(fit)
#r2 = 0.066, p=0.37

#ratio of outer canopy surface area to ground surface area 
ggplot(combo16sub, aes(x = rumple.aop, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

fit <- lm(exotic_SR ~ rumple.aop, data = combo16sub)
summary(fit)
#r2 = 0.00034, p=0.94

#####################################################
#linear regression with structural diversity and cover data

#external heterogeneity
ggplot(combo16sub, aes(x = top.rugosity.aop, y = log(exotic_cov +1), label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")
#expect to see negative relationship here

fit <- lm(log(exotic_cov +1) ~ top.rugosity.aop, data = combo16sub)
summary(fit)
#r2 = 0.028, p=0.567

#internal heterogeneity
ggplot(combo16sub, aes(x = sd.sd.aop, y = log(exotic_cov +1), label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(log(exotic_cov +1) ~ sd.sd.aop, data = combo16sub)
summary(fit)
#r2 = 0.005, p=0.81

#mean canopy height
ggplot(combo16sub, aes(x = mean.max.canopy.ht.aop, y = log(exotic_cov +1), label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(log(exotic_cov +1) ~ mean.max.canopy.ht.aop, data = combo16sub)
summary(fit)
#r2 = 0.0003, p=0.95

#gap fraction
ggplot(combo16sub, aes(x = deepgap.fraction.aop, y = log(exotic_cov +1), label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(log(exotic_cov +1) ~ deepgap.fraction.aop, data = combo16sub)
summary(fit)
#r2 = 0.0049, p=0.81

#max canopy height
ggplot(combo16sub, aes(x = max.canopy.ht.aop, y = log(exotic_cov +1), label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

fit <- lm(log(exotic_cov +1) ~ max.canopy.ht.aop, data = combo16sub)
summary(fit)
#r2 = 0.238, p=0.0764

#ratio of outer canopy surface area to ground surface area 
ggplot(combo16sub, aes(x = rumple.aop, y = log(exotic_cov +1), label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

fit <- lm(log(exotic_cov +1) ~ rumple.aop, data = combo16sub)
summary(fit)
#r2 = 0.01186, p=0.711







#####################################################
#linear regression with spectral data and exotic SR and cover
ggplot(combo16sub, aes(x = specCV, y = exotic_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")
#expect to see negative relationship here

fit <- lm(exotic_SR ~ specCV, data = combo16sub)
summary(fit)
#r2 = 0.005, p=0.7974


ggplot(combo16sub, aes(x = specCV, y = exotic_cov, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
#expect to see negative relationship here

fit <- lm(exotic_cov ~ specCV, data = combo16sub)
summary(fit)
#r2 = 0.09527, p=0.263



#####################################################
#multiple regression with structural diversity and specCV

#external heterogeneity
fit <- lm(exotic_SR ~ top.rugosity.aop + specCV, data = combo16sub)
summary(fit)
#r2 = 0.1156, p=0.5088

#internal heterogeneity
fit <- lm(exotic_SR ~ sd.sd.aop + specCV, data = combo16sub)
summary(fit)
#r2 = 0.0199, p=0.8952

#mean canopy height
fit <- lm(exotic_SR ~ mean.max.canopy.ht.aop + specCV, data = combo16sub)
summary(fit)
#r2 = 0.0157, p=0.9167

#gap fraction
fit <- lm(exotic_SR ~ deepgap.fraction.aop + specCV, data = combo16sub)
summary(fit)
#r2 = 0.02167, p=0.8865

#max canopy height
fit <- lm(exotic_SR ~ max.canopy.ht.aop + specCV, data = combo16sub)
summary(fit)
#r2 = 0.075, p=0.6476

#ratio of outer canopy surface area to ground surface area 
fit <- lm(exotic_SR ~ rumple.aop + specCV, data = combo16sub)
summary(fit)
#r2 = 0.02473, p=0.8713



#####################################################
#multiple regression with structural diversity, native_SR and specCV

#external heterogeneity
fit <- lm(exotic_SR ~ top.rugosity.aop + native_SR + specCV, data = combo16sub)
summary(fit)
#r2 = 0.12, p=0.7076

#internal heterogeneity
fit <- lm(exotic_SR ~ sd.sd.aop + native_SR + specCV, data = combo16sub)
summary(fit)
#r2 = 0.0244, p=0.9675

#mean canopy height
fit <- lm(exotic_SR ~ mean.max.canopy.ht.aop + native_SR + specCV, data = combo16sub)
summary(fit)
#r2 = 0.0199, p=0.9758

#gap fraction
fit <- lm(exotic_SR ~ deepgap.fraction.aop + native_SR + specCV, data = combo16sub)
summary(fit)
#r2 = 0.028, p=0.9595

#max canopy height
fit <- lm(exotic_SR ~ max.canopy.ht.aop + native_SR + specCV, data = combo16sub)
summary(fit)
#r2 = 0.1163, p=0.73

#ratio of outer canopy surface area to ground surface area 
fit <- lm(exotic_SR ~ rumple.aop + native_SR + specCV, data = combo16sub)
summary(fit)
#r2 = 0.02687, p=0.9628



################################################
#spectral diversity vs. structural diversity

#external heterogeneity
fit <- lm(specCV ~ top.rugosity.aop, data = combo16sub)
summary(fit)
#r2 = 0.1655, p = 0.1489

#internal heterogeneity
fit <- lm(specCV ~ sd.sd.aop, data = combo16sub)
summary(fit)
#r2= 0.122, p = 0.2205

#mean canopy height
fit <- lm(specCV ~ mean.max.canopy.ht.aop, data = combo16sub)
summary(fit)
#r2=0.2863, p=0.04863

#mean canopy height and exotic SR
fit <- lm(specCV ~ mean.max.canopy.ht.aop + exotic_SR, data = combo16sub)
summary(fit)
#r2=0.2962, p=0.1449

#mean canopy height and exotic cov
fit <- lm(specCV ~ mean.max.canopy.ht.aop + exotic_cov, data = combo16sub)
summary(fit)
#r2=0.4067, p=0.05664

#exotic cov
fit <- lm(specCV ~ exotic_cov, data = combo16sub)
summary(fit)
#r2=0.09527, p=0.263

#gap fraction
fit <- lm(specCV ~ deepgap.fraction.aop, data = combo16sub)
summary(fit)
#r2=0.2967, p=0.0440

#gap fraction and exotic SR
fit <- lm(specCV ~ deepgap.fraction.aop + exotic_SR, data = combo16sub)
summary(fit)
#r2=0.3119, p=0.128

#gap fraction and exotic cov
fit <- lm(specCV ~ deepgap.fraction.aop + exotic_cov, data = combo16sub)
summary(fit)
#r2=0.4368, p=0.04251

#spec div and exotic
fit <- lm(specCV ~ exotic_cov, data = combo16sub)
summary(fit)
#r2=0.09527, p=0.263

#max canopy height
fit <- lm(specCV ~ max.canopy.ht.aop, data = combo16sub)
summary(fit)
#r2=0.011, p=0.72

#ratio of outer canopy surface area to ground surface area 
fit <- lm(specCV ~ rumple.aop, data = combo16sub)
summary(fit)
#r2=0.28, p=0.05157

#choose most dominant cover type
combo16sub$dominantforest <- gsub( " .*$", "", combo16sub$Dominant.NLCD.Classes)

#mean canopy height
ggplot(combo16sub, aes(x = mean.max.canopy.ht.aop, y = specCV, label = Site.ID, color = dominantforest))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")

#gap fractionn
ggplot(combo16sub, aes(x = deepgap.fraction.aop, y = specCV, label = Site.ID, color = dominantforest))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")

#ratio of outer canopy surface area to ground surface area
ggplot(combo16sub, aes(x = rumple.aop, y = specCV, label = Site.ID, color = dominantforest))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")

################################################
#linear mixed models
#getting this error: Error: number of levels of each grouping factor must be < number of observations
library(lme4)
library(lmerTest)

final_mods <- list()
final_mods$exotic_cov <- lmer(exotic_cov ~ top.rugosity.aop * native_SR * specCV *Dominant.NLCD.Classes + (1|Site.ID), 
                       data=combo16sub, REML = TRUE)




################################################
###
###these are less useful
ggplot(combo16sub, aes(x = mean.max.canopy.ht.aop, y = exotic_SR/all_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo16sub, aes(x = max.canopy.ht.aop, y = exotic_SR/all_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo16sub, aes(x = rumple.aop, y = exotic_SR/all_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo16sub, aes(x = deepgap.fraction.aop, y = exotic_SR/all_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)

ggplot(combo16sub, aes(x = top.rugosity.aop, y = exotic_SR/all_SR, label = Site.ID))+
  geom_point()+
  geom_text(aes(label=Site.ID),hjust=0, vjust=0)
