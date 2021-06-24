#mosaic lidar images

# load packages
library(lidR)
library(gstat)
library(neondiversity)
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
library(gdata)

wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

#files <- list.files(path="DP1.30003.001/2019/FullSite/D17/2019_SOAP_4/L1/DiscreteLidar/ClassifiedPointCloud", pattern="*.laz", full.names=TRUE, recursive=FALSE)


ctg <- readLAScatalog("DP1.30003.001/2019/FullSite/D17/2019_SOAP_4/L1/DiscreteLidar/ClassifiedPointCloud")

plot(ctg)
head(ctg)

tot_table_plots_en <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/tot_table_plots_en.csv')

ggplot()+
  geom_polygon(data = ctg, x = decimalLatitude, y = decimalLongitude)



str_table = data.frame()

#files <- list.files(path="DP1.30003.001/2019/FullSite/D17/2019_SOAP_4/L1/DiscreteLidar/ClassifiedPointCloud", pattern="*.laz", full.names=TRUE, recursive=FALSE)
files <- list.files(path="DP1.30003.001", pattern="*.laz", full.names=TRUE, recursive=TRUE)


files 


#for (q in files)  {
  #i <- readLAS(q) # load each file
  # apply function to remove noise

opt_output_files(ctg) <- "DP1.30003.001/2019_SOAP_classified"
  i <- lidR::classify_noise(ctg, sor(15,3))
  i <- filter_poi(i, Classification != LASNOISE)
  
  
  
  #find extent of tile
  minx <- min(i$X) 
  maxx <- max(i$X)
  miny <- min(i$Y)
  maxy <- max(i$Y)
  
  #select only the plot centroids that are found in this tile
  matches<-filter(plots, plots$easting <= maxx & plots$easting >= minx & plots$northing <= maxy & plots$northing >= miny) 
  
  #loop through the matches
  foreach(x = matches$easting, y = matches$northing) %do% {
    
    #subset to 40 m x 40 m (1600m2) subtile
    data.40m <- clip_rectangle(i, 
                               xleft = (x - 20), ybottom = (y - 20),
                               xright = (x + 20), ytop = (y + 20))
    
    
    #creates rasterized digital terrain model
    dtm <- grid_terrain(data.40m, 1, kriging(k = 10L))
    
    #normalize the data based on the digital terrain model (correct for elevation)
    data.40m <- normalize_height(data.40m, dtm)
    
    #replaces really short veg with zeros
    data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
    
    
    
    #Calculate 13 structural metrics in a single function 
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
    
    #create table that contains 1 row for each plot centroid
    new_structural_diversity <- structural_diversity_metrics(data.40m)
    str_table <- rbind(str_table, new_structural_diversity)
  }  



str_table



#join cover table with structure table by easting and northing
tot_table_plots_en_str <- tot_table_plots_en %>%
  left_join(str_table)



