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
