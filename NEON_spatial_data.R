#working with NEON spatial data to spatially align field sampling plots with AOP data

library(sf)
library(rgdal)
library(ggplot2)

wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

#site_locations <- st_read("data/terrestrialSamplingBoundaries.shp")
site_locations <- readOGR(dsn="data/", layer = "terrestrialSamplingBoundaries")
plot_locations <- readOGR(dsn="data/", layer = "NEON_TOS_Plot_Centroids")
states <- readOGR(dsn="data/", layer = "s_11au16")

tot_table_plots <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/cover_by_plot.csv')




#initial plot
plot(site_locations)
#this looks like it works
plot(plot_locations)
#this does not look like it works
plot(states)
#takes a long time to load



ggplot() + 
  geom_polygon(data = site_locations, aes(x = lat, y = long))
#Regions defined for each Polygons

#does not have the site polygon, obviously, but avoids the Regions problem from geom_polygon
ggplot() + 
  geom_point(data = site_locations, aes(x = lat, y = long))

ggplot() + 
  geom_point(data = plot_locations, aes(x = lat, y = long))
#Error: `data` must be a data frame, or other object coercible by `fortify()`, not an S4 object with class SpatialPointsDataFrame

summary(site_locations)
head(site_locations)


#bring in shapefile of US states








###
#bring in kmz of AOP flight boundaries

file <- "data/Burrows_et_al_Nature_traj_ocean_NH1.kmz"
SST_start = readOGR(dsn = file, layer = "SST_start") 




###plot plot cover data
ggplot() + 
  geom_point(data = tot_table_plots, aes(x = decimalLatitude, y = decimalLongitude, color = Dominant.NLCD.Classes))

