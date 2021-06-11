#working with NEON spatial data to spatially align field sampling plots with AOP data


library(sf)
library(rgdal)
library(ggplot2)

#site_locations <- st_read("data/terrestrialSamplingBoundaries.shp")
site_locations <- readOGR(dsn="data/", layer = "terrestrialSamplingBoundaries")
plot_locations <- readOGR(dsn="data/", layer = "NEON_TOS_Plot_Centroids")

ggplot() + 
  geom_polygon(data = site_locations, aes(x = lat, y = long))+
  geom_point(data = plot_locations, aes(x = latitude, y = longitude))

summary(site_locations)
head(site_locations)

plot(site_locations)
plot(plot_locations)




file <- "Burrows_et_al_Nature_traj_ocean_NH1.kmz"

SST_start = readOGR(file,"SST_start") 