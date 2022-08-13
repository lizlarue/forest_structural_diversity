library(lidR)
library(gstat)
library(neonUtilities)
library(dplyr)
source("helper.R")

data <- read.csv("NEON_sites_dates_for_cover.csv", colClasses = "character")
coverdata <- read.csv("prelim_cover.csv")
all_neon_plot_centroids_data <- read.csv("All_NEON_TOS_Plot_Centroids_V8.csv", colClasses = "character")

plot_data_table = data.frame()

slicedcoverdata <- unique(coverdata %>% select(siteID, monthyear, sitemonthyear, plotID))
for (i in 1:nrow(slicedcoverdata)){
  slice <- all_neon_plot_centroids_data[all_neon_plot_centroids_data$plotID == slicedcoverdata$plotID[i] & all_neon_plot_centroids_data$subtype == "basePlot",]
  if (nrow(slice) > 0){
    row <- cbind(slicedcoverdata[i,], slice$latitude, slice$longitude, slice$easting, slice$northing)
    plot_data_table <- rbind(plot_data_table, row)
  }
  else{
    slice <- coverdata[coverdata$plotID == slicedcoverdata$plotID[i],][1,]
    slice <- slice %>% rename(latitude = "decimalLatitude", longitude = "decimalLongitude")
    utm_coords <- longlat_to_UTM(slice$longitude, slice$latitude)
    row <- cbind(slicedcoverdata[i,], slice$latitude, slice$longitude, "slice$easting" = utm_coords$easting, "slice$northing" = utm_coords$northing)
    plot_data_table <- rbind(plot_data_table, row)
  }
}
plot_data_table <- plot_data_table %>% rename(latitude = "slice$latitude", 
                                              longitude = "slice$longitude",
                                              easting = "slice$easting",
                                              northing = "slice$northing")

write.table(plot_data_table, file = "plot_data_table.csv", sep = ",", row.names = FALSE)
