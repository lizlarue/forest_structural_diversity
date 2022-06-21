library(lidR)
library(gstat)
library(neonUtilities)
library(dplyr)

data <- read.csv("NEON_sites_dates_for_cover.csv", colClasses = "character")
coverdata <- read.csv("prelim_cover.csv")
all_neon_plot_centroids_data <- read.csv("All_NEON_TOS_Plot_Centroids_V8.csv", colClasses = "character")

plot_data_table = data.frame()
extra <- c()

slicedcoverdata <- unique(coverdata %>% select(siteID, monthyear, sitemonthyear, plotID))
for (i in 1:nrow(slicedcoverdata)){ #nrow(t1)
  slice <- all_neon_plot_centroids_data[all_neon_plot_centroids_data$plotID == slicedcoverdata$plotID[i] & all_neon_plot_centroids_data$subtype == "basePlot",]
  if (nrow(slice) > 0){
    row <- cbind(slicedcoverdata[i,], slice$latitude, slice$longitude, slice$easting, slice$northing)
    plot_data_table <- rbind(plot_data_table, row)
  }
  else{
    extra <- c(extra, as.character(slicedcoverdata[i,]$plotID))
  }
}
plot_data_table <- plot_data_table %>% rename(latitude = "slice$latitude", 
                                              longitude = "slice$longitude",
                                              easting = "slice$easting",
                                              northing = "slice$northing")
extra <- unique(extra)
#"OSBS_051" "OSBS_018" "OSBS_023" "SCBI_033" "SCBI_041" "SCBI_044"
#"SCBI_034" "SCBI_039" "SCBI_042" "SCBI_012" "MLBS_031" "BLAN_002"
#"DELA_010" "DELA_007" "DELA_019" "DELA_024"

write.table(plot_data_table, file = "plot_data_table.csv", sep = ",", row.names = FALSE)
