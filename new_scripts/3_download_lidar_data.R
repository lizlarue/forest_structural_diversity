library(lidR)
library(gstat)
library(neonUtilities)
library(dplyr)

lastools_bin_path <- "C:/Users/Lenovo/LAStools/LAStools/bin/"
root <- paste0(getwd(), "/")

# if folder does not exist, create one
if (!dir.exists(paste0(root, "LiDAR_Data/merged_LiDAR_files"))){
  dir.create(paste0(root, "LiDAR_Data/merged_LiDAR_files"))
}
if (!dir.exists(paste0(root, "LiDAR_Data/plot_level_LiDAR"))){
  dir.create(paste0(root, "LiDAR_Data/plot_level_LiDAR"))
} 

merged_folder_path <- paste0(root, "LiDAR_Data/merged_LiDAR_files/")
plot_level_lidar_path <- paste0(root, "LiDAR_Data/plot_level_LiDAR/")

plot_data_table <- read.csv(paste0(root, "plot_data_table.csv"), colClasses = "character")

dpID="DP1.30003.001"
buffer=200

uniquesitemonthyear <- c(unique(plot_data_table$sitemonthyear))

for (sitemonthyear in uniquesitemonthyear){
  sliced_plot_data_table <- plot_data_table[plot_data_table$sitemonthyear == sitemonthyear,]
  
  notexist <- c()
  notexistdf <- data.frame()
  for (i in 1:nrow(sliced_plot_data_table)) {
    row <- sliced_plot_data_table[i,]
    site <- row$siteID
    monthyear <- row$monthyear
    year <- substr(monthyear,1,4)
    easting <- row$easting
    northing <- row$northing
    
    fname <- paste0(site, "_", monthyear, "_", easting, "_", northing, ".laz")
    fpath <- paste0(plot_level_lidar_path, fname)
    
    if (!file.exists(fpath)){
      notexist <- c(notexist, fname)
      notexistdf <- rbind(notexistdf, row)
    }
  }
  if (length(notexist) > 0){
    print(sitemonthyear)
    print(c(notexistdf$easting))
    print(c(notexistdf$northing))
    
    flag <- TRUE
    tryCatch(               
      # download Lidar data
      expr = {                     
        byTileAOP(dpID=dpID, site=site, 
                  year=year, easting=notexistdf$easting, northing=notexistdf$northing,
                  buffer=(buffer/2), check.size = FALSE,
                  savepath = paste0('./LiDAR_Data/',site,'/'))
      },
      error = function(e){ 
        flag <<- FALSE
        print(paste("There was an error message:- ", e))
      },
      warning = function(w){   
        flag <<- FALSE
        print(paste("There was a warning message:- ", w))
      }
    )
    
    print(flag)
    # code to merge .laz files
    if (flag){
      pt <- paste0(root, "LiDAR_Data/", site, "/", dpID, "/neon-aop-products/" ,year ,"/FullSite/")
      f <- list.files(pt)
      pt <- paste0(pt, f, "/")
      f <- list.files(pt)
      pt <- paste0(pt, f, "/L1/DiscreteLidar/ClassifiedPointCloud/")
      
      if (length(pt) > 1){pt <- pt[1]}
      
      merged_cmd <- paste0(lastools_bin_path, 
                           "lasmerge -i ", 
                           pt, 
                           "*.laz -o ", 
                           merged_folder_path, site, "_", monthyear, "_merged.laz")
      system(merged_cmd)
      
      # code to remove folder structure
      dir_path = paste0(root, 'LiDAR_Data/',site)
      shell( glue::glue("rmdir /s /q \"{dir_path}\" ") )
      
      
      # code to crop plots required from the merged plot
      for (i in 1:nrow(notexistdf)) {
        row <- notexistdf[i,]
        easting <- as.numeric(row$easting)
        northing <- as.numeric(row$northing)
        
        print(paste(site, monthyear, easting, northing))
        
        merged_cmd <- paste0(lastools_bin_path, 
                             "las2las -i ", 
                             paste0(merged_folder_path, site, "_", monthyear, "_merged.laz"), 
                             " -o ", 
                             paste0(plot_level_lidar_path, site, "_", monthyear, "_", easting, "_", northing, ".laz"), 
                             " -keep_xy ",
                             paste(as.character(easting - 100), as.character(northing - 100), as.character(easting + 100), as.character(northing + 100)))
        
        system(merged_cmd)
        
      }
      print("All plots cropped!")
      file.remove(paste0(merged_folder_path, site, "_", monthyear, "_merged.laz"))
    }
  }
}
