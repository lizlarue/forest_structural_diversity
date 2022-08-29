# This script will extract structural metrics data from plot level LiDAR files

library(lidR)
library(gstat)
library(neonUtilities)
library(dplyr)
library(rgeos)

lastools_bin_path <- "C:/Users/Lenovo/LAStools/LAStools/bin/"
root <- paste0(getwd(), "/")
plot_level_lidar_path <- paste0(root, "LiDAR_Data/plot_level_LiDAR/")

# read required files
plot_data_table <- read.csv(paste0(root, "plot_data_table.csv"), colClasses = "character")
plot_data_table$fname <- paste0(plot_data_table$siteID, "_", plot_data_table$monthyear, "_", plot_data_table$easting, "_", plot_data_table$northing, ".laz")

# create an empty dataframe
lidar_structural_metrics <- data.frame()
lidar_structural_metrics_fname <- c()
if (file.exists("lidar_structural_metrics.csv")){
  lidar_structural_metrics <- read.csv(paste0(root, "lidar_structural_metrics.csv"), colClasses = "character")
  lidar_structural_metrics_fname <- paste0(lidar_structural_metrics$siteID, "_", lidar_structural_metrics$monthyear, "_", lidar_structural_metrics$easting, "_", lidar_structural_metrics$northing, ".laz")
}

buffer_size <- 200 
plot_size <- 40

get_dtm_and_normalized_lidar_data <- function(cropped_lidar_data, res_level = 1, fn = kriging(k = 10L), use_class = c(1L, 2L), ct = 3){
  # gets cropped plot with buffer as input and returns normalized data as output
  flag <- TRUE
  loop <- 1
  while (flag & loop<(ct+1)){
    tryCatch(               
      expr = {       
        # referred from: https://github.com/r-lidar/lidR/issues/184
        # make a raster that encompass the point cloud
        layout <- raster(extent(cropped_lidar_data))
        res(layout) <- res_level
        # used:-  use_class = c(1L, 2L), to solve issue:- "Error: No ground points found. Impossible to compute a DTM."
        # referred from: https://github.com/r-lidar/lidR/issues/350
        dtm <<- grid_terrain(cropped_lidar_data, res = layout, kriging(k = 10L), use_class = c(1L, 2L))
        normalized_lidar_data <<- normalize_height(cropped_lidar_data, dtm)
        flag <- FALSE
      },
      error = function(e){ 
        #print(paste("There was an error message:- ", e))
        res_level <<- res_level + 1 
        loop <<- loop + 1
      }
    )
  }
}

get_cropped_plot_and_chm <- function(normalized_lidar_data, new_plot_size, easting, northing, ct){
  # get normalized cropped plot with buffer as input and returns normalized cropped plot without buffer with chm
  flag <- TRUE
  loop <- 1
  while (flag & loop<(ct+1)){
    tryCatch(               
      expr = {       
        cropped_plot <<- clip_rectangle(normalized_lidar_data, 
                                        xleft = (easting - (new_plot_size/2)), ybottom = (northing - (new_plot_size/2)),
                                        xright = (easting + (new_plot_size/2)), ytop = (northing + (new_plot_size/2)))
        cropped_plot@data$Z[cropped_plot@data$Z <= 0] <- NA 
        cropped_plot <<- cropped_plot
        chm <<- grid_canopy(cropped_plot, res = 1, dsmtin())  
        flag <- FALSE
      },
      error = function(e){ 
        #print(paste("There was an error message:- ", e))
        new_plot_size <<- new_plot_size + 10
        loop <<- loop + 1
      },
      warning = function(w){ 
        #print(paste("Warning on loop:- ", loop, w)) 
        new_plot_size <<- new_plot_size + 10
        loop <<- loop + 1
      }
    )
  }
}

structural_diversity_metrics <- function(normalized_lidar_data, new_plot_size, easting, northing) {
  # extract structural metrics from the ormalized data and the chm plot
  outerflag <- TRUE
  tryCatch(               
    expr = {                     
      
      get_cropped_plot_and_chm(normalized_lidar_data, new_plot_size, easting, northing, ct = 3)
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
      vert.sd <- cloud_metrics(lidar_data, sd(Z, na.rm = TRUE)) 
      sd.1m2 <- grid_metrics(lidar_data, sd(Z), 1) 
      sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
      Zs <- lidar_data@data$Z
      Zs <- Zs[!is.na(Zs)]
      entro <- entropy(Zs, by = 1) 
      gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
      GFP.AOP <- mean(gap_frac$gf) 
      LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
      VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
      VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
      out.plot <- data.frame(
        matrix(c(easting, northing, mean.max.canopy.ht,max.canopy.ht, 
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
      return(out.plot)
    },
    error = function(e){ 
      outerflag <<- FALSE
      print(paste("There was an error message:- ", e))
    }
  )
  
  if (!outerflag){
    return("error")
  }
}

error <- c()
not_available <- c()
for (i in 1:nrow(plot_data_table)){
  row <- plot_data_table[i,]
  if(!row$fname %in% lidar_structural_metrics_fname){
    # get plot data from the plot_data_table
    fn <- paste0(plot_level_lidar_path, row$fname)
    easting <- as.numeric(row$easting)
    northing <- as.numeric(row$northing)
    
    print(paste(i, row$fname, row$plotID))

    if (file.exists(fn))
    {
      # read file and figure out the outlier values
      lidar_data <- readLAS(fn, filter = "-drop_class 7")
      median_data <- median(lidar_data@data$Z)
      sd_data <- sd(lidar_data@data$Z)
      lower_limit <- median_data - (2*sd_data)
      upper_limit <- median_data + (3*sd_data)
      filter_data <- paste("-drop_z_below ",lower_limit," -drop_z_above ",upper_limit, " -drop_class 7")
      
      # read the file again with filtered values
      lidar_data <- readLAS(fn, filter = filter_data)
      
      # crop data for the plot along with buffer
      cropped_lidar_data <- clip_rectangle(lidar_data,
                                           xleft = (easting - (buffer_size/2)), ybottom = (northing - (buffer_size/2)),
                                           xright = (easting + (buffer_size/2)), ytop = (northing + (buffer_size/2)))
      
      # get normalized data from the function defined above
      get_dtm_and_normalized_lidar_data(cropped_lidar_data)
      # extract structural metrics from the data using the function defined above
      output <- structural_diversity_metrics(normalized_lidar_data, plot_size, easting, northing)
      
      plot_details <- row %>% select(siteID, monthyear, sitemonthyear, plotID)
      
      if (!output == "error"){
        record <- cbind(plot_details, output)
      }else{
        error <- c(error, row$fname)
        print(paste("Error in file:-", row$fname, "with PlotID:-", row$plotID))
        output <- data.frame(matrix(c(easting, northing, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), ncol = 15)) 
      }
      # write entries in the dataframe created above
      record <- cbind(plot_details, output)
      lidar_structural_metrics <- rbind(lidar_structural_metrics, record)
    
      }else{
      not_available <- c(not_available, row$fname)
      print(paste("File", fn, "not available"))
    }
  }
}
# write dataframe into a csv
write.table(lidar_structural_metrics, file = "lidar_structural_metrics.csv", sep = ",", row.names = FALSE)
