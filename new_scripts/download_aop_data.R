library(neonUtilities)
library(rhdf5)
library(plyr)
library(ggplot2)

root <- paste0(getwd(), "/")

# if folder does not exist, create one
if (!dir.exists(paste0(root, "AOP_Data"))){
  dir.create(paste0(root, "AOP_Data"))
}
if (!dir.exists(paste0(root, "AOP_Data/plot_level_AOP"))){
  dir.create(paste0(root, "AOP_Data/plot_level_AOP"))
}

# read plot data table
plot_data_table <- read.csv(paste0(root, "plot_data_table.csv"), colClasses = "character")

dpID = "DP3.30006.001"
nBands <- 426
buffer = 200

plot_level_aop_path <- paste0(root,'AOP_Data/plot_level_AOP/')

all_file_names <- c()
download_error <- c()
other_errors <- c()

for (i in 1:nrow(plot_data_table)){ 
  row <- plot_data_table[i,]
  site <- row$siteID
  monthyear <- row$monthyear
  year <- substr(monthyear,1,4)
  easting <- row$easting
  northing <- row$northing
  plotId <- row$plotID
  
  fname <- paste0(site, "_", monthyear, "_", easting, "_", northing, ".h5")
  fpath <- paste0(plot_level_aop_path, fname)
  dir_path = paste0(root, 'AOP_Data/', site)
  print(paste(i, fname))
  all_file_names <- c(all_file_names, fname)
  
  easting <- as.character(round(as.numeric(easting))) # rounding off
  northing <- as.character(round(as.numeric(northing))) # rounding off
  
  # if file does not exist, download from NEON
  if (!file.exists(fpath)){
    
    flag <- TRUE
    tryCatch(               
      
      expr = {                     
        byTileAOP(dpID=dpID, site=site, 
                  year=year, easting=easting, northing=northing,
                  buffer=(buffer/2), check.size = FALSE,
                  savepath = paste0('./AOP_Data/',site))
        
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
    
    if (flag){
      
      # get file location
      pt <- paste0(root, "AOP_Data/", site, "/", dpID, "/neon-aop-products/" ,year ,"/FullSite/")
      f <- list.files(pt)
      pt <- paste0(pt, f, "/")
      f <- list.files(pt)
      pt <- paste0(pt, f, "/")
      f <- list.files(pt)
      pt <- paste0(pt, f[1], "/")
      f <- list.files(pt)
      pt <- paste0(pt, f, "/Reflectance/")
      f <- list.files(pt)
      
      tryCatch(               
        
        expr = {
                  # if more than one file are downloaded, get a list of easting and northing 
                  # values for each file based on the file names, get unique values and sort them.
                  if (length(f) > 1){
                    
                    easting_list <- c()
                    northing_list <- c()
                    
                    for (i in 1:length(f)){
                      fn <- strsplit(f[i], "_")
                      e <- fn[[1]][5]
                      n <- fn[[1]][6]
                      easting_list <- c(easting_list, e)
                      northing_list <- c(northing_list, n)
                    }
                    
                    easting_list <- sort(unique(easting_list))
                    northing_list <- sort(unique(northing_list))
                    
                    # create an array of size buffer_size*buffer_size*total_no_of_bands filled with zero
                    bands_data <- array(0, dim=c(buffer,buffer,nBands))
                    
                    # now for each file, we figure out the area in the tile which is present in our buffer zone
                    # and then we slice the data from the entire tile for that particular area
                    for (i in 1:length(f)){
                      print(f[i])
                      fn <- strsplit(f[i], "_")
                      site <- fn[[1]][3]
                      e <- fn[[1]][5]
                      n <- fn[[1]][6]
                      
                      xmin <- round(max(c(as.numeric(e), (as.numeric(easting)-(buffer/2)))) - as.numeric(e)) + 1
                      xmax <- round(min(c(as.numeric(e)+1000, (as.numeric(easting)+(buffer/2)))) - as.numeric(e))
                      ymin <- round(max(c(as.numeric(n), (as.numeric(northing)-(buffer/2)))) - as.numeric(n)) + 1
                      ymax <- round(min(c(as.numeric(n)+1000, (as.numeric(northing)+(buffer/2)))) - as.numeric(n))
                      
                      #print(paste(xmin, xmax, ymin, ymax))
                      
                      cropped_bands_data <- h5read(paste0(pt, f[i]), paste0("/",site,"/Reflectance/Reflectance_Data"), index=list(NULL,xmin:xmax,ymin:ymax))
                      h5closeAll()
                      
                      # rearrange the matrix to shape that would be consistent with buffer_size*buffer_size*total_no_of_bands
                      cropped_bands_data <- aperm(cropped_bands_data, c(2,3,1))
                      
                      # find the appropriate position for this data to be filled in the zero array created above 
                      if (as.numeric(e) < (as.numeric(easting)-(buffer/2))){
                        x1 = 1
                        x2 = xmax-xmin+1
                      }else{
                        x1 = buffer-xmax+1
                        x2 = buffer
                      }
                      
                      if (as.numeric(n) < (as.numeric(northing)-(buffer/2))){
                        y1 = 1
                        y2 = ymax-ymin+1
                      }else{
                        y1 = buffer-ymax+1
                        y2 = buffer
                      }
                      
                      bands_data[x1:x2, y1:y2,] <- cropped_bands_data
                    }
                    
                    pt <- paste0(pt, f[1])
                  
                  }else{
                    
                    # if only one file is downloaded, read the data file, and get the following data from that file
                    pt <- paste0(pt, f)
                    site <- strsplit(f, "_")[[1]][3]
                    
                    reflInfo <- h5readAttributes(pt,  paste0("/",site, "/Reflectance/Reflectance_Data"))
                    nRows <- reflInfo$Dimensions[1]
                    nCols <- reflInfo$Dimensions[2]
                    nBands <- reflInfo$Dimensions[3]
                    
                    xmin <- as.numeric(reflInfo$Spatial_Extent_meters[1])
                    xmax <- as.numeric(reflInfo$Spatial_Extent_meters[2])
                    ymin <- as.numeric(reflInfo$Spatial_Extent_meters[3])
                    ymax <- as.numeric(reflInfo$Spatial_Extent_meters[4])
                    
                    x = round(as.numeric(easting) - xmin)
                    y = round(as.numeric(northing) - ymin)
                    
                    # find out the coordinates for the plot centroid and the buffer region and slice data from the downloaded file
                    bands_data <- h5read(pt, paste0("/",site, "/Reflectance/Reflectance_Data"),index=list(NULL,(x-(buffer/2)+1):(x+(buffer/2)),(y-(buffer/2)+1):(y+(buffer/2)))) 
                    bands_data <- aperm(bands_data, c(2,3,1))
                  }
                  
                  # by this time, we have the data for the required plot centroid along with the buffer data
                  # we now get the wavelength array, and the scale factor.
          
                  wavelengths <- h5read(pt, paste0("/",site, "/Reflectance/Metadata/Spectral_Data/Wavelength"))
                  reflInfo <- h5readAttributes(pt,  paste0("/",site, "/Reflectance/Reflectance_Data"))
                  reflMetadata <- h5readAttributes(pt, paste0("/",site,"/Reflectance"))
                  h5closeAll()
                  
                  scaleFact <- reflInfo$Scale_Factor
                  # we scale the reflectance data matrix by this scale factor
                  bands_data <- bands_data/as.vector(scaleFact)
                  
                  # remove reflectance data for the atmospheric absorption zone, which has noisy values
                  ab1 <- reflMetadata$Band_Window_1_Nanometers
                  ab2 <- reflMetadata$Band_Window_2_Nanometers
                  
                  bands_data[,,wavelengths>ab1[1] & wavelengths<ab1[2]] <- NA
                  bands_data[,,wavelengths>ab2[1] & wavelengths<ab2[2]] <- NA
                  
                  # Write the final reflectance matrix into a .hf format file with all relevant data
                  h5write(bands_data, fpath, "/reflectance")
                  h5write(wavelengths, fpath, "/wavelengths")
                  h5write(scaleFact, fpath, "/scaleFact")
                  h5write(easting, fpath, "/easting")
                  h5write(northing, fpath, "/northing")
                  h5write(paste(site, monthyear), fpath, "/site_monthyear")
                  h5write(plotId, fpath, "/plotId")
                  
                  # code to remove folder structure
                  unlink(dir_path, recursive=TRUE)
        },
        error = function(e){ 
          #flag <<- FALSE
          other_errors <- c(other_errors, fname)
          print(paste("There was an error message:- ", e))
        },
        warning = function(w){   
          #flag <<- FALSE
          other_errors <- c(other_errors, fname)
          print(paste("There was a warning message:- ", w))
        }
      )
    }else{
      print(paste('Download error for file: ', fname))
      download_error <- c(download_error, fname)
    }
  }
}

#-------------------------------------------------------------------------

# To plot reflectance for x,y pixel in a matrix
reflectance_plot_by_pixel <- function(reflectance_data_matrix, x, y){
  aPixeldf <- adply(reflectance_data_matrix[x,y,],c(1))
  aPixeldf <- aPixeldf[2]
  aPixeldf$Wavelength <- wavelengths
  names(aPixeldf) <- c('ScaledReflectance','Wavelength')
  ggplot(data=aPixeldf)+
    geom_line(aes(x=Wavelength, y=ScaledReflectance))+
    xlab("Wavelength (nm)")+
    ylab("Reflectance")
}

# Example of plotting the reflecatnce values for a random coordinate in the reflectance matrix

#random_file_number <- 8
#temp_h5_path <- paste0(root, 'AOP_Data/plot_level_AOP/', all_file_names[random_file_number])
#View(h5ls(temp_h5_path,all=T))

#temp_reflectance_matrix <- h5read(temp_h5_path, "/reflectance")
#reflectance_plot_by_pixel(temp_reflectance_matrix, 100, 80)
