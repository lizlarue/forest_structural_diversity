library(dplyr)
library(sp)
library(tibble)

#--------------------------------------------------------------------------------------------------------------------

# https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm/48838123#48838123

find_UTM_zone <- function(longitude, latitude) {
  
  # Special zones for Svalbard and Norway
  if (latitude >= 72.0 && latitude < 84.0 ) 
    if (longitude >= 0.0  && longitude <  9.0) 
      return(31);
  if (longitude >= 9.0  && longitude < 21.0)
    return(33)
  if (longitude >= 21.0 && longitude < 33.0)
    return(35)
  if (longitude >= 33.0 && longitude < 42.0) 
    return(37)
  
  (floor((longitude + 180) / 6) %% 60) + 1
}


find_UTM_hemisphere <- function(latitude) {
  
  ifelse(latitude > 0, "north", "south")
}

# returns a DF containing the UTM values, the zone and the hemisphere
longlat_to_UTM <- function(long, lat, units = 'm') {
  
  easting <- c()
  northing <- c()
  for(i in 1:length(long)){
    df <- data.frame(
      id = seq_along(long[i]), 
      x = long[i], 
      y = lat[i]
    )
    sp::coordinates(df) <- c("x", "y")
    
    hemisphere <- find_UTM_hemisphere(lat[i])
    zone <- find_UTM_zone(long[i], lat[i])
    
    sp::proj4string(df) <- sp::CRS("+init=epsg:4326") 
    CRSstring <- paste0(
      "+proj=utm +zone=", zone,
      " +ellps=WGS84 +datum=WGS84",
      " +", hemisphere,
      " +units=", units)
    if (dplyr::n_distinct(CRSstring) > 1L) 
      stop("multiple zone/hemisphere detected")
    
    res <- sp::spTransform(df, sp::CRS(CRSstring[1L])) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        zone = zone,
        hemisphere = hemisphere
      )
    easting <- c(easting, res$x)
    northing <- c(northing, res$y)
  
  }
  return(data.frame(easting, northing))
}

UTM_to_longlat <- function(utm_df, zone, hemisphere) {
  
  CRSstring <- paste0("+proj=utm +zone=", zone, " +", hemisphere)
  utmcoor <- sp::SpatialPoints(utm_df, proj4string = sp::CRS(CRSstring))
  longlatcoor <- sp::spTransform(utmcoor, sp::CRS("+init=epsg:4326"))
  tibble::as_data_frame(longlatcoor)
}

#--------------------------------------------------------------------------------------------------------------------

