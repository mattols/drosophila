#
# Area Calculations
#
#

# # # # # # # #
# area of pixel
area_of_pixel <- function(pixel_size, center_lat) {
  # Calculate m^2 area of a wgs84 square pixel.
  # https://gis.stackexchange.com/questions/127165/more-accurate-way-to-calculate-area-of-rasters
  
  
  # Parameters:
  # pixel_size (float): length of side of pixel in degrees.
  # center_lat (float): latitude of the center of the pixel. Note this
  #     value +/- half the `pixel-size` must not exceed 90/-90 degrees
  #     latitude or an invalid area will be calculated.
  
  # Returns:
  # Area of square pixel of side length `pixel_size` centered at
  # `center_lat` in m^2.
  
  a <- 6378137  # meters
  b <- 6356752.3142  # meters
  e <- sqrt(1 - (b/a)^2)
  area_list <- c()
  for (f in c(center_lat + pixel_size/2, center_lat - pixel_size/2)) {
    zm <- 1 - e * sin(deg2rad(f))
    zp <- 1 + e * sin(deg2rad(f))
    area_list <- c(
      area_list,
      pi * b^2 * (
        log(zp/zm) / (2*e) +
          sin(deg2rad(f)) / (zp*zm))
    )
  }
  return(pixel_size / 360 * (area_list[1] - area_list[2]))
}

deg2rad <- function(deg) {
  return(deg * pi / 180)
}


## TEST ##
# Example usage:
pixel_size <- 0.01  # degrees
center_lat <- 45.0  # degrees
area <- area_of_pixel(pixel_size, center_lat)
print(area)
#
# for 0.5 grid along latitude
lat_ls <- seq(1,90, 5)
area_km2 <- unlist(lapply(lat_ls, function(x) area_of_pixel(0.5, x)/1000**2))
plot(area_km2~lat_ls, type="b")
# 
# compare projections
# wgs 84 vs equal area projection
path_name = "archive/data/wc2.1_10m/wc2.1_10m_bio_5.tif"
tmax_test = rast(path_name)
plot(tmax_test)
res(tmax_test)
tmax_test
#
# transform to km cellSize
cs <- cellSize(tmax_test, mask=TRUE, unit="km")
plot(cs)
# test
lats_test = data.frame(area = cellSize(tmax_test, unit="km")[,200])
# get coords - centroids (could get y max?)
lat_coords = crds(tmax_test, na.rm=FALSE)
lat01 = rast(tmax_test); lat01[] <- lat_coords[,2];names(lat01) = "Latitude"
lon01 = rast(tmax_test); lon01[] <- lat_coords[,1];names(lon01) = "Longitude"
plot(c(lat01,lon01))
# compare
# area_km2 <- unlist(lapply(lat_ls, function(x) area_of_pixel(0.5, x)/1000**2))
lats_test$lat = lat01[,200]$Latitude
head(lats_test)
plot(lats_test[,2:1], main = "Latitudinal distortion")
#
# Eckert IV - https://projectionwizard.org/#
"+proj=eck4 +lon_0=0 +datum=WGS84 +units=m +no_defs"

# terra
# area function test