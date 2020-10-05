library(sf)
library(sp)
library(rgdal)
library(rnaturalearth)


#' @export
#'
#' @title Build a map grid
#'
#' @description The function download from NaturalEarth database a worldmap raster.
#'  A grid is built from this raster.
#'  Then the list of centroids of each site of the grid is returned.
#'
#' @param siteSize size of each cell of the grid in meter
#' @param projectionCRS Interface class to the PROJ projection and transformation system.
#' Default value is "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
#'
#' @return a spatialpoint object of all sites of the grid into the projection `projectionCRS` given as argument
grid_coordinates <- function (siteSize=260000,
                              projectionCRS="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") {
  # download the data
  coastlines <- rnaturalearth::ne_download(scale = 110, type = "coastline", category="physical", load=TRUE)
  coastlines <- rgdal::readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
  coastlines.proj <- sp::spTransform(coastlines,
                                     sp::CRS(projectionCRS))
  coast.bbox <- sf::st_bbox(coastlines.proj)
  worldtangle <- sf::st_sfc(sf::st_polygon(list(rbind(c(coast.bbox[1],coast.bbox[2]),
                                              c(coast.bbox[3],coast.bbox[2]),
                                              c(coast.bbox[3],coast.bbox[4]),
                                              c(coast.bbox[1],coast.bbox[4]),
                                              c(coast.bbox[1],coast.bbox[2])))))
  grid <- sf::st_make_grid(worldtangle, cellsize=siteSize)
  gridCentroid <- sf::st_centroid(grid)
  gridCentroid.sp <- as_Spatial(gridCentroid)
  sp::proj4string(gridCentroid.sp) <- projectionCRS
  return(gridCentroid.sp)
}
