library(sf)
library(sp)
library(rgdal)
library(rnaturalearth)


#' @export
#'
#' @title Build a map grid spatialpolygons
#'
#' @description The function download from NaturalEarth database a worldmap raster.
#'  A grid is built from this raster.
#'  Then it returns a spatialpolygons object from `sp` package of each site of the grid.
#'
#' @param siteSize size of each cell of the grid in meter
#' @param projectionCRS Interface class to the PROJ projection and transformation system.
#' Default value is "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
#'
#' @return a spatialpolygons object of all sites of the grid into the projection `projectionCRS` given as argument
grid_spatialpolygons <- function (siteSize=260000,
                                projectionCRS="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") {
  # download the data
  coastlines <- rnaturalearth::ne_download(scale = 110, type = "coastline", category="physical", load=TRUE)
  coastlines.proj <- sp::spTransform(coastlines,
                                     sp::CRS(projectionCRS))
  coast.bbox <- sf::st_bbox(coastlines.proj)
  worldtangle <- sf::st_sfc(sf::st_polygon(list(rbind(c(coast.bbox[1],coast.bbox[2]),
                                                      c(coast.bbox[3],coast.bbox[2]),
                                                      c(coast.bbox[3],coast.bbox[4]),
                                                      c(coast.bbox[1],coast.bbox[4]),
                                                      c(coast.bbox[1],coast.bbox[2])))))
  grid <- sf::st_make_grid(worldtangle, cellsize=siteSize)
  grid.sp <- as_Spatial(grid)
  sp::proj4string(grid.sp) <- projectionCRS
  return(grid.sp)
}
