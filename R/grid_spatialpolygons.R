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
  #coastlines <- rnaturalearth::ne_download(scale = 110, type = "coastline", category="physical", load=TRUE)
  #coastlines.proj <- sp::spTransform(coastlines, behrman)
  #coast.bbox <- sf::st_bbox(coastlines.proj)
  coast.bbox <- data.frame(xmin =-17367530,
                           xmax =17367530,
                           ymin = -7320487,
                           ymax = 7296713)
  worldtangle <- sf::st_sfc(sf::st_polygon(list(rbind(c(coast.bbox$xmin,coast.bbox$ymin),
                                                      c(coast.bbox$xmax,coast.bbox$ymin),
                                                      c(coast.bbox$xmax,coast.bbox$ymax),
                                                      c(coast.bbox$xmin,coast.bbox$ymax),
                                                      c(coast.bbox$xmin,coast.bbox$ymin)
                                                      ))))
  grid <- sf::st_make_grid(worldtangle, cellsize=siteSize)
  grid.sp <- sf::as_Spatial(grid)
  sp::proj4string(grid.sp) <- sp::CRS(SRS_string = projectionCRS)
  #slot(grid.sp, 'proj4string')@projargs <- projectionCRS
  return(grid.sp)
}
