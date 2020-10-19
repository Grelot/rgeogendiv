library(sp)
library(sf)
library(rgdal)
library(rgeos)

#' @export
#' @title nucleotide diversity by site on the worldmap grid
#' @description The function assign to each polygons of a grid the value of nucleotide diversity
#' as it is indicated in the data.frame `nucdivSites`.
#' @param nucdivSites a data.frame of nucleotide diversity by site and site grid's ids
#' @param grid.sp a spatialpolygons object of all sites of the grid
#' @return a data.frame of nucleotide diversity by site
nucleotide_diversity_grid <- function(nucdivSites, grid.sp) {
  grid.sp.nc <- as(grid.sp[nucdivSites$site.ids,], "SpatialPolygonsDataFrame")
  grid.sp.nc@data = data.frame(nucdivMean = nucdivSites$nucDivMean,
                               id = nucdivSites$site.ids.grid)
  return(grid.sp.nc)
}
