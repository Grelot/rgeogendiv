library(sp)
library(rgeos)



#' @export
#' @title specimen spatial points within spatial polygons site of grid
#' @param specimen.df a data.frame of specimen
#' @param grid.sp a spatialpolygons object of all sites of the grid
#' @param projectionCRS Interface class to the PROJ projection and transformation system.
#' Default value is the projection of the `grid.sp` parameter.
#' @return a matrix of presence/absence TRUE/FALSE of a specimen within a site
#' with column as site ID
#' and row as specimen spatial points ID
specimen_intersect_site <- function (specimen.df,
                                  grid.sp,
                                  projectionCRS= slot(grid.sp, 'proj4string')@projargs,
                                  fishbaseValidation = FALSE) {
  grid.sp.proj <- sp::spTransform(grid.sp,
                            sp::CRS(SRS_string = projectionCRS))
  specimen.sp.proj <- spatialpoints_bold_df(specimen.df,
                                       projectionCRS = projectionCRS,
                                       fishbaseValidation = fishbaseValidation)
  specimenIntersectSites <- rgeos::gIntersects(grid.sp.proj, specimen.sp.proj, byid = TRUE)
  return(specimenIntersectSites)

}
