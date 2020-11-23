library(sp)

#' @export
#' @title Transform curated BOLD data.frame into spatialpoints object from `sp` package
#' @description Given a BOLD individual georeferenced sequences dataframe and a projection CRS,
#' this function will return a spatialpoints object with all coordinates of BOLD specimen.
#' @param dfBold A data.frame with specimen sequences as row and 81 descriptors as column.
#' Mandatory fields:
#' lon,
#' lat,
#' species_name,
#' fishbase_species_name,
#' genus_name,
#' family_name,
#' order_name,
#' class_name,
#' sequence
#' @param projectionCRS Interface class to the PROJ projection and transformation system.
#' Default value is "+init=epsg:3347"
#' @return a spatialpoints object with taxonomic information and DNA sequence
#'
spatialpoints_bold_df <- function (dfBold, projectionCRS="+init=epsg:3347", fishbaseValidation=FALSE) {
  ## define coordinates of BOLD points as spatialpoints and projection in meter unit
  bold.coo <- data.frame(lon=dfBold$lon,
                         lat= dfBold$lat
  )
  ## define information species/sequence related to each point BOLD
  if(fishbaseValidation) {
    bold.info <- data.frame(species_name=dfBold$species_name,
                           fishbase_species_name=dfBold$fishbase_species_name,
                           genus_name=dfBold$genus_name,
                           family_name=dfBold$family_name,
                           order_name=dfBold$order_name,
                           class_name=dfBold$class_name,
                           sequence=dfBold$sequence
     )
  } else {
    bold.info <- data.frame(species_name=dfBold$species_name,
                            genus_name=dfBold$genus_name,
                            family_name=dfBold$family_name,
                            order_name=dfBold$order_name,
                            class_name=dfBold$class_name,
                            sequence=dfBold$sequence
    )
  }
  bold.pts <- sp::SpatialPointsDataFrame(bold.coo,
                                         data=bold.info,
                                         proj4string=sp::CRS("+init=epsg:4326")
  )
  bold.sp <- sp::spTransform(bold.pts,
                             sp::CRS(projectionCRS)
  )
  return(bold.sp)
}

