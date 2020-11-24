library(tidyverse)


#' @export
#' @title Cluster classification
#' @description Number of potentially informative intraspecific clusters
#' @param bold_res a data.frame with same fields than resBold$data$markercode and a supplementary column with DNA sequences as string
#' @return a data.frame of number of species and number of individual by species
cluster_classification_bold_res <- function(bold_res) {
  countIndvSpecies <- bold_res %>% group_by(species_name) %>% count
  table(cut_interval(countIndvSpecies$n,10))
  grid.sp <- grid_spatialpolygons(siteSize=260000,
                                  projectionCRS="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
  )
  ggplot2::cut_interval(1:100, 10)

}

