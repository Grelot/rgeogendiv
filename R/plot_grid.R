library(ggplot2)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(rnaturalearth)



#' @export
#' @title labels of Intervals
#' @description A list of label of intervals is produced based on a
#' numeric vector of probabilities with values in [0,1].
#' The smallest observation corresponds to a probability of 0
#'  and the largest to a probability of 1.
#' @param gd list of genetic diversity values
#' @param intervalles numeric vector of probabilities with values in [0,1].
#' @return a list of interval labels
legend_qt_GD <- function(gd, intervalles=c(0,0.1,0.2,0.3,0.4,0.6,1)) {
  gd.qt <- quantile(gd, probs=intervalles)
  gd.qt.format <- lapply(gd.qt, FUN = function(X) {
                                                    format(as.double(X[1]),
                                                           nsmall=4,
                                                           scientific = FALSE)
                                                  }
                        )
  gd.qt.ch <- as.character(gd.qt.format)
  gd.intervalles = paste("< ",gd.qt.ch[1])
  for( i in 1:(length(gd.qt.ch)-1)) {
    gd.intervalle = paste(gd.qt.ch[i],gd.qt.ch[i+1],sep=" - ")
    gd.intervalles=c(gd.intervalles, gd.intervalle)
  }
  gd.intervalles= c(gd.intervalles, paste("> ",gd.qt.ch[length(gd.qt.ch)]))
  return(gd.intervalles)
}


#' @export
#' @title labels of Limits
#' @description A list of label of intervals is produced based on a
#' list of limits.
#' For instance if you provide list of limits [0.1, 0.2, 0.3, 1]
#' It will returns a list ['< 0.1', '0.1 - 0.2', '0.2 - 0.3', '> 1'].
#' @param limites numeric vector of numeric limits
#' default value = c(0.001,0.002,0.003,0.005,0.01,0.025)
#' @return a list of intervals labels
legend_limit_GD <- function(limites=c(0.001,0.002,0.003,0.005,0.01,0.025)) {
  limites.ch <- as.character(limites)
  limites.intervalles = paste("< ",limites.ch[1])
  for( i in 1:(length(limites.ch)-1)) {
    limites.intervalle = paste(limites.ch[i],limites.ch[i+1],sep=" - ")
    limites.intervalles=c(limites.intervalles, limites.intervalle)
  }
  limites.intervalles= c(limites.intervalles, paste("> ",limites.ch[length(limites.ch)]))
  return(limites.intervalles)
}

#' @export
#' @title Colors of numeric value
#' @description attributes a hex color value to a numeric value based on
#' intervals of numeric value
#' @param gd list of numeric values
#' @param limites numeric vector of numeric limits values
#' default value = c(0.001,0.002,0.003,0.005,0.01,0.025)
#' @return a list of color hex value corresponding to the interval of each numeric value
color_GD <- function(gd, limites=c(0.001,0.002,0.003,0.005,0.01,0.025)) {
  limites.intervalles <- legend_limit_GD(limites)
  colfunc <- colorRampPalette(c("#3333A2","#3333FF","#33CBFF","#33FFFF","#FFDF33","#FFA333","#FF3333"))
  color.intervalles <- colfunc(length(limites))
  gdc=gd
  gdc[which(gd<limites[1])] = color.intervalles[1]
  for( i in 1:(length(limites)-1)) {
    gdc[which(gd>=limites[i] & gd<limites[i+1])] = color.intervalles[i]
  }
  gdc[which(gd>limites[length(limites)])] = color.intervalles[length(limites)]
  return(gdc)
}

#' @export
#' @title Labels of numeric value
#' @description attributes a label string to a numeric value based on
#' intervals of numeric value
#' @param gd list of numeric values
#' @param limites numeric vector of numeric limits values
#' default value = c(0.001,0.002,0.003,0.005,0.01,0.025)
#' @return a list of labels corresponding to the interval of each numeric value
tag_GD <- function(gd, limites=c(0.001,0.002,0.003,0.005,0.01,0.025)) {
  limites.intervalles <- legend_limit_GD(limites)
  gdc=gd
  gdc[which(gd<limites[1])] = limites.intervalles[1]
  for( i in 1:(length(limites)-1)) {
    gdc[which(gd>=limites[i] & gd<limites[i+1])] = limites.intervalles[i]
  }
  gdc[which(gd>limites[length(limites)])] =limites.intervalles[length(limites)]
  return(gdc)
}

#' @export
#' @title Print worldmap of mean species genetic diversity
#' @description downloads coastlines, rivers and lakes shapefile from `naturalearth`.
#' print sites of the grid worldmap with colors and labels of their genetic values.
#' @param nucdivGrid a data.frame of nucleotide diversity by site
#' @param limites numeric vector of numeric limits values
#' default value = c(0.001,0.002,0.003,0.005,0.01,0.025)
#' @return a `ggplot2` object of the worldmap of mean species genetic diversity
plot_grid <- function(nucdivGrid, limites=c(0.001,0.002,0.003,0.005,0.01,0.025)) {
  coastlines <- rnaturalearth::ne_download(scale = 'medium',
                                           type = "coastline",
                                           category="physical",
                                           load=TRUE)
  coastlines.proj <- sp::spTransform(coastlines,
                                     sp::CRS(sp::proj4string(nucdivGrid)))
  coastlines.fort <- ggplot2::fortify(coastlines.proj)
  riverlines <-  rnaturalearth::ne_download(scale = 'medium',
                                            type = "rivers_lake_centerlines",
                                            category="physical",
                                            load=TRUE)
  riverlines.proj <- sp::spTransform(riverlines,
                                     sp::CRS(sp::proj4string(nucdivGrid)))
  riverlines.fort <- ggplot2::fortify(riverlines.proj)

  lakelines <-  rnaturalearth::ne_download(scale = 'medium',
                                            type = "lakes",
                                            category="physical",
                                            load=TRUE)
  lakelines.proj <- sp::spTransform(lakelines,
                                     sp::CRS(sp::proj4string(nucdivGrid)))
  lakelines.fort <- ggplot2::fortify(lakelines.proj)
  sites <- ggplot2::fortify(nucdivGrid, region = "id")
  sites.df <- merge(sites, nucdivGrid@data, by = "id")
  #sites.df$nucDivColor <- color_GD(nucdivGrid@data$nucdivMean, limites)
  sites.df$nucDivInterval <- tag_GD(nucdivGrid@data$nucdivMean, limites)
  namedVal <- color_GD(nucdivGrid@data$nucdivMean, limites)
  names(namedVal) <- tag_GD(nucdivGrid@data$nucdivMean, limites)
  #levels(sites.df$nucDivInterval) <- sites.df$nucDivColor
  gg <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data=coastlines.fort, ggplot2::aes(x=long,
                                                             y=lat,
                                                             group=group),
                                        fill="grey", color= "black",size=0.25) +
    ggplot2::geom_polygon(data=lakelines.fort, ggplot2::aes(x=long,
                                                            y=lat,
                                                            group=group),
                                        fill="white", color= "black",size=0.25) +
    ggplot2::geom_path(data=riverlines.fort, ggplot2::aes(x=long,
                                                          y=lat,
                                                          group=group),
                                        color="white",size=0.35) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())+
    ggplot2::geom_polygon(data =sites.df, ggplot2::aes(x=long,
                                                       y=lat,
                                                       group = group,
                                     fill = factor(nucDivInterval)),
                                     colour="white")+
    ggplot2::scale_fill_manual(values = namedVal)+
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Genetic diversity"))
  return(gg)
}




