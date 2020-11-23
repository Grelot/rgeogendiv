library(ggplot2)
library(ggpubr)

#' @export
#' @title Preliminary Analysis
#' @description From bold prepared query, generates a figure as ggplot object
#' Distribution of markercode
#' Accumulation of georeferenced GenBank accessions and frequency of submissions with mostly georeferenced data
#'
#' @param bold_res a data.frame with same fields than resBold$data$markercode and a supplementary column with DNA sequences as string
#' @return a ggplot object
preliminary_analysis <- function (bold_res) {
  markerdf <- marker_bold_res(bold_res)
  yeardf <- date_bold_res(bold_res)
  markergg <- ggplot2::ggplot(markerdf, aes(x = 2, y = markerscount, fill = markersname)) +
    ggplot2::geom_bar(stat = "identity", color = "white") +
    ggplot2::coord_polar(theta = "y", start = 0)+
    ggplot2::geom_text(aes(y = markerscount, label =  markerscount), color = "white")+
    ggplot2::theme_void()+
    ggplot2::xlim(0.5, 2.5)+
    ggplot2::ggtitle("(a) Distribution of markers")
  yeargg <- ggplot2::ggplot(yeardf, aes(x=year, y=cumsum(yearcount))) +
    ggplot2::geom_line() +
    ggplot2::geom_point()+
    ggplot2::ylab("Cumulative number of sequences")+
    ggplot2::theme_classic()+
    ggplot2::ggtitle("(b) Accumulation of accessions")
  gg <- ggpubr::ggarrange(markergg, yeargg, nrow = 2, ncol=2)
 return(gg)
}
