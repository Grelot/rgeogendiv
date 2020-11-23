#' @export
#' @title Return a data.frame of markers and count from BOLD dataset
#'
#' @description after preparing bold dataset query using `prepare_bold_res` function, this function will extract a data.frame of markers.
#' @param bold_res a data.frame with same fields than resBold$data$markercode and a supplementary column with DNA sequences as string.
#' @return a table of makercode and number of occurences into the BOLD dataset
marker_bold_res <- function (bold_res) {
  tablemarker <- table(bold_res$markercode)
  markernames <- names(tablemarker)
  markercount <- as.vector(tablemarker)
  dfmarker = data.frame("markersname" = markernames,
                  "markerscount" = markercount)
  return(dfmarker)
}
