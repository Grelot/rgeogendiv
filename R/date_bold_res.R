

#' @export
#' @title extract date from a string
#' @description From a string formated as follow: "2006-01-26 19:24:04|2006-01-27 01:29:30". This function return the first year.
#' @param datestr a string of a list of dates
#' @return a string of a year
extract_date_from_str <- function (datestr) {
  if (datestr != "")  {
    listDates <- unlist(strsplit(unlist(datestr), "|", fixed=TRUE))
    splitDates <- unlist(strsplit(listDates[1],"-"))
    year <- splitDates[1]
    return(year)
  } else {
    return(NA)
  }
}

#' @export
#' @title run date of bold specimen
#' @description from the prepared data.fram of query BOLD dataset, the function returns a data.frame of number of cumulative contribtion
#' of new specimen for each year
#' @param bold_res a data.frame with same fields than resBold$data$markercode and a supplementary column with DNA sequences as string.
#' @return a data.frame with years and cumulative number of BOLD specimenby year
date_bold_res <- function (bold_res) {
  yearList <- unlist(lapply(bold_res$run_dates, FUN = function(X) { extract_date_from_str(X) }))
  countNA <- length(which(is.na(yearList)))
  tableYear <- table(yearList)
  years <- as.numeric(names(tableYear))
  allYears <- seq(min(years), max(years),1)
  allYearsCount <- c()
  for(y in allYears) {
    idYear <- which(years == y)
    if (length(idYear) > 0) {
      allYearsCount <- c(allYearsCount, as.integer(tableYear[idYear]))
    } else {
      allYearsCount <- c(allYearsCount, 0)
    }
  }
  dfYear <- data.frame("year" = allYears,
                       "yearcount" = allYearsCount)
  return(dfYear)
}
