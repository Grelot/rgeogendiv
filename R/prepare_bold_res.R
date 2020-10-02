library("tidyverse")

#' @export
#' @title Prepare raw BOLD specimen and DNA sequences dataset
#'
#' @description Filter and mutate BOLD dataset to produce a curated data.frame
#'  with rows as individual specimen and columns as specimen information.
#'  It adds a new column sequence with fasta sequences as string.
#' @param bold_res A list of lists returned by bold_seqspec command from bold package. They are 2 lists:
#' * Specimen information (spatial coordinates, taxonomy...)
#' * DNA barcode sequences
#' @param marker_code If not empty, only specimen with field data$markercode matching marker_code names are kept. By default, markercode filtering is not applied.
#' @param species_names a logical value indicating whether specimen with no species name information should be removed.
#' @param coordinates a logical value indicating whether specimen with no latitude or longitude spatial coordinates information should be removed.
#' @param ambiguities a logical value indicating whether specimen with DNA sequence containing IUPAC ambiguities should be removed.
#' @param min_length numeric. Minimum length of a DNA sequence. Default value: 0 bp.
#' @param max_length numeric. Maximum length of a DNA sequence. Default value: 800 bp
#' @return a data.frame with same fields than resBold$data$markercode and a supplementary column with DNA sequences as string.
prepare_bold_res <- function (bold_res, marker_code="",
                              species_names=TRUE,
                              coordinates=TRUE,
                              ambiguities=TRUE,
                              min_length=0,
                              max_length=800) {
  if (marker_code != "") {
    if(!marker_code %in% bold_res$data$markercode) {
      print(paste0("WARNING: markercode ",marker_code ," not found in colon 'markercode' of ",bold_res))
    }
  }
  filteredBold <- bold_res$data %>%
    { if (species_names) dplyr::filter(., species_name != "") } %>%
    { if (marker_code != "") dplyr::filter(., markercode == marker_code) } %>%
    { if (coordinates) dplyr::filter(., !is.na(lat)) %>% dplyr::filter(., !is.na(lon)) } %>%
    dplyr::mutate(sequence=as.character(bold_res$fasta[processid])) %>%
    dplyr::filter(., nchar(sequence) > min_length) %>%
    dplyr::filter(., nchar(sequence) < max_length) %>%
    { if (ambiguities) dplyr::filter(., !grepl("I",sequence)) }
  return(filteredBold)
}
