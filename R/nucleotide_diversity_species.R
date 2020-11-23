library(Biostrings)
library(tidyverse)
library(ape)
library(pegas)
library(future.apply)
library(data.table)

#' @export
#' @title function to align sequences with DECIPHER
align_sequence <- function(sequences, MinimumNumberOfSequencesBySpecies) {
  seqStringSet <- Biostrings::DNAStringSet(as.vector(sequences))
  ## check if at least 2 sequences
  if(length(seqStringSet@ranges) < MinimumNumberOfSequencesBySpecies) {
    return(NA)
  } else {
    #alignement <- muscle::muscle(seqStringSet, quiet=TRUE)
    seqStringSetNoGap <- DECIPHER::RemoveGaps(seqStringSet,
                                              removeGaps = "all",
                                              processors = 1)
    alignement <- DECIPHER::AlignSeqs(seqStringSetNoGap, verbose= FALSE)

    #alignement <- msa::msaMuscle(seqStringSet)
    return(alignement)
  }
}

#' @export
#' @title calculate nucleotide diversity
nuc_diversity_alignment <- function(alignment) {
  if(typeof(alignment) != "logical") {
    return(pegas::nuc.div(ape::as.DNAbin(alignment)))
  } else {
    return(NA)
  }
}


#' nucleotide diversity by species
#' arguments:
#' * specimen.df a dataframe with taxonomy+sequence dna information
#' * buffer.Within a dataframe with row as indiv seq and col as rls pts
#' with TRUE/FALSE presence/absence id of sequence specimen.df within a given buffer around a given rls pts
#' * rls.pts.index index of the RLS pts
#' it calculates
#' * alignment of sequences by species
#' * number of sequences by species
#' * nucleotide diversity
#' it returns
#' a dataframe with taxonomy and nucleotide diversity and number of sequences by species
#' species with not enough sequences to calculate nucleotide diversity are removed
#' @export
#' @title Nucleotide diversity by species
#' @param site.ids index of the site of the grid
#' @param specimen.df a data.frame with taxonomy+sequence dna information
#' @param sequenceIntersectSites a matrix of presence/absence TRUE/FALSE of a specimen within a site
#' with column as site ID
#' and row as specimen spatial points ID
#' @param MinimumNumberOfSequencesBySpecies an integer minimum number of sequences representatives of a species in a site
#' @return a data.frame with number of individuals and nucleotide diversity for each species in the site.
species_nucleotide_diversity_site <- function(site.ids, specimen.df, sequenceIntersectSites, MinimumNumberOfSequencesBySpecies, fishbaseValidation=FALSE) {
  seqboldwithin <- specimen.df[which(sequenceIntersectSites[, site.ids] == TRUE),]
  seqboldwithinGroupSpecies <- seqboldwithin %>% dplyr::group_by(species_name)
  seqboldwithinUniqSpecies <- seqboldwithinGroupSpecies %>% dplyr::slice(1) %>% dplyr::ungroup()
  numberSequencesAli <- seqboldwithinGroupSpecies %>% dplyr::tally()
  groupSeq <- seqboldwithinGroupSpecies %>% dplyr::group_map(~Biostrings::DNAStringSet(as.vector(.x$sequence)))
  alignSeq <- purrr::map(groupSeq, ~ align_sequence(.x, MinimumNumberOfSequencesBySpecies))
  nucDivAll <- alignSeq %>% purrr::map(nuc_diversity_alignment) %>% unlist()
  if(fishbaseValidation) {
    seqboldwNucDiv <- data.table::data.table(species_name=seqboldwithinUniqSpecies$species_name,
                                             fishbase_species_name=seqboldwithinUniqSpecies$fishbase_species_name,
                                             genus_name=seqboldwithinUniqSpecies$genus_name,
                                             family_name=seqboldwithinUniqSpecies$family_name,
                                             order_name=seqboldwithinUniqSpecies$order_name,
                                             class_name=seqboldwithinUniqSpecies$class_name,
                                             nucDiv=nucDivAll,
                                             number_individuals=numberSequencesAli$n,
                                             site.ids=rep(site.ids,length(nucDivAll))
  )
  } else {
    seqboldwNucDiv <- data.table::data.table(species_name=seqboldwithinUniqSpecies$species_name,
                                             genus_name=seqboldwithinUniqSpecies$genus_name,
                                             family_name=seqboldwithinUniqSpecies$family_name,
                                             order_name=seqboldwithinUniqSpecies$order_name,
                                             class_name=seqboldwithinUniqSpecies$class_name,
                                             nucDiv=nucDivAll,
                                             number_individuals=numberSequencesAli$n,
                                             site.ids=rep(site.ids,length(nucDivAll))
  )
  }
  seqboldwNucDiv <- seqboldwNucDiv[which(!is.na(seqboldwNucDiv$nucDiv)),]
  return(seqboldwNucDiv)
}


#'
#' @export
#' @title Number of sequences in a site
#' @param sequenceIntersectSites a matrix of presence/absence TRUE/FALSE of a specimen within a site
#' with column as site ID
#' @return a list with the number of sequences within a site
number_sequences_by_points <- function(sequenceIntersectSites) {
  numberOfSequencesBySite <- data.frame(apply(sequenceIntersectSites, 2, function(x) length(which(x))))
  names(numberOfSequencesBySite) <- "numberOfSequencesSite"
  return(numberOfSequencesBySite)
}


#' @export
#' @title nucleotide diversity by species applied to all sites
#' @param specimen.df a data.frame of specimen
#' @param sequenceIntersectSites a matrix of presence/absence TRUE/FALSE of a specimen within a site
#' with column as site ID
#' and row as specimen spatial points ID
#' @param MinimumNumberOfSequencesBySpecies an integer minimum number of sequences representatives of a species in a site
#' @return a dataframe with taxonomy and nucleotide diversity and number of sequences by species
#' species with not enough sequences to calculate nucleotide diversity are removed
nucleotide_diversity_species <- function(specimen.df, sequenceIntersectSites, MinimumNumberOfSequencesBySpecies, fishbaseValidation= FALSE) {
  ## keep only buffer with at least 3 sequences
  if(fishbaseValidation) {
    specimen.df.selectedCol <- data.table::data.table(species_name=specimen.df$species_name,
                                                fishbase_species_name=specimen.df$fishbase_species_name,
                                                genus_name=specimen.df$genus_name,
                                                family_name=specimen.df$family_name,
                                                order_name=specimen.df$order_name,
                                                class_name=specimen.df$class_name,
                                                sequence=specimen.df$sequence
    )
  } else {
    specimen.df.selectedCol <- data.table::data.table(species_name=specimen.df$species_name,
                                                      genus_name=specimen.df$genus_name,
                                                      family_name=specimen.df$family_name,
                                                      order_name=specimen.df$order_name,
                                                      class_name=specimen.df$class_name,
                                                      sequence=specimen.df$sequence
    )
  }
  numberOfSequencesBySite <- number_sequences_by_points(sequenceIntersectSites)
  ## keep sites with at least 2 sequences
  site.ids <- as.vector(which(numberOfSequencesBySite > 2))
  ## calculate nucdiv for each points
  nucDivPts <- lapply(site.ids,
                      FUN = function(i) {
                        species_nucleotide_diversity_site(i,
                                                         specimen.df.selectedCol,
                                                         sequenceIntersectSites,
                                                         MinimumNumberOfSequencesBySpecies
                        )
                      })


  nucDivPts.df <- dplyr::bind_rows(nucDivPts)
  return(nucDivPts.df)
}





