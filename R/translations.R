#' internal function: function factory for 'translating' numbers to their
#'  corresponding strings/messages
#'  
#' @param translation character vector specifying the strings/messages
#' @param setNAZero boolean specifying whether NA values should be converted to
#'  zero
#' @param zeroIsFirst boolean specifying whether the first element is zero or 1
#'  (for dealing with zero values). If TRUE then all values are increased by 1
#'  
#' @return a function which converts a (vector of) number(s) to the specified
#'  strings/messages
#' @noRd
translateInfo <- function(translation, setNAZero = TRUE, zeroIsFirst = TRUE){
  if (zeroIsFirst){
    if (setNAZero){
      function(info){
        info[is.na(info)] <- 0
        return(translation[info+1])
      }
    } else {
      function(info){
        return(translation[info+1])
      }
    }
  } else {
    if (setNAZero){
      function(info){
        info[is.na(info)] <- 0
        return(translation[info])
      }
    } else {
      function(info){
        return(translation[info])
      }
    }
  }
}

#' function for 'translation' of the psmAmbiguity values (1..5) in the psmTable
#'  to words (like in Proteome Discoverer). <...> --> means not encountered/
#'  undefined/no inference
#' 
#' @param info integer vector to be 'translated'
#' 
#' @return character vector (the translation)
#' @export
psmAmbiguity <- translateInfo(
  translation = c("<Not Considered>",
                  "Rejected",
                  "<Ambiguous>",
                  "Selected",
                  "Unambiguous"),
  setNAZero = FALSE, zeroIsFirst = FALSE)

#' function for 'translation' of the isMasterProtein values (0..4) in the
#'  proteinTable to words (like in Proteome Discoverer).
#' 
#' @param info integer vector to be 'translated'
#' 
#' @return character vector (the translation)
#' @export
isMasterProtein <- translateInfo(
  translation = c("Master Protein",
                  "Master Protein Candidate",
                  "<Undefined>",
                  "Master Protein Rejected",
                  "Master Protein Rejected"),
  setNAZero = FALSE, zeroIsFirst = TRUE)

#' function for 'translation' of the QuanInfoDetails values in the
#'  QuanSpectrumInfo table to words (like in Proteome Discoverer).
#' 
#' @param info integer vector to be 'translated'
#' 
#' @return character vector (the translation)
#' @export
quanInfoDetails <- translateInfo(
  translation = c("N/A",
                  "Filtered by Isolation Interference",
                  "Filtered by average S/N"),
  setNAZero = TRUE, zeroIsFirst = TRUE)

#' function for 'translation' of the QuanInfo values in the QuanSpectrumInfo
#'  table to words (like in Proteome Discoverer).
#' 
#' @param info integer vector to be 'translated'
#' 
#' @return character vector (the translation)
#' @export
quanInfo <- translateInfo(
  translation = c("N/A",
                  "No Quan Values",
                  "Missing Values",
                  "Rejected by Method",
                  "Shared"),
  setNAZero = TRUE, zeroIsFirst = TRUE)

#' function for translation of the QuanInfos values in the psms & peptide
#'  tables to words (like in Proteome Discoverer).
#' 
#' @param info integer vector to be 'translated'
#' 
#' @return character vector (the translation)
#' @export
pQuanInfo <- translateInfo(
  translation = c("Unknown",
                  "No Quan Values",
                  "Shared",
                  "Redundant",
                  "No Proteins",
                  "Filtered Out",
                  "No Quan Labels",
                  "Inconsistently Labelled",
                  "Incompatible Labels",
                  "Unique",
                  "Exluded by Method",
                  "Indistinguishable Channels",
                  "Breaks Strict Parsimony",
                  "Excluded Modification",
                  "Not Reliable",
                  "Mandatory Modification Missing"),
  setNAZero = FALSE, zeroIsFirst = TRUE)