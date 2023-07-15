#' internal helper function to determine if object == whichClass or a descendant of
#'  whichClass
#' 
#' @param object a data object of some class
#' @param whichClass character string: class name to be tested
#' @return TRUE or FALSE
#' @noRd
is.Class <- function(object, whichClass){
  return(whichClass %in% class(object))
}

#' ifelse replacement for properly returning all datatypes.
#' internal helper function
#'
#' @param logicValue variable or expression resulting in TRUE or FALSE,
#'  if missing or not logical then the function will return NULL.
#' @param ifTrue variable or expression to be returned when logicValue == TRUE
#' @param ifFalse variable or expression to be returned when logicValue == FALSE
#'
#' @returns depending on logicValue, ifTrue ir ifFalse.
#' @note not vectorized
#' @noRd
ifelseProper <- function(logicValue = NULL, ifTrue = NULL, ifFalse = NULL){
  if (missing(logicValue)){
    return(NULL)
  } else {
    if (!is.logical(logicValue)){
      return(NULL)
    } else {
      if (logicValue){
        return(ifTrue)
      } else {
        return(ifFalse)
      }
    }
  }
}

#' converts character string date into date/time format
#' 
#' @param theDate character string to be converted
#' (can be vectorized)
#' @param dateFormat function that defines the output date/time format, default
#'  is lubridate::mdy_hms
#' @returns date in mdy hms format
#' @export
thermo.date <- function(theDate, dateFormat = lubridate::mdy_hms){
  return(dateFormat(theDate))
}

#' converts character string date into date/time format
#' 
#' @param theDate character string to be converted
#' (can be vectorized)
#' @param dateFormat function that defines the output date/time format, default
#'  is lubridate::ymd_hms
#' @returns date 
#' @export
system.date <- function(theDate, dateFormat = lubridate::ymd_hms){
  return(dateFormat(theDate))
}

#' fake converter for times when no conversion is wanted/needed
#' 
#' @param theDate character string
#' (can be vectorized)
#' @returns theDate (original character string)
#' @export
na.date <- function(theDate){
  return(theDate)
}

#' internal helper function to prevent having to remember the somewhat
#' long names of the most used tables
#' 
#' @param whichTable can be either "proteins","peptides","psms" or "consensus"
#'  character do not need to be lower or upper case (all are converted to upper
#'  case). If another string is used as a parameter, the function will return
#'  NA
#' @return a string containing the protein discoverer table name corresponding
#'  to the parameter whichTable
#' @export
tableNames <- function(whichTable = "proteins"){
  return(switch(toupper(toString(whichTable)),
                "PROTEINS"  = "TargetProteins",
                "PEPTIDES"  = "TargetPeptideGroups",
                "PSMS"      = "TargetPsms",
                "CONSENSUS" = "ConsensusFeatures",
                NA))
}