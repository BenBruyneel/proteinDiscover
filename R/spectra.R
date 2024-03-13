#' transforms a spectrum from the table 'MassSpectrumItems' into a R compatible
#'  list
#'
#' @param spectrumObject must be of class 'raw'
#'
#' @return a list object containing info on the spectrum (object). This list
#'  object can be further translated via the function 'translateSpectrumInfo'
#' 
#' @note this functions writes a temporary file tp disk, which is unzipped, read
#'  and deleted again
#' 
#' @export
transformSpectrumRaw <- function(spectrumObject){
  if (is.Class(spectrumObject, "raw")){
    currentPath <- getwd()
    tempPath <- tempdir()
    setwd(tempPath)
    writeBin(spectrumObject, "tempfile.zip")
    spectrum <- unzip("tempfile.zip", list = TRUE)
    if (nrow(spectrum) > 1){
      warning("Raw spectrum object contains more than one spectrum! Only the first spectrum is extracted.")
      spectrum <- spectrum[1,]
    } else {
      unzip("tempfile.zip", files = spectrum$name, overwrite = T)
      result <- XML::xmlToList(spectrum$Name)
    }
    setwd(currentPath)
    unlink(paste0(tempPath,"/tempfile.zip"))
    unlink(paste(c(tempPath, "/", spectrum$Name), collapse = ""))
    return(result)
  }
  warning("Spectrum object must be class 'raw'")
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum header
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return data.frame
#' 
#' @export
spectrum.header <- function(spectrum){
  if ("Header" %in% names(spectrum)){
    result <- as.data.frame(spectrum$Header[-which(names(spectrum$Header) == "SpectrumIdentifiers")])
    result2 <- t(as.data.frame(spectrum$Header[which(names(spectrum$Header) == "SpectrumIdentifiers")]))
    rownames(result2) <- NULL
    return(cbind(result, result2))
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum scan event
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return data.frame
#' 
#' @export
spectrum.scanEvent <- function(spectrum){
  if ("ScanEvent" %in% names(spectrum)){
    result <- as.data.frame(spectrum$ScanEvent)
    rownames(result) <- NULL
    return(result)
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum centroided spectrum
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return data.frame
#' 
#' @export
spectrum.centroid <- function(spectrum){
  if ("PeakCentroids" %in% names(spectrum)){
    result <- as.data.frame(purrr::map_df(spectrum$PeakCentroids, ~.x))
    result$X <- as.numeric(result$X)
    result$Y <- as.numeric(result$Y)
    result$Z <- as.integer(result$Z)
    result$R <- as.integer(result$R)
    result$SN <- as.numeric(result$SN)
    colnames(result) <- c("mz","intensity","charge", "resolution","signalToNoise")
    return(result)
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum profile spectrum
#' 
#' @note this is a convenience function, type of data not observed
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return NULL or NA
#' 
#' @export
spectrum.profile <- function(spectrum){
  if ("ProfilePoints" %in% names(spectrum)){
    warning("Profile Not implemented, Function returns object")
    result <- spectrum$ProfilePoints
    return(result)
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum parent header
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return data.frame
#' 
#' @export
spectrum.precursor.header <- function(spectrum){
  if ("PrecursorInfo" %in% names(spectrum)){
    spectrum <- spectrum$PrecursorInfo
    if ("SpectrumHeader" %in% names(spectrum)){
      result <- as.data.frame(spectrum$SpectrumHeader[-which(names(spectrum$SpectrumHeader) == "SpectrumIdentifiers")])
      result2 <- t(as.data.frame(spectrum$SpectrumHeader[which(names(spectrum$SpectrumHeader) == "SpectrumIdentifiers")]))
      rownames(result2) <- NULL
      return(cbind(result, result2))
    }
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum parent scan event
#'
#' @param spectrum list object containing info on a spectrum
#' @param returnRaw logical vector, if TRUE them the data is returned as a list. if FALSE (default)
#'  then a data.frame of all character-type data is returned' 
#'
#' @return data.frame or list
#' 
#' @export
spectrum.precursor.scanEvent <- function(spectrum, returnRaw = FALSE){
  if ("PrecursorInfo" %in% names(spectrum)){
    spectrum <- spectrum$PrecursorInfo
    if ("ScanEvent" %in% names(spectrum)){
      if (returnRaw) {
        return(spectrum$ScanEvent)
      } else {
        result <- as.data.frame(    # remove all non character items
          spectrum$ScanEvent[-which(unname(unlist(lapply(spectrum$ScanEvent,
                                                         class))) != "character")]
        )
        rownames(result) <- NULL
        return(result)
      }
    }
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum parent monoisotopic peak
#'
#' @param spectrum list object containing info on a spectrum
#' @param measured logical vector, if TRUE then the measured data is returned
#'
#' @return data.frame
#' 
#' @export
spectrum.precursor.info  <- function(spectrum, measured = TRUE){
  if ("PrecursorInfo" %in% names(spectrum)){
    spectrum <- spectrum$PrecursorInfo
    result <- data.frame()
    if (!measured){
      if ("MonoisotopicPeakCentroids" %in% names(spectrum)){
        result <- as.data.frame(t(spectrum$MonoisotopicPeakCentroids[[1]]))
      }
    } else {
      if ("MeasuredMonoisotopicPeakCentroids" %in% names(spectrum)){
        result <- as.data.frame(t(spectrum$MeasuredMonoisotopicPeakCentroids[[1]]))
      }
    }
    if (nrow(result) < 1){
      result <- NA
    } else {
      result$X <- as.numeric(result$X)
      result$Y <- as.numeric(result$Y)
      result$Z <- as.integer(result$Z)
      result$R <- as.integer(result$R)
      result$SN <- as.numeric(result$SN)
      colnames(result) <- c("mz","intensity","charge", "resolution","signalToNoise")
    }
    return(result)
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum parent centroided spectrum
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return data.frame
#' 
#' @export
spectrum.precursor.centroid <- function(spectrum){
  if ("PrecursorInfo" %in% names(spectrum)){
    spectrum <- spectrum$PrecursorInfo
    if ("IsotopeClusterPeakCentroids" %in% names(spectrum)){
      result <- as.data.frame(purrr::map_df(spectrum$IsotopeClusterPeakCentroids, ~.x))
      result$X <- as.numeric(result$X)
      result$Y <- as.numeric(result$Y)
      result$Z <- as.integer(result$Z)
      result$R <- as.integer(result$R)
      result$SN <- as.numeric(result$SN)
      colnames(result) <- c("mz","intensity","charge", "resolution","signalToNoise")
      return(result)
    }
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum parent profile spectrum
#' 
#' @note this is a convenience function, type of data not observed
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return NULL or NA
#' 
#' @export
spectrum.precursor.profile <- function(spectrum){
  if ("PrecursorInfo" %in% names(spectrum)){
    spectrum <- spectrum$PrecursorInfo
    if ("IsotopeClusterProfilePoints" %in% names(spectrum)){
      warning("Profile Not implemented, Function returns object")
      result <- spectrum$IsotopeClusterProfilePoints
      return(result)
    }
  }
  return(NA)
}

#' gets the info in the list object coming from the function
#' 'transformSpectrumRaw': spectrum parent additonal info
#'
#' @param spectrum list object containing info on a spectrum
#'
#' @return data.frame
#' 
#' @export
spectrum.precursor.additionalInfo <- function(spectrum){
  if ("PrecursorInfo" %in% names(spectrum)){
    spectrum <- spectrum$PrecursorInfo
    if (".attrs" %in% names(spectrum)){
      result <- as.data.frame(t(spectrum[[".attrs"]]))
      return(result)
    }
  }
  return(NA)
}