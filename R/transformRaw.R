#' Converts a raw vector into its numeric counterpart(s)
#'
#' @param rawVector   A vector of type 'raw' (blob)
#' @param minimumSize A number (integer) indicating how many numbers (numeric)
#'                    are in the rawVector. Default = 1. This catches potential
#'                    rawVectors that are empty/NA
#' @return a list of numbers (numeric) or NA (if rawVector is empty)
#' @note numeric vectors are 9 bytes long, the first 8 are the actual
#'       number, if the last byte (9) == 0 then the result is NA
#'
convertRawNumeric <- function(rawVector, minimumSize = 1){
  lengthVector <- length(rawVector) %/% 9
  if (lengthVector == 0){
    return(rep(as.list(NA), minimumSize))
  } else {
    return(
      lapply(0:(lengthVector-1),function(x){
        if (!is.null(rawVector[((x*9)+1):((x+1)*9)])){
          if (rawVector[((x*9)+1):((x+1)*9)][[9]] != as.raw(0)){
            return(readBin(rawVector[((x*9)+1):((x+1)*9)],"numeric"))
          }
        }
        return(NA)
      })
    )
  }
}

#' Converts a raw vector into its integer counterpart(s)
#'
#' @param rawVector   A vector of type 'raw' (blob)
#' @param minimumSize A number (integer) indicating how many numbers (numeric)
#'                    are in the rawVector. Default = 1. This catches potential
#'                    rawVectors that are empty/NA
#' @return a list of numbers (integer) or NA (if rawVector is empty)
#' @note integer vectors are 5 bytes long, the first 4 are the actual
#'       number, if the last byte (5) == 0 then the result is NA
convertRawInteger <- function(rawVector, minimumSize = 1){
  lengthVector <- length(rawVector) %/% 5
  if (lengthVector == 0){
    return(rep(as.list(NA), minimumSize))
  } else {
    return(
      lapply(0:(lengthVector-1),function(x){
        if (!is.null(rawVector[((x*5)+1):((x+1)*5)])){
          if (rawVector[((x*5)+1):((x+1)*5)][[5]] != as.raw(0)){
            return(readBin(rawVector[((x*5)+1):((x+1)*5)],"integer"))
          }
        }
        return(NA)
      })
    )
  }
}

#' Converts a raw vector into its boolean counterpart(s)
#'
#' @param rawVector   A vector of type 'raw' (blob)
#' @param minimumSize A number (integer) indicating how many booleans are
#'                    in the rawVector. Default = 1. This catches potential
#'                    rawVectors that are empty/NA
#' @return a list of booleans or NA (if rawVector is empty)
#' @note boolean vectors are 5 bytes long, the first 4 are the actual
#'       boolean (anything but 0 is TRUE), if the last byte (5) == 0 then
#'       the result is NA
convertRawSpecial <- function(rawVector, size, minimumSize = 1){
  if (identical(rawVector,NULL)){
    return((rep(list(NA),minimumSize)))
  } else {
    lengthVector <- length(rawVector) %/% size
    lapply(0:(lengthVector-1),function(x){
      if (!is.null(rawVector[((x*size)+1):((x+1)*5)])){
        if (rawVector[((x*size)+1):((x+1)*size)][[size]] != as.raw(0)){
          return(ifelse(rawVector[((x*size)+1)] == as.raw(0),FALSE,TRUE))
        }
      }
      return(NA)
    })
  }
}

#' specials are not numeric or integer, but have chunks of a certain size
#'
#' @format data.frame with columns 'names' and 'size'
#' @note each chunk exists of two bytes, first one is logical (boolean):
#' zero = FALSE, otherwise TRUE. Second byte = also logical: determines if
#' value is NA (1) or not (0)
columnSpecials <- data.frame(names = c("AspectBiologicalProcess",
                                       "AspectCellularComponent",
                                       "AspectMolecularFunction"),
                             size = c(2,2,2))

#' function that converts a data.frame column with raw vectors into
#' one or more columns containing the integer/numeric/... counterparts
#'
#' @param columnVector data.frame column, eg df[,1] or df$column1
#' @param what         "integer" or "numeric" (can only be one of these two!)
#' @param columnName   character string, can be new name, can be the original
#'                     column's name. Must be suitable name
#' @param minimumSize  A number (integer) indicating how many values are
#'                     in the rawVectors that make up the columnVector.
#'                     Default = 1. This catches potential rawVectors
#'                     that are empty/NA
#' @param forceBlob    must be data.frame with 3 columns: name (columnName),
#'                     what (type)& minimumSize (number of values in a cell)
#'                     Default = NA. forceBlob exists to override automatic
#'                     conversion which has a potential for mistakes when
#'                     determining the type of a rawVector in a columnVector.
#'                     If 'what' in the forceBlob data.frame = NA, then the
#'                     columnVector will not be converted, but returned as it
#'                     is!
#' @param allowSpecials  allow for the 'special' cpnversions of columns with
#'                       the names specified om the specials data.frame.
#'                       Default = TRUE. Note: it is best to use forceBlob
#'                       with 'what' set to NA to prevent conversion problems
#' @return a data.frame (!) with one or more of the converted (or not) columns
#' @note: if a columnVector contains multiple values per 'cell' then
#'        the column will get split into an equal number of columns with
#'        the name 'columnName'+"_"+number, eg column1 becomes:
#'        column1_1, column1_2, etc
convertRawColumn <- function(columnVector, what, columnName, minimumSize = 1,
                             forceBlob = NA, allowSpecials = TRUE){
  if (allowSpecials & (columnName %in% columnSpecials$names)){
    specialSize <- columnSpecials[columnSpecials$names == columnName,]$size
    converted <- unlist(lapply(columnVector,
                               function(x){
                                 convertRawSpecial(x,specialSize,
                                                   minimumSize = minimumSize)}))
  } else {
    if (!identical(forceBlob,NA)){
      if (columnName %in% forceBlob$name){
        what <- forceBlob[forceBlob$name == columnName,]$what
        if (is.na(what)){
          converted <- data.frame(tempName = columnVector)
          colnames(converted) <- columnName
          return(columnVector) # no conversion is (to be) performed!
        }
        minimumSize <- forceBlob[forceBlob$name == columnName,]$minimumSize
      }
    }
    if (what == "integer"){
      converted <- unlist(lapply(columnVector, function(x){convertRawInteger(x, minimumSize = minimumSize)}))
    } else {
      if (what == "numeric"){
        converted <- unlist(lapply(columnVector,function(x){convertRawNumeric(x, minimumSize = minimumSize)}))
      }
    }
  }
  numberColumns <- minimumSize  #length(converted) %/% length(columnVector)
  if (numberColumns > 1){
    tempdf <- data.frame()
    for (counter in 1:numberColumns){
      mask <- rep(FALSE,numberColumns)
      mask[counter] <- TRUE
      if (nrow(tempdf) == 0){
        tempdf <- data.frame(tempName = converted[mask])
      } else {
        tempdf <- dplyr::bind_cols(tempdf, data.frame(tempName = converted[mask]))
      }
      colnames(tempdf)[counter] <- paste(c(columnName,"_",toString(counter)),collapse = "")
    }
    return(tempdf)
  } else {
    converted <- data.frame(tempName = converted)
    colnames(converted) <- columnName
    return(converted)
  }
}

#' Determine the type of raw vector a rawVector is
#'
#' @param rawVector  A vector of type 'raw' (blob)
#' @param number     boolean to specifiy if the return value should be
#'                   be character (FALSE) or numeric (TRUE). Default = TRUE
#' @return character string ("integer", "numeric") or a numeric value (5 or 9)
#' @note this method has a weakness, if raw vector length = 45 or similar,
#'       then the type returned will be integer, even though it may very well
#'       be numeric  45 %% 5 = 0. but 45 %% 9 = 0 too!
whichRaw <- function(rawVector, number = FALSE){
  if (identical(rawVector,NA) | identical(rawVector,NULL)){
    return(NA)
  }
  if (((length(rawVector) %% 5) == 0)){
    return(ifelse(number, 5, "integer"))
  } else {
    if (((length(rawVector) %% 9) == 0)){
      return(ifelse(number, 9, "numeric"))
    } else {
      return(ifelse(number, as.numeric(NA), as.character(NA)))
    }
  }
}

#' Determine the type of a list of raw vectors
#'
#' @param blobList A list of raw vectors, usually the column of a data.frame
#' @param naValue  A default value to be used, when it the type cannot be
#'                 determined (when a column contains only NA's)
#' @return a list of (type = , number = ) where type is the type of raw vector
#'         and number is the number of items per cell
#' @note Some of the data in the tables that come out of the database system
#' in # proteome discoverer are in the <blob> raw Vector format. So some
#' columns are lists of raw Vectors. To convert them, one has to know the type
#' of raw vector (integer/numeric). Because some elements in the 'blob' list
#' are NA, this function seeks the first non-NA element and tries to determine
#' it's type. If it the column contains only NA's then the naValue will be
#' returned (default = NA). This was included to force a certain type
#' even if the column (temporarily?) contains no data
#' @note only properly works with "integer" or "numeric" raw vectors
whichRawList <- function(blobList, naValue = NA){
  anyValidBlobs <- sum(!is.na(blobList))
  if (anyValidBlobs < 1){
    theType <- naValue
  } else {
    firstBlob <- blobList[[which(!is.na(blobList))[1]]]
    theType <- whichRaw(firstBlob)
  }
  numberValues <- ifelse(theType == "integer",length(firstBlob)/5,
                         ifelse(theType == "numeric", length(firstBlob)/9,
                                length(firstBlob)))
                                # unknowns will give length of the blob itself
  return(list(type = theType, number = numberValues))
}

#' Determine the type of a list of raw vectors that is a columnSpecial
#'
#' @param blobList A list of raw vectors, usually the column of a data.frame
#' @param naValue  A default value to be used, when it the type cannot be
#'                 determined (when a column contains only NA's)
#' @return a list of (type = , number = ) where type is the type of raw vector
#'         and number is the number of items per cell
#' @note the type is always either NA or "special". This function is actually
#' more for determining the number of items per cell, which is a list
#' of boolean values
whichRawListSpecial <- function(blobList, naValue = NA, specialName){
  anyValidBlobs <- sum(!is.na(blobList))
  if (anyValidBlobs < 1){
    theType = naValue
    numberValues <- 0
  } else {
    firstBlob <- blobList[[which(!is.na(blobList))[1]]]
    theType <- "special"
    numberValues <-
      length(firstBlob)/
              columnSpecials[columnSpecials$names == specialName,]$size
      # unknowns will give length of the blob itself
  }
  return(list(type = theType, number = numberValues))
}


#' automatically attempts determines the types (numeric or integer) of the
#' raw vector columns in a data.frame
#'
#' @param df             data.frame containing ONLY blob/rawVector columns
#' @param allowSpecials  allow for the 'special' conversions of columns with
#'                       the names specified om the specials data.frame.
#'                       Default = TRUE. Note: it is best to force the type
#'                       if this is set to FALSE
#' @return a list of the raw types
determineRawTypes <- function(df, allowSpecials = TRUE){
  rawType = list()
  for (counter in 1:(ncol(df))){
    blobList <- df[,counter][[1]]
    # as.list(df[,counter])[1] --> doesn't work when first item = NA ?!
    if (allowSpecials & (colnames(df[,counter]) %in% columnSpecials$names)){
      # special
      rawType[[counter]] <- whichRawListSpecial(blobList,
                                          specialName = colnames(df[,counter]))
    } else {
      rawType[[counter]] <- whichRawList(blobList)
    }
  }
  return(rawType)
}

#' df_transform_raws(): converts raw columns in a data.frame to the correct
#' data types
#'
#' @param df   data.frame coming from a table from a Proteome Discoverer
#'             database (eg .pdResult files)
#' @param forceBlob  must be data.frame with 3 columns: name (columnName),
#'                   what (type)& minimumSize (number of values in a cell)
#'                   default = NA. forceBlob exists to override automatic
#'                   conversion (which has a potential for mistakes when
#'                   determining the type of a rawVector in a columnVector).
#'                   If 'what' in the forceBlob data.frame = NA, then the
#'                   columnVector will not be converted, but returned as it
#'                   is!
#' @param allowSpecials  allow for the 'special' cpnversions of columns with
#'                       the names specified om the specials data.frame.
#'                       Default = TRUE. Note: it is best to use forceBlob
#'                       with 'what' set to NA to prevent conversion problems
#' @return data.frame with all raw vector ('blob') columns converted to more
#'         more regular data types
#' @note the tables/data.frame's coming from a Proteome Discoverer database
#' (eg .pdResult files) have columns of the type raw vecotr (blob).
#' These can be converted automatically or semi-automatically by this function
#' @note If there are no raw vector columns, then this function has no use and
#' may even trigger errors/warnings
#' @note This function can only do integer & numeric blob columns (and the
#' specials) at the moment
#' @note some raw vector columns are actually two (or possibly more) columns
#' in one. In those cases each element/cell of the column is two (or more)
#' values. This function splits these columns into two seperate ones.
#' @export
dfTransformRaws <- function(df, forceBlob = NA, allowSpecials = TRUE){
  # figure out which columns are raw vector (class 'blob')
  blobColumns <- which(
    unlist(lapply(unname(lapply(df,
                                function(x){class(x)})),
                  function(x){return("blob" %in% x)})))
  # determine the type of each blob column
  colClass <- determineRawTypes(df[blobColumns])
  # start a new data.frame w/o the blob columns
  newdf <- df[,-blobColumns]
  blobNames <- names(df[,blobColumns])
  for (counter in seq_along(blobColumns)){
    if (!identical(colClass[counter],NA)){
      newColumns <- convertRawColumn(
        df[,blobColumns[counter]][[1]],
        what = colClass[[counter]]$type,
        columnName = blobNames[counter],
        minimumSize = colClass[[counter]]$number,
        forceBlob = forceBlob
        )
      newdf <- dplyr::bind_cols(newdf, newColumns)
    } else {
      warningMessage <- paste(
        c("Warning: cannot automatically convert column '",
          blobNames[counter],
          "' "),
        collapse = "")
      warning(warningMessage)
      newdf <- dplyr::bind_cols(newdf, df[,blobColumns[counter]])
    }
  }
  return(newdf)
}
