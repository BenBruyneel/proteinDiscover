library(DBI)
library(RSQLite)
library(pool)
library(dplyr)
library(stringr)
library(purrr)

# ---- raw vector (blob) related functions ----

#' Converts a raw vector into its numeric counterpart(s)
#'
#' @param rawVector   A vector of type 'raw' (blob)
#' @param minimumSize A number (integer) indicating how many numbers (numeric)
#'                    are in the rawVector. Default = 1. This catches potential
#'                    rawVectors that are empty/NA
#' @return a list of numbers (numeric) or NA (if rawVector is empty)
#' @note numeric vectors are 9 bytes long, the first 8 are the actual
#'       number, if the last byte (9) == 0 then the result is NA
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

#' specials are not numeric or integer, but have chunks of size
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
    specialSize <- columnSpecials[columnSecials$names == columnName,]$size
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
        tempdf <- bind_cols(tempdf, data.frame(tempName = converted[mask]))
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
      length(firstBlob) / specials[specials$names == specialName,]$size
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
    if (allowSpecials & (colnames(df[,counter]) %in% specials$names)){ # special
      rawType[[counter]] <- whichRawListSpecial(blobList,
                                          specialName = colnames(df[,counter]))
    } else {
      rawType[[counter]] <- whichRawList(blobList)
    }
  }
  return(rawType)
}

#' converts raw columns in a data.frame to the correct data types
#' the tables/data.frame's coming from a Proteome Discoverer database
#' (eg .pdResult files) have columns of the type raw vecotr (blob).
#' These can be converted automatically or semi-automatically by this function
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
#' @note If there are no raw vector columns, then this function has no use and
#' may even trigger errors/warnings
#' @note This function can only do integer & numeric blob columns (and the
#' specials) at the moment
#' @note some raw vector columns are actually two (or possibly more) columns
#' in one. In those cases each element/cell of the column is two (or more)
#' values. This function splits these columns into two seperate ones.
#'
#' @export
df_transform_raws <- function(df, forceBlob = NA, allowSpecials = TRUE){
  # figure out which columns are raw vector (class 'blob')
  blobColumns <- which(
    unlist(
      lapply(
        unname(
          lapply(df,
                 function(x){class(x)})),function(x){return("blob" %in% x)})))
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
      newdf <- bind_cols(newdf, newColumns)
    } else {
      warningMessage <- paste(
        c("Warning: cannot automatically convert column '",
          blobNames[counter],
          "' "),
        collapse = "")
      warning(warningMessage)
      newdf <- bind_cols(newdf, df[,blobColumns[counter]])
    }
  }
  return(newdf)
}


# ---- database related functions ----

# Note: proteome discoverer .pdResult files are in the SQLite format
# functions in this section deal with the database. Some of the
# functions work with database (query) commands, which is the preferred
# way. A lot of these functions are merely wrappers to make programming
# a little easier
#
# Note: for future developments when buidling a shiny app it's recommended
# to use the 'pool' & 'dbplyr' packages to streamline things

#' Wrapper around pool::dbPool(): opens a databse
#'
#' @param  fileName  a character vector specifying the name and location
#'                   of the database
#' @param  drv defines database connection type, default = RSQLite::SQLite()
#' @param  ... to pass on additional parameters to pool::dbPool, exmples are
#'             host = "shiny-demo.csa7qlmguqrf.us-east-1.rds.amazonaws.com"
#'             username = "guest"
#'             password = "guest"
#' @return database access 'handle'
#' @note if no file with the name 'fileName' exists, then it will be created
#' (but obviously it will be empty, so most further commands will fail)
#' @note if fileName == ":memory:" the database will be an in-memory database
#' @export
db_open <- function(fileName, drv = RSQLite::SQLite(), ...){
  return(pool::dbPool(drv = drv,
                      dbname = fileName,
                      ...))
}

#' Wrapper around pool::pooClose():  closes an open database
#' (normally opened earlier via eg db_open())
#'
#'
#' @param db   database access 'handle' to be closed
#'
#' @export
db_close <- function(db){
  pool::poolClose(db)
}

#' find the number of rows in a table of an (open) database
#'
#' @param db    database access 'handle'
#' @param tabl  name of the table
#' @return number of rows in the table
#'
#' @export
db_nrow <- function(db,tabl){
  return(pool::dbGetQuery(db,paste(c("SELECT COUNT(*) FROM ",tabl),
                                   collapse = ""))[1,1])
}

#' get the classes of the different columns in a database table
#'
#' @param db    database access 'handle'
#' @param tabl  name of the table
#' @return a data.frame with columns name & type (= class)
#'
#' @export
db_columnInfo <- function(db, tabl){
  # get first row of the table to determine names & classes
  temp <- pool::dbGetQuery(db1, "SELECT * FROM TargetProteins LIMIT 1")
  # determine classes, note that in case of class lists
  # (as with type OOP descendants) this statement only take the first name
  # of the list
  temp <- map_chr(temp,~class(.x)[1])
  return(data.frame(name = names(temp), type = unname(temp)))
}

#' gathers info on all tables in the database
#'
#' @param db    database access 'handle'
#' @return a data.frame with the following columns: name, nrows, info
#' every entry in the column info contains a data.frame with the names
#' of the columns in the table and their type (see also db_columnInfo)
#' @export
db_tbl_def <- function(db){
  datazz <- pool::dbListTables(db)
  dftbl <- data.frame(name = datazz,
                      nrows = unlist(lapply(datazz,
                                            function(x){db_nrow(db,x)})),
                      info = NA)
  for (counter in 1: length(datazz)){
    dftbl$info[counter] <- list(db_columnInfo(db,datazz[counter]))
  }
  return(dftbl)
}

