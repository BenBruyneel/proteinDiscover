#' Converts a raw vector into its numeric counterpart(s)
#'
#' @param rawVector   A vector of type 'raw' (blob)
#' @param minimumSize A number (integer) indicating how many numbers (numeric)
#'                    are in the rawVector. Default = 1. This catches potential
#'                    rawVectors that are empty/NA
#' @return a list of numbers (numeric) or NA (if rawVector is empty)
#' 
#' @note numeric vectors are 9 bytes long, the first 8 are the actual
#'       number, if the last byte (9) == 0 then the result is NA
#' @note internal function
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
#' 
#' @note integer vectors are 5 bytes long, the first 4 are the actual
#'       number, if the last byte (5) == 0 then the result is NA
#' @note internal function
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
#' @param specialSize length of the result, ie the number of columns in the raw.
#'                    This is relevant for the so called specials. It indicates
#'                    how many booleans are in the rawVector
#' @return a list of booleans or NA (if rawVector is empty)
#' 
#' @note boolean vectors are 5 bytes long, the first 4 are the actual
#'       boolean (anything but 0 is TRUE), if the last byte (5) == 0 then
#'       the result is NA
#' @note internal function
convertRawSpecial <- function(rawVector, specialSize){
  if (identical(rawVector,NULL)){
    return((rep(list(NA), specialSize)))
  } else {
    lengthVector <- specialSize
    lapply(0:(lengthVector-1),function(x){
      if (!is.null(rawVector[((x*2)+1):((x+1)*5)])){
        if (rawVector[((x*2)+1):((x+1)*2)][[2]] != as.raw(0)){
          return(ifelse(rawVector[((x*2)+1)] == as.raw(0),FALSE,TRUE))
        }
      }
      return(NA)
    })
  }
}

#' Specials are not numeric or integer, but have chunks of a certain size
#' All encountered in Proteome Discoverer are actually booleans with a value
#' 0 (FALSE), 1 (TRUE) or NA
#' 
#' @return data.frame with columns 'names' and 'size'
#' @note each chunk consists of two bytes, first one is logical (boolean):
#' zero = FALSE, otherwise TRUE. Second byte = also logical: determines if
#' value is NA (1) or not (0)
#' @export
columnSpecials <- function(){
  return(data.frame(names = c("AspectBiologicalProcess",
                              "AspectCellularComponent",
                              "AspectMolecularFunction"),
                    size = c(2,2,2)))
}

#' function that converts a data.frame column with raw vectors into
#' one or more columns containing the integer/numeric/... counterparts
#'
#' @param columnVector data.frame column, eg df[,1] or df$column1
#' @param blobDF       must be data.frame with 3 columns: name (columnName),
#'                     what (type) & minimumSize (number of values in a cell)
#'                     Default = NA. forceBlob exists to override automatic
#'                     conversion which has a potential for mistakes when
#'                     determining the type of a rawVector in a columnVector.
#'                     If 'what' in the forceBlob data.frame = NA, then the
#'                     columnVector will not be converted, but returned as it
#'                     is!
#' @return a data.frame (!) with one or more of the converted (or not) columns
#' 
#' @note: if a columnVector contains multiple values per 'cell' then
#'        the column will get split into an equal number of columns with
#'        the name 'columnName'+"_"+number, eg column1 becomes:
#'        column1_1, column1_2, etc
#' @note internal function
convertRawColumn <- function(columnVector, blobDF){
  if (colnames(columnVector)[1] %in% blobDF$name){
    blobDF <- blobDF[blobDF$name == colnames(columnVector)[1],]
    if (is.na(blobDF$what[1])){
      return(columnVector) # no conversion is (to be) performed!
     } 
  }
  if (blobDF$what[1] == "integer"){
    converted <- unlist(lapply(columnVector[,1], function(x){convertRawInteger(x, minimumSize = blobDF$minimumSize[1])}))
  } else {
    if (blobDF$what[1] == "numeric"){
      converted <- unlist(lapply(columnVector[,1],function(x){convertRawNumeric(x, minimumSize = blobDF$minimumSize[1])}))
    } else {
      if (blobDF$what[1] == "special"){
        converted <- unlist(lapply(columnVector[,1],
                                   function(x){
                                     convertRawSpecial(x,
                                         specialSize = blobDF$minimumSize[1])}))
      }
    }
  }
  numberColumns <- blobDF$minimumSize[1]
  # numberColumns <- minimumSize  #length(converted) %/% length(columnVector)
  tempdf <- data.frame(matrix(converted, ncol = numberColumns, byrow = TRUE))
  colnames(tempdf) <- paste(blobDF$name[1],"_",1:numberColumns,sep = "")
  return(tempdf)
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
      colnames(tempdf)[counter] <- paste(c(blobDF$name[1],"_",toString(counter)),collapse = "")
    }
    return(tempdf)
  } else {
    converted <- data.frame(tempName = converted)
    colnames(converted) <- blobDF$name[1]
    return(converted)
  }
}

#' detemines which columns in a table are of the blob (raw) type
#' 
#' @param theTable the table containing the data
#' @return a data.frame with two columns: name = colum name) and type (which
#'  should always be 'blob')
#'  
#' @note meant for use in debugging problems
#' @export
getBlobs <- function(theTable){
  temp <- purrr::map_chr(theTable,~class(.x)[1])
  return(data.frame(name = names(temp), type = unname(temp)) %>%
           dplyr::filter(.data$type == "blob"))
}

#' attempts to determine the length (in bytes) of the individual elements of a
#'  blob-type column of a data.frame. It should (!) return an integer value of
#'  course (as all elements are supposed to have the same length). Also: if all
#'  elements of the column are NA, the the result will be NaN
#'  
#' @param blobList one column of a data.frame (as a list) of blob (raw) element
#'  type elements
#' @return the length of the elements in the data.frame (or list) column. Again:
#'  this should be an integer
#'  
#' @note meant for use in debugging problems
#' @export
blobLength <- function(blobList){
  return(mean(unlist(lapply((lapply(blobList,length)),
                            function(x){ifelse(x==0,NA,x)})),
              na.rm = TRUE))
}

#' determines the length of the blob (raw) type columns in a data.frame from a
#'  table (possibly from a database) 
#'  
#' @param blobDF the result from getBlobs
#' @param theTable the table to which blobDF refers (the table (data.frame) used
#'  for getBlobs)
#' @return blobDF with a single column (named 'length') added which contains
#'  the length (number of bytes) of each blob column
#'  
#' @note internal function 
determineBlobLengths <- function(blobDF, theTable){
  blobDF$length <- unlist(lapply(1:nrow(blobDF),
                        function(x){blobLength(theTable[,blobDF$name[x]])}))
  return(blobDF)
}

#' function which converts the blob length to the appropiate (R) datatype
#' 
#' @param blobLength currently only 5 (integer) & 9 (numeric) are supported.
#'  Anything else will result in NA
#' @return a string (either "integer" or "numeric") or NA
#' 
#' @note internal function
determineBlobTypeRaw <- function(blobLength){
  return(unlist(lapply(blobLength, function(x){
                                              if (x %% 9 == 0){
                                                return("numeric")
                                              } else {
                                                if (x %% 5 == 0){
                                                  return("integer")
                                                } else {
                                                  return(NA)
                                                }
                                              }
                                            })))
}

#' function that attempts to assign a type to the blob (raw) length  as found
#'  by dtermineBlobLengths
#'  
#' @note this function works only with single numbers and is meant to be used
#'  primarily by the function blobEstimateTypes
#'  
#' @param blobLength the actual length (number of bytes) of the element we wish
#'  to assign a type to
#' @param minimumNumber this defines the minimum number of columns a
#'  blob/raw type column should be split into. In TMT10plex experiments, the
#'  minimumNumber will usually be 10, becauseyou have 10 channels/abundances
#' @param numberOfGroups this defines how many 'groups' are present in the data.
#'  Taking Abundances as an example: Proteone Discoverer has both the original
#'  columns (say Abundances_1 through Abundances_2), but also columns where the
#'  abundances, that 'belong' together, are eg averaged or some other
#'  (statistical) measure  is calculated over a number of columns. You may have
#'  eg 10 'Abundance channels' which are 5 samples total, each in duplo. This
#'  means that some columns in the resulting table will need to be split in 10
#'  different columns (the original 'Abundances') while 'grouped' columns should
#'  be split into 5 different columns (eg the calculated means or variations of
#'  the 'abundances' columns). Note that although not enforced by the code, the
#'  numberOfGroups should always be equal or less than the  minimumNumber
#'  parameter. Default value = minimumNumber
#' @param ratioNumberOfGroups when ratios between groups are calculated we get
#'  columns (ratio columns) that need to be split into numberOfGroups - 1
#'  (which is the efault value)
#' @return a single row data.frame with columns what (type) and minimumSize
#'  (number of variables in the blob)
#'  
#' @note this function does not deal properly with specials, their types/
#'  translations are resolved in a different way
#' @note there are two ways to see potential problems with the type assignments:
#'  the columns may contain NA values
#' @note internal function
determineBlobType <- function(blobLength, minimumNumber = 1,
                              numberOfGroups = minimumNumber,
                              ratioNumberOfGroups = numberOfGroups - 1){
  if (blobLength %in% c(5,9)){
    return(data.frame(what = determineBlobTypeRaw(blobLength), minimumSize = 1))
  }
  if ((blobLength %% minimumNumber) == 0){
    if (minimumNumber == 1){
      result <- data.frame(what = NA,
                           minimumSize = NA)
    } else {
      result <- data.frame(
        what = determineBlobTypeRaw(blobLength = blobLength %/% minimumNumber),
        minimumSize = minimumNumber)
    }
  } else {
    if ((blobLength %% numberOfGroups) == 0){
      if (numberOfGroups == 1){
        result <- data.frame(what = NA,
                             minimumSize = NA)
      } else {
        result <- data.frame(
        what = determineBlobTypeRaw(blobLength = blobLength %/% numberOfGroups),
        minimumSize = numberOfGroups)
      }
    } else {
      if ((blobLength %% ratioNumberOfGroups) == 0){
        if (ratioNumberOfGroups == 1){
          result <- data.frame(what = NA,
                               minimumSize = NA)
        } else {
          result <- data.frame(
  what = determineBlobTypeRaw(blobLength = blobLength %/% ratioNumberOfGroups),
  minimumSize = ratioNumberOfGroups)
        }
      } else {
        result <- data.frame(what = NA, minimumSize = NA)
      }
    }
  }
  if (identical(result$what[1],NA)){
    # this section is problematic with eg blobLengths of eg 45,
    # since this is integer divisble by both 5 & 9
    if (blobLength %% 9 == 0){
      return(data.frame(what = "numeric",
                        minimumSize = blobLength %/% 9))
    } else {
      if (blobLength %% 5 == 0){
        return(data.frame(what = "integer",
                          minimumSize = blobLength %/% 5))
      } else {
        return(data.frame(what = NA, minimumSize = NA))
      }
    }
  }
  return(result)
}

#' function that attempts to assign a type to the blob (raw) lengths as found
#'  by dtermineBlobLengths
#'  
#' @note this function works with single numbers and multiple numbers
#'  
#' @param blobLengths the actual lengths (number of bytes) of the elements we
#'  wish to assign types to. Can be 1 or more lengths
#' @param minimumNumber this defines the minimum number of columns a
#'  blob/raw type column should be split into. In TMT10plex experiments, the
#'  minimumNumber will usually be 10, becauseyou have 10 channels/abundances
#' @param numberOfGroups this defines how many 'groups' are present in the data.
#'  Taking Abundances as an example: Proteone Discoverer has both the original
#'  columns (say Abundances_1 through Abundances_2), but also columns where the
#'  abundances, that 'belong' together, are eg averaged or some other
#'  (statistical) measure  is calculated over a number of columns. You may have
#'  eg 10 'Abundance channels' which are 5 samples total, each in duplo. This
#'  means that some columns in the resulting table will need to be split in 10
#'  different columns (the original 'Abundances') while 'grouped' columns should
#'  be split into 5 different columns (eg the calculated means or variations of
#'  the 'abundances' columns). Note that although not enforced by the code, the
#'  numberOfGroups should always be equal or less than the  minimumNumber
#'  parameter. Default value = minimumNumber
#' @param ratioNumberOfGroups when ratios between groups are calculated we get
#'  columns (ratio columns) that need to be split into numberOfGroups - 1
#'  (which is the efault value)
#' @return a data.frame with columns what (type) and minimumSize (number of
#'  variables in the blob)  
#'  
#' @note this function does not deal properly with specials, their types/
#'  translations are resolved in a different way
#' @note there are two ways to see potential problems with the type assignments:
#'  the columns may contain NA values
#' @note internal function 
blobEstimateTypes <- function(blobLengths, minimumNumber,
                              numberOfGroups = minimumNumber,
                              ratioNumberOfGroups = numberOfGroups - 1){
  return(
    dplyr::bind_rows(
      lapply(blobLengths, function(x){determineBlobType(blobLength = x,
                                                        numberOfGroups = numberOfGroups,
                                                        minimumNumber = minimumNumber,
                                                        ratioNumberOfGroups = ratioNumberOfGroups)})
    )
  )
}

#' function that attempts to assign types and sizes to the blob type columns
#'  in a table. The result from this function can be used in the dfTransformRaws 
#'  function
#' 
#' @param theTable a data.frame with blob Columns (if no blobColumns are
#'  present, then NA is returned)
#' @param minimumNumber this defines the minimum number of columns a
#'  blob/raw type column should be split into. In TMT10plex experiments, the
#'  minimumNumber will usually be 10, becauseyou have 10 channels/abundances
#' @param numberOfGroups this defines how many 'groups' are present in the data.
#'  Taking Abundances as an example: Proteone Discoverer has both the original
#'  columns (say Abundances_1 through Abundances_2), but also columns where the
#'  abundances, that 'belong' together, are eg averaged or some other
#'  (statistical) measure  is calculated over a number of columns. You may have
#'  eg 10 'Abundance channels' which are 5 samples total, each in duplo. This
#'  means that some columns in the resulting table will need to be split in 10
#'  different columns (the original 'Abundances') while 'grouped' columns should
#'  be split into 5 different columns (eg the calculated means or variations of
#'  the 'abundances' columns). Note that although not enforced by the code, the
#'  numberOfGroups should always be equal or less than the  minimumNumber
#'  parameter. Default value = minimumNumber
#' @param ratioNumberOfGroups when ratios between groups are calculated we get
#'  columns (ratio columns) that need to be split into numberOfGroups - 1
#'  (which is the efault value)
#' @param blobDF essentially the result from either getBlobs; if NA then it will
#'  be generated by the getBlobs function with theTable as an argument
#' @param specials default is TRUE, means that specials will be taken care of
#' @return a data.frame with the name of the blob columns, their lengths,
#'  what (type) and minimumSize (number of variables in the blob)
#'  
#' @note this function does not deal properly with specials, their types/
#'  translations are resolved in a different way
#' @note there are two ways to see potential problems with the type assignments:
#'  the columns may contain NA values
#' @export
determineBlobTypes <- function(theTable, minimumNumber = 1,
                               numberOfGroups = minimumNumber,
                               ratioNumberOfGroups = numberOfGroups - 1,
                               blobDF = NA, specials = TRUE){
  blobDF <- getBlobs(theTable = theTable)
  if (nrow(blobDF) == 0){
    return(NA) # no blobs
  }
  if (!specials){
   blobDF <- blobDF %>% dplyr::filter(!(.data$name %in% columnSpecials()$names))
  }
  blobDF <- determineBlobLengths(blobDF = blobDF, theTable = theTable)
  blobDF <- dplyr::bind_cols(blobDF, blobEstimateTypes(
                                                blobLengths = blobDF$length,
                                                minimumNumber = minimumNumber,
                                                numberOfGroups = numberOfGroups,
                                     ratioNumberOfGroups = ratioNumberOfGroups))
  if (specials){
    blobSpecials <- columnSpecials()
    if (sum(grepl(blobDF$name, pattern = "Aspect")) > 0){
      for (counter in 1:nrow(blobDF)){
        if (blobDF$name[counter] %in% blobSpecials$names){
          blobDF$what[counter] <- "special"
          blobDF$minimumSize[counter] <- 
            blobDF$length[counter] %/%
            (blobSpecials %>%
               dplyr::filter(names == blobDF$name[counter]))$size
        }
      }
    }
  }
  return(blobDF)
}

#' df_transform_raws(): converts raw columns in a data.frame to the correct
#' data types
#'
#' @param df   data.frame coming from a table from a Proteome Discoverer
#'  database (eg .pdResult files)
#' @param blobDF  must be data.frame with 3 columns: name (columnName),
#'  what (type) & minimumSize (number of values in a cell) default = NA.
#'  If 'what' in the data.frame = NA, then the columnVector will not be
#'  converted, but returned as it is
#' @param minimumNumber this defines the minimum number of columns a
#'  blob/raw type column should be split into. In TMT10plex experiments, the
#'  minimumNumber will usually be 10, becauseyou have 10 channels/abundances
#' @param numberOfGroups this defines how many 'groups' are present in the data.
#'  Taking Abundances as an example: Proteone Discoverer has both the original
#'  columns (say Abundances_1 through Abundances_2), but also columns where the
#'  abundances, that 'belong' together, are eg averaged or some other
#'  (statistical) measure  is calculated over a number of columns. You may have
#'  eg 10 'Abundance channels' which are 5 samples total, each in duplo. This
#'  means that some columns in the resulting table will need to be split in 10
#'  different columns (the original 'Abundances') while 'grouped' columns should
#'  be split into 5 different columns (eg the calculated means or variations of
#'  the 'abundances' columns). Note that although not enforced by the code, the
#'  numberOfGroups should always be equal or less than the  minimumNumber
#'  parameter. Default value = minimumNumber
#' @param ratioNumberOfGroups when ratios between groups are calculated we get
#'  columns (ratio columns) that need to be split into numberOfGroups - 1
#'  (which is the efault value)
#' @param specials default is TRUE, means that specials will be taken care of
#' 
#' @return data.frame with all raw vector ('blob') columns converted to more
#'         more regular data types
#'         
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
dfTransformRaws <- function(df, blobDF = NA,
                            minimumNumber = 1,
                            numberOfGroups = minimumNumber,
                            ratioNumberOfGroups = numberOfGroups-1,
                            specials = TRUE){
  # determine the type of each blob column
  if (identical(blobDF,NA)){
    blobDF <- determineBlobTypes(theTable =  df,
                                 minimumNumber = minimumNumber,
                                 numberOfGroups = numberOfGroups,
                                 ratioNumberOfGroups = ratioNumberOfGroups,
                                 specials = specials)
  }
  if (identical(blobDF,NA)){
    return(df)
  }
  # start a new data.frame w/o the blob columns
  newdf <- df %>% dplyr::select(-blobDF$name)
  for (counter in 1:nrow(blobDF)){
    if (!identical(blobDF$type,NA)){
      newColumns <- convertRawColumn(
        df %>% dplyr::select(blobDF$name[counter]),
        blobDF = blobDF)
      newdf <- dplyr::bind_cols(newdf, newColumns)
    } else {
      warningMessage <- paste(
        c("Warning: cannot (automatically) convert column '",
          blobDF$name[counter],
          "' "),
        collapse = "")
      warning(warningMessage)
    newdf <- dplyr::bind_cols(newdf, df %>% dplyr::select(blobDF$name[counter]))
    }
  }
  return(newdf)
}