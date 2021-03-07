#' get a table from a .pdResult file
#'
#' @param db database access 'handle'
#' @param tableName used to pass on the name of the table containing the data
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @filtering allows for " WHERE <expression>" additions to the SQL statement
#'  default = " " (no filtering). Note: always put a space (" ") before any
#'  statement
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid value is a character character
#'  vector of columnNames to be used for sorting string (with "ASC" or "DESC"
#'  if needed) 
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return a data.frame containing requested data from a database table or
#'  a character string specifying a SQL query
#' @export
dbGetTable <- function(db,
                       tableName,
                       columnNames = NA, 
                       filtering = " ",
                       sortOrder = NA,
                       SQL = FALSE){
  if (!identical(sortOrder,NA)){
    sortColumns = paste(c(sortOrder),collapse = ", ")
  } else {
    sortOrder = NA
  }
  query <- paste(c("SELECT ",
                   ifelse(identical(columnNames,NA),
                          "*",
                          paste(columnNames,
                                collapse = ", ")),
                   " FROM ", tableName,
                   filtering, 
                   ifelse(identical(sortOrder,NA),
                          "",
                          paste(c(" ORDER BY ",sortColumns),
                                collapse = "")
                   )
  ), collapse = "")
  if (!SQL){
    return(pool::dbGetQuery(db, query))
  } else (
    return(query)
  )
}

#' get the protein table from a .pdResult file
#' (essentially a wrapper around db_getTable())
#'
#' @param db database access 'handle'
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param masterProtein use the IsMasterProtein column to be zero,
#'  default == TRUE. If more advanced filtering is needed, use db_getTable()
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid values are a single character
#'  string ("ASC" or "DESC") or a character vector of the same length as the
#'  columnNames vector containing a series of "ASC" or "DESC"
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return a data.frame containing requested data from the protein table or
#'  a character string specifying a SQL query
#' @export
dbGetProteinTable <- function(db,
                               columnNames = NA,
                               masterProtein = TRUE,
                               sortOrder = NA,
                               SQL = FALSE){
  return(dbGetTable(
    db = db,
    tableName = "TargetProteins",
    columnNames = columnNames,
    filtering = ifelse(masterProtein,
                       " WHERE IsMasterProtein = 0",""),
    sortOrder = sortOrder,
    SQL = SQL))
}

#' get the peptideID's from (a set of) proteinGroupIDs
#' 
#' @param db database access 'handle'
#' @param proteinGroupIDs the proteinGroupIDs usually come from the
#'  TargetProtein Table. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the
#'  TargetProteinGroupsTargetPeptideGroups table or a character string
#'  specifying a SQL query
#' @note to get the proteinpeptidelink (table =
#'  "TargetProteinGroupsTargetPeptideGroups"). In goes "ProteinGroupID" from
#'  the table "TargetProteins" (Note: it's possible to use a c(,,,) to get
#'  the result for a number of proteins at the same time). The result is a
#'  list of numbers which are the "TargetProteinGroupsProteinGroupID" in the
#'  "TargetPeptideGroups" table
#'  @export
dbGetPeptideIDs <- function(db, ProteinGroupIDs, SQL = FALSE){
  if (isClass(ProteinGroupIDs,"data.frame")){
    # if so then assumed to be output from dbGetProteinTable
    # for speed set columnNames = "ProteinGroupID"
    ProteinGroupIDs <-
      as.character(ProteinGroupIDs$ProteinGroupID)
  } else {
    if (!is.character(ProteinGroupIDs)){
      ProteinGroupIDs <- as.character(ProteinGroupIDs)
    }
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetProteinGroupsTargetPeptideGroups",
    columnNames = "TargetPeptideGroupsPeptideGroupID",
    filtering = paste(c(" WHERE TargetProteinGroupsProteinGroupID IN (",
                        paste(c("'",
                                paste(ProteinGroupIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

#' helper function to determine if object == whichClass or a descendant of
#'  whichClass
#' 
#' @param object a data object of some class
#' @param whichClass character string: class name to be tested
#' @return TRUE or FALSE
isClass <- function(object, whichClass){
  return(whichClass %in% class(object))
}

#' get the paptide table belonging defined by PeptideIDs
#'
#' @param db database access 'handle'
#' @param PeptideIDs the peptideIDs to be retrieved. This can be in numeric or
#'  character vector format OR the output from the dbGetPeptideIDs function
#'  (a data.frame with column "TargetPeptideGroupsPeptideGroupID")
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param masterProtein use the IsMasterProtein column to be zero,
#'  default == TRUE. If more advanced filtering is needed, use db_getTable()
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid values are a single character
#'  string ("ASC" or "DESC") or a character vector of the same length as the
#'  columnNames vector containing a series of "ASC" or "DESC"
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return a data.frame containing requested data from the peptide table or
#'  a character string specifying a SQL query
#' @export
dbGetPeptideTable <- function(db,
                              PeptideIDs = NA,
                              columnNames = NA,
                              masterProtein = TRUE,
                              sortOrder = NA,
                              SQL = FALSE){
  if (identical(PeptideIDs,NA)){
    return(
      dbGetTable(
        db = db,
        tableName = "TargetPeptideGroups",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE PeptideGroupID IN ",
            "(SELECT TargetPeptideGroupsPeptideGroupID FROM ",
            "TargetProteinGroupsTargetPeptideGroups WHERE ",
            "TargetProteinGroupsProteinGroupID IN (SELECT ",
            "ProteinGroupIDs FROM TargetProteins ",
            ifelse(masterProtein,
                   "WHERE IsMasterProtein = 0",""),
            "))"),
          collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
  } else {
    if (isClass(PeptideIDs,"data.frame")){
      # if so then assumed to be output from dbGetPeptideIDs
      PeptideIDs <- as.character(PeptideIDs$TargetPeptideGroupsPeptideGroupID)
    } else {
      if (!is.character(PeptideIDs)){
        PeptideIDs <- as.character(PeptideIDs)
      }
    }
    return(
      dbGetTable(
        db = db,
        tableName = "TargetPeptideGroups",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE PeptideGroupID IN ",
            "(", paste(PeptideIDs, collapse = ","),")"), collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
  }
}

#' get the PsmID's from (a set of) PeptideGroupIDs
#' 
#' @param db database access 'handle'
#' @param PeptideGroupIDs the PeptideGroupIDs usually come from the
#'  TargetPeptideGroups Table. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the
#'  TargetPsmsTargetPeptideGroups  table or a character string specifying
#'  a SQL query
#'  @export
dbGetPsmIDs <- function(db, PeptideGroupIDs, SQL = FALSE){
  if (isClass(PeptideGroupIDs,"data.frame")){
    # if so then assumed to be output from dbGetTable
    # for speed set columnNames = "PeptideGroupID"
    PeptideGroupIDs <- as.character(PeptideGroupIDs$PeptideGroupID)
  } else {
    if (!is.character(PeptideGroupIDs)){
      PeptideGroupIDs <- as.character(PeptideGroupIDs)
    }
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetPsmsTargetPeptideGroups",
    columnNames = "TargetPsmsPeptideID",
    filtering = paste(c(" WHERE TargetPeptideGroupsPeptideGroupID IN (",
                        paste(c("'",
                                paste(PeptideGroupIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}


#' get the PSM table belonging to the PsmIDs
#'
#' @param db database access 'handle'
#' @param PsmIDs the PsmIDs to be retrieved. This can be in numeric or
#'  character vector format OR the output from the dbGetPsmIDs function
#'  (a data.frame with column "TargetPsmsPeptideID")
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param masterProtein use the IsMasterProtein column to be zero,
#'  default == TRUE. If more advanced filtering is needed, use db_getTable()
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid values are a single character
#'  string ("ASC" or "DESC") or a character vector of the same length as the
#'  columnNames vector containing a series of "ASC" or "DESC"
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return a data.frame containing requested data from the peptide table or
#'  a character string specifying a SQL query
#' @export
dbGetPsmTable <- function(db,
                          PsmIDs = NA,
                          columnNames = NA,
                          masterProtein = TRUE,
                          sortOrder = NA,
                          SQL = FALSE){
  if (identical(PsmIDs,NA)){
    return(
      dbGetTable(
        db = db,
        tableName = "TargetPsms",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE PeptideID IN ",
            "(SELECT TargetPsmsPeptideID FROM TargetPsmsTargetPeptideGroups",
            " WHERE TargetPeptideGroupsPeptideGroupID IN ",
            "(SELECT PeptideGroupID FROM TargetPeptideGroups",
            " WHERE PeptideGroupID IN ",
            "(SELECT TargetPeptideGroupsPeptideGroupID FROM ",
            "TargetProteinGroupsTargetPeptideGroups WHERE ",
            "TargetProteinGroupsProteinGroupID IN (SELECT ",
            "ProteinGroupIDs FROM TargetProteins ",
            ifelse(masterProtein,
                   "WHERE IsMasterProtein = 0",""),
            "))))"),
          collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
    # alternative: (not full SQL)
    # dbGetPsmTable(db = db1,
    #               PsmIDs = dbGetPeptideTable(db = db1,
    #                          PeptideIDs = dbGetProteinTable(db = db1,
    #                                   columnNames = "ProteinGroupIDs") %>%
    #           dbGetPeptideIDs(db = db1), columnNames = "PeptideGroupID") %>%
    #          dbGetPsmIDs(db = db1)
    # )
    # )
  } else {
    if (isClass(PsmIDs,"data.frame")){
      # if so then assumed to be output from dbGetPsmIDs
      PsmIDs <- as.character(PsmIDs$TargetPsmsPeptideID)
    } else {
      if (!is.character(PsmIDs)){
        PsmIDs <- as.character(PsmIDs)
      }
    }
    return(
      dbGetTable(
        db = db,
        tableName = "TargetPsms",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE PeptideID IN ",
            "(", paste(PsmIDs, collapse = ","),")"), collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
  }
}


#' get the ConsensusID's from (a set of) PeptideGroupIDs
#' 
#' @param db database access 'handle'
#' @param PeptideGroupIDs the PeptideGroupIDs usually come from the
#'  TargetPeptideGroups Table. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the
#'  TargetPsmsTargetPeptideGroups  table or a character string specifying
#'  a SQL query
#'  @export
dbGetConsensusIDs <- function(db, PeptideGroupIDs, SQL = FALSE){
  if (isClass(PeptideGroupIDs,"data.frame")){
    # if so then assumed to be output from dbGetTable
    # for speed set columnNames = "PeptideGroupID"
    PeptideGroupIDs <- as.character(PeptideGroupIDs$PeptideGroupID)
  } else {
    if (!is.character(PeptideGroupIDs)){
      PeptideGroupIDs <- as.character(PeptideGroupIDs)
    }
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetPeptideGroupsConsensusFeatures",
    columnNames = "ConsensusFeaturesId",
    filtering = paste(c(" WHERE TargetPeptideGroupsPeptideGroupID IN (",
                        paste(c("'",
                                paste(PeptideGroupIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

#' get the Consensus Features table belonging to the ConsensusIDs
#'
#' @param db database access 'handle'
#' @param ConsensusIDs the PsmIDs to be retrieved. This can be in numeric or
#'  character vector format OR the output from the dbGetConsensusIDs function
#'  (a data.frame with column "ConsensusFeaturesId")
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param masterProtein use the IsMasterProtein column to be zero,
#'  default == TRUE. If more advanced filtering is needed, use db_getTable()
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid values are a single character
#'  string ("ASC" or "DESC") or a character vector of the same length as the
#'  columnNames vector containing a series of "ASC" or "DESC"
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return a data.frame containing requested data from the peptide table or
#'  a character string specifying a SQL query
#' @export
dbGetConsensusTable <- function(db,
                          ConsensusIDs = NA,
                          columnNames = NA,
                          masterProtein = TRUE,
                          sortOrder = NA,
                          SQL = FALSE){
  if (identical(ConsensusIDs,NA)){
    return(
      dbGetTable(
        db = db,
        tableName = "ConsensusFeatures",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE Id IN ",
       "(SELECT ConsensusFeaturesId FROM TargetPeptideGroupsConsensusFeatures",
            " WHERE TargetPeptideGroupsPeptideGroupID IN ",
            "(SELECT PeptideGroupID FROM TargetPeptideGroups",
            " WHERE PeptideGroupID IN ",
            "(SELECT TargetPeptideGroupsPeptideGroupID FROM ",
            "TargetProteinGroupsTargetPeptideGroups WHERE ",
            "TargetProteinGroupsProteinGroupID IN (SELECT ",
            "ProteinGroupIDs FROM TargetProteins ",
            ifelse(masterProtein,
                   "WHERE IsMasterProtein = 0",""),
            "))))"),
          collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
      # alternative: (not full SQL)
      # dbGetConsensusTable(db = db1, ConsensusIDs =  
      #                       dbGetPeptideTable(db = db1,
      #                                         PeptideIDs = dbGetProteinTable(db = db1,
      #                                                                        columnNames = "ProteinGroupIDs") %>%
      #                                           dbGetPeptideIDs(db = db1), columnNames = "PeptideGroupID") %>%
      #                       dbGetConsensusIDs(db = db1)
      # )
    )
  } else {
    if (isClass(ConsensusIDs,"data.frame")){
      # if so then assumed to be output from dbGetConsensusIDs
      ConsensusIDs <- as.character(ConsensusIDs$ConsensusFeaturesId)
    } else {
      if (!is.character(ConsensusIDs)){
        ConsensusIDs <- as.character(ConsensusIDs)
      }
    }
    return(
      dbGetTable(
        db = db,
        tableName = "ConsensusFeatures",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE Id IN ",
            "(", paste(ConsensusIDs, collapse = ","),")"), collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
  }
}

#' get the names of the identification types used in the database
#'
#' @param db database access 'handle'
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return data.frame with a single column: "GroupName"
#' @export
proteinIDTypes <- function(db, SQL = FALSE){
  return(
    dbGetTable(
      db = db,
      tableName = "ProteinIdentificationGroups",
      columnNames = "GroupName",
      SQL = SQL
    )
  )
}

#' converts character string date into date/time format
#' (mdy hms format)
#' 
#' @param theDate character string to be converted
#' (can be vectorized)
#' @returns date in mdy hms format
#' @export
thermo.date <- function(theDate){
  return(lubridate::mdy_hms(theDate))
}

#' converts character string date into date/time format
#' (ymd hms format)
#' 
#' @param theDate character string to be converted
#' (can be vectorized)
#' @returns date 
#' @export
system.date <- function(theDate){
  return(lubridate::ymd_hms(theDate))
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

#' get the table with info on the files used in the search from the database
#'
#' @param db database access 'handle'
#' @param type allows for selection of the FileTypes
#'  default = "XCaliburRawFile"
#' @param dates allows transformation of the date/time strings from te database
#'  to be transformed into proper data/time fields. Default function used is
#'  thermo.date. If no transformation is required, use na.date
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return data.frame
#' @export
MSfileInfo <- function(db, type = "XcaliburRawfile",
                       dates = thermo.date, SQL = FALSE){
  if (nchar(type) == 0){
    if (SQL){
      return(
        dbGetTable(
          db = db,
          tableName = "WorkFlowInputFiles",
          SQL = TRUE))
    } else {
      return(
        dbGetTable(
          db = db,
          tableName = "WorkFlowInputFiles") %>%
          dplyr::mutate(CreationDate = dates(CreationDate)))
    }
  } else {
    if (SQL){
      return(
        dbGetTable(
          db = db,
          tableName = "WorkFlowInputFiles",
          filtering = paste(c(" WHERE FileType IN (",
                              paste(c("'",
                                      paste(type,collapse = "','"),
                                      "'"),
                                    collapse = ""),
                              ")"),
                            collapse = ""),
          SQL = TRUE
        )
      )
    } else {
      return(
        dbGetTable(
          db = db,
          tableName = "WorkFlowInputFiles",
          filtering = paste(c(" WHERE FileType IN (",
                              paste(c("'",
                                      paste(type,collapse = "','"),
                                      "'"),
                                    collapse = ""),
                              ")"),
                            collapse = "")) %>%
          dplyr::mutate(CreationDate = dates(CreationDate)))
    }
  }
}

#' get the table with info on the search itself from the database
#'
#' @param db database access 'handle'
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return data.frame
#' @export
SearchInfo <- function(db, SQL = FALSE){
  return(
    dbGetTable(
      db = db,
      tableName = "WorkflowMessages",
      SQL = SQL))
}

#' get the total search time  from the database
#'
#' @param db database access 'handle'
#' @return numeric: search time in seconds
#' @export
totalSearchTime <- function(db){
  temptbl <-
    pool::dbGetQuery(db,
                     "SELECT * FROM (SELECT Time FROM WorkflowMessages
                     ORDER BY Time DESC LIMIT 1) UNION
                     SELECT * FROM (SELECT Time FROM WorkflowMessages
                     ORDER BY Time ASC LIMIT 1)")
  return((temptbl$Time[2] - temptbl$Time[1])/10000000)
}

#' for translation of the value (1..5) in the psmTable to words
#' (like in Proteome Discoverer). <...> --> means not encountered/undefined/
#' no inference
#' 
#' @format character vector
PsmAmbiguity <- c("<Not Considered>","Rejected","<Ambiguous>",
                  "Selected","Unambiguous")


#' for translation of the value (0..4) in the proteinTable to words
#' (like in Proteome Discoverer). <...> --> means not encountered/undefined/
#' no inference
#' @format character vector
IsMasterProtein <- c("Master Protein","Master Protein Candidate","<Undefined>",
                     "Master Protein Rejected","Master Protein Rejected")
