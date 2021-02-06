
# ---- database (related) functions ----
# ---- Proteome Discoverer specific code ----

# Function to transform the table target proteins into a proper data.frame
# note: the argument dbtble must be the in tbl() format. Under normal circumstances
# it is the "TargetProteins" table in the .pdResult file
# using IsMasterProtein = TRUE will give the table that proteome discoverer automatically
# displays (if filters are default)


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
#'  default = NA, (no sorting). Other valid values are a single character
#'  string ("ASC" or "DESC") or a character vector of the same length as the
#'  columnNames vector containing a series of "ASC" or "DESC"
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
  if (!identical(sortOrder,NA) & !identical(columnNames,NA)){
    sortColumns = paste(unlist(purrr::map2(columnNames,
                                          sortOrder,
                                          ~paste(.x,.y,collapse = ""))),
                        collapse = ", ")
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
#'  TargetProteinTable. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the protein table or
#' a character string specifying a SQL query
#' @note to get the proteinpeptidelink (table =
#'  "TargetProteinGroupsTargetPeptideGroups"). In goes "ProteinGroupID" from
#'  the table "TargetProteins" (Note: it's possible to use a c(,,,) to get
#'  the result for a number of proteins at the same time). The result is a
#'  list of numbers which are the "TargetProteinGroupsProteinGroupID" in the
#'  "TargetPeptideGroups" table
#'  @export
dbGetPeptideIDs <- function(db, ProteinGroupIDs, SQL = FALSE){
  if (!is.character(ProteinGroupIDs)){
    ProteinGroupIDs <- as.character(ProteinGroupIDs)
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

# to get the peptidepsmlink (table "TargetPsmsTargetPeptideGroups")
# in goes "PeptideGroupID" from the table "TargetPeptideGroups" (Note: it's possible to use a c(,,,)
# result ia a list of numbers which are the "PeptideID" in the "TargetPsms" table
db_get_psmIDs <- function(db, PeptideGroupIDs){
  return((tbl(db, "TargetPsmsTargetPeptideGroups") %>% filter(TargetPeptideGroupsPeptideGroupID %in% PeptideGroupIDs) %>% collect())$TargetPsmsPeptideID)
}

# to get the psm table belonging to a peptide
# in goes the psmIDs coming from ef db_get_psmIDs
# out comes a data.frame which is a subset of "TargetPsms"
db_get_psmTable <- function(db, PsmIDs){
  return(tbl(db,"TargetPsms") %>% filter(PeptideID %in% PsmIDs) %>% collect())
}

# to get the peptideconsensusfeatureslink (table "TargetPeptideGroupsConsensusFeatures")
# in goes "PeptideGroupID" from the table "TargetPeptideGroups" (Note: it's possible to use a c(,,,)
# result ia a list of numbers which are the "PeptideID" in the "TargetConsensusFeatures" table
db_get_consensusIDs <- function(db, PeptideGroupIDs){
  return((tbl(db, "TargetPeptideGroupsConsensusFeatures") %>% filter(TargetPeptideGroupsPeptideGroupID %in% PeptideGroupIDs) %>% collect())$ConsensusFeaturesId)
}

# to get the consensus features table belonging to a peptide
# in goes the psmIDs coming from ef db_get_consensusIDs
# out comes a data.frame which is a subset of "TargetConsensusFeatures"
db_get_consensusTable <- function(db, consensusIDs){
  return(tbl(db,"ConsensusFeatures") %>% filter(Id %in% consensusIDs) %>% collect() %>% df_transform_raws())
}

