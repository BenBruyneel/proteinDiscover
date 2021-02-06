
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
#' @return a data.frame containing requested data from the protein table or
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
#' @param tableName used to pass on the name of the table containing the data,
#'  default = "TargetProteins"
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
                               tableName = "TargetProteins",
                               columnNames = NA,
                               masterProtein = TRUE,
                               sortOrder = NA,
                               SQL = FALSE){
  return(dbGetTable(
    db = db,
    tableName = tableName,
    columnNames = columnNames,
    filtering = ifelse(masterProtein,
                       " WHERE IsMasterProtein = 0",""),
    sortOrder = sortOrder,
    SQL = SQL))
}

# Function to trasnform the table target peptide groups into a proper data.frame
# note: the argument dbtble must be the in tbl() format. Under normal circumstances
# it is the "TargetPeptideGroups" table in the .pdResult file
# note: removeNA is to remove a lot of rows present in the 'raw' table
# during development of the code it was discovered that there are a lot
# of rows with (a lot of ) NA values in their fields. Removing these will
# give a table which is displayed by default when using Proteome Discoverer
db_get_peptides <- function(dtbl, removeNA = TRUE){
  if (removeNA){
    return(df_transform_raws(dtbl %>% filter(!is.na(MasterProteinAccessions)) %>% collect()))
  } else {
    return(df_transform_raws(dtbl %>% collect()))
  }
}

# Function to trasnform the table target psms into a proper data.frame
# note: the argument dbtble must be the in tbl() format. Under normal circumstances
# it is the "TargetPsms" table in the .pdResult file
# note: no raw vector columns were  discovered in this type of table (so far)
# so df_transform_raws() is not used
# note: removeNA is to remove a lot of rows present in the 'raw' table
# during development of the code it was discovered that there are a lot
# of rows with (a lot of ) NA values in their fields. Removing these will
# give a table which is displayed by default when using Proteome Discoverer
db_get_PSMs <- function(dtbl, removeNA = TRUE){
  if (removeNA){
    return(dtbl %>% filter(!is.na(MasterProteinAccessions)) %>% collect())
  } else {
    return(dtbl %>% collect())
  }
}

# to get the proteinpeptidelink (table = "TargetProteinGroupsTargetPeptideGroups")
# in goes "ProteinGroupID" from the table "TargetProteins" (Note: it's possible to use a c(,,,) to get the result for a number of proteins at the same time)
# result ia a list of numbers which are the "TargetProteinGroupsProteinGroupID" in the "TargetPeptideGroups" table
db_get_peptideIDs <- function(db, ProteinGroupIDs){
  return((tbl(db, "TargetProteinGroupsTargetPeptideGroups") %>% filter(TargetProteinGroupsProteinGroupID %in% ProteinGroupIDs) %>% collect())$TargetPeptideGroupsPeptideGroupID)
}

# to get the peptide table belonging to a protein
# in goes the PeptideIDs coming from eg db_get_peptideIDs
# out comes a data.frame which is a subset of "TargetPeptideGroups"
db_get_peptideTable <- function(db, PeptideIDs){
  return(tbl(db,"TargetPeptideGroups") %>% filter(PeptideGroupID %in% PeptideIDs) %>% collect() %>% df_transform_raws())
}

# to get the complete peptide table = all peptides belonging to all proteins
# in goes only the database, out comes the table WITH raw vector columns still intact
# advantage of this function is that is 'pure' SQL
# note: out comes tibble() to make sure for compatibility with df_transform_raws()
db_get_all_peptideTable <- function(db){
  return(tibble(dbGetQuery(db,"SELECT * FROM TargetPeptideGroups WHERE PeptideGroupID IN
                    (SELECT TargetPeptideGroupsPeptideGroupID FROM TargetProteinGroupsTargetPeptideGroups
                     WHERE  TargetProteinGroupsProteinGroupID IN
                      (SELECT ProteinGroupIDs FROM TargetProteins WHERE IsMasterProtein = 0))"))
  )
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

