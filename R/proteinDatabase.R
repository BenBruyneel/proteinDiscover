# ---- General ----

#' Wrapper around pool::dbPool(): opens a database
#'
#' @param  fileName  a character vector specifying the name and location
#'                   of the database
#' @param  drv defines database connection type, default = RSQLite::SQLite()
#' @param  ... to pass on additional parameters to pool::dbPool, exmples are
#'             host = "shiny-demo.csa7qlmguqrf.us-east-1.rds.amazonaws.com"
#'             username = "guest"
#'             password = "guest"
#'
#' @return database access 'handle'
#' @note if no file with the name 'fileName' exists, then it will be created
#' (but obviously it will be empty, so most further commands will fail)
#' @note if fileName == ":memory:" the database will be an in-memory database
#' @export
dbOpen <- function(fileName, drv = RSQLite::SQLite(), ...){
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
dbClose <- function(db){
  pool::poolClose(db)
}

#' get a table from a .pdResult file
#'
#' @param db database access 'handle'
#' @param tableName used to pass on the name of the table containing the data
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param filtering allows for " WHERE <expression>" additions to the SQL
#'  statement default = " " (no filtering). Note: always put a space (" ")
#'  before any statement
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid value is a character character
#'  vector of columnNames to be used for sorting string (with "ASC" or "DESC"
#'  if needed) 
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#' @return a data.frame containing requested data from a database table or
#'  a character string specifying an SQL query
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

#' get the names of the identification types (sequest HT etc) used in the
#'  database
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
      tempResult <- dbGetTable(
        db = db,
        tableName = "WorkFlowInputFiles")
    }
  } else {
    tempResult <- dbGetTable(
      db = db,
      tableName = "WorkFlowInputFiles",
      filtering = paste(c(" WHERE FileType IN (",
                          paste(c("'",
                                  paste(type,collapse = "','"),
                                  "'"),
                                collapse = ""),
                          ")"),
                        collapse = ""), SQL = SQL)
  }
  if (!is.Class(tempResult,"character")){
    tempResult$CreationDate__ <- dates(tempResult$CreationDate)
    tempResult <- tempResult %>%
      dplyr::select(-dplyr::all_of(c("CreationDate")))
    colnames(tempResult)[which(colnames(tempResult) == "CreationDate__")] <- "CreationDate"
  }
  return(tempResult)
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

#' function to retrieve the acquisition date & time of the files used to
#'  generate the pdResult file
#'
#' @param db database access 'handle'
#' @param useAmPm logical, influences what default format is used. Ignored if a
#'  format is specified
#' @param format character vector specifying the format of the resulting
#'  POSIXct/POSIXt object. See \code{\link[base]{strptime}} for more info
#' 
#' @note this function is essentially a wrapper around 
#'  \code{\link[proteinDiscover]{studyDefinitionFileSets}}
#'                               
#' @returns one or more POSIXct/POSIXt object(S)
#' @export
getAcquistionDateTime <- function(db, useAmPm = TRUE,
                                  format = ifelse(useAmPm,
                                                  "%m/%d/%Y %I:%M:%S %p",
                                                  "%m/%d/%Y %H:%M:%S")){
  return(studyDefinitionFileSets(analysisDefinition(db = db))$FileTime %>%
           lubridate::as_datetime(format = format))
}

#' function to retrieve the acquisition date of the files used to generate the
#'  pdResult file
#'
#' @param db database access 'handle'
#' 
#' @note this function is essentially a wrapper around 
#'  \code{\link{getAcquistionDateTime}}
#'  
#' @returns one or more POSIXct/POSIXt object(S)
#' @export
getAcquistionDate <- function(db){
  return(getAcquistionDateTime(db = db, format = "%m/%d/%Y"))
}

# ---- Proteins ----

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
    tableName = tableNames("proteins"),
    columnNames = columnNames,
    filtering = ifelse(masterProtein,
                       " WHERE IsMasterProtein = 0",""),
    sortOrder = sortOrder,
    SQL = SQL))
}

#' Function to get protein information from the TargetProteins table on the
#'  basis of their UniqueSequenceID
#'
#' @param db database access 'handle'
#'
#' @param UniqueSequenceIDs specifies from which proteins to get info
#' @param columns character vector, specifies which columns to retrieve
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame#'
#'
#' @return a data.frame or a character vector (SQL)
#' @export
dbGetProteins <- function(db, UniqueSequenceIDs, columns = NA, SQL = FALSE){
  if (is.Class(UniqueSequenceIDs,"data.frame")){
    if ("UniqueSequenceID" %in% colnames(UniqueSequenceIDs)){
      UniqueSequenceIDs <- UniqueSequenceIDs$UniqueSequenceID
    } else {
      UniqueSequenceIDs <- UniqueSequenceIDs[,1]
    }
  }
  return(
    dbGetTable(db = db, tableName = tableNames("proteins"),
               columnNames = columns,
               filtering = paste(c(" WHERE UniqueSequenceID IN ('",
                                   paste(UniqueSequenceIDs, collapse = "','"),
                                   "') "),
                                 collapse = "")))
}

#' A bit more advanced version of \code{\link{dbGetProteinTable}} which allows
#'  for filtering (via SQL). Note that filtering raw columns (BLOB's) will
#'  not work properly
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
#' @param filtering SQL statement to be used for filtering of the query. The
#'  IsMasterProtein column is already covered when masterProtein is set to TRUE
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame containing requested data from the protein table or
#'  a character string specifying an SQL query
#' @export
dbGetProteinFiltered <- function(db, columnNames = NA, masterProtein = FALSE,
                                 sortOrder = NA, filtering = NA, SQL = FALSE){
  return(dbGetTable(db = db, tableName = tableNames("proteins"), 
                    columnNames = columnNames,
            filtering = paste(c(ifelse(!identical(filtering,NA) | masterProtein,
                                               " WHERE ", " "),
                                        ifelse(masterProtein,
                                               " IsMasterProtein = 0"," "),
                                ifelse(!identical(filtering,NA) & masterProtein,
                                               " AND ", " "),
                                        ifelse(identical(filtering, NA),
                                               " ", filtering), " "),
                                      collapse = ""),
                    sortOrder = sortOrder,
                    SQL = SQL))
}

#' Function to retrieve the UniqueSequenceID's based on the accession field of
#'  the proteinTable. Essentially a wrapper for
#'  \code{\link{dbGetProteinFiltered}}
#'
#' @param db database access 'handle'
#' @param accession accession(s) of the proteins
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame or a character vector (SQL)
#' 
#' @export
dbGetProteinUniqueSequenceIDs <- function(db, accession = NA, SQL = FALSE){
  return(dbGetProteinFiltered(db = db,
                              columnNames = "UniqueSequenceID",
                              filtering = paste(c("Accession IN ('",
                                                  paste(accession,
                                                        collapse = "','"),
                                                  "')"),
                                                collapse = ""),
                              SQL = SQL))
}

# ---- Protein Grouping ----

#' Gets the ProteinGroup information from the TargetProteinGroups table
#'
#' @param db database access 'handle'
#' @param proteinGroupIDs specifies which protein groups to get, these values
#'  can come from eg the protein table
#' @param columns character vector, specifies which columns to retrieve
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame#'
#'
#' @return a data.frame or a character vector (SQL)
#' 
#' @export
dbGetProteinGroups <- function(db, proteinGroupIDs, columns = NA, SQL = FALSE){
  if (is.Class(proteinGroupIDs,"data.frame")){
    if ("proteinGroupIDs" %in% colnames(proteinGroupIDs)){
      proteinGroupIDs <- proteinGroupIDs$proteinGroupID
    } else {
      proteinGroupIDs <- proteinGroupIDs[,1]
    }
  }
  if (sum(grepl(proteinGroupIDs, pattern = ";")) > 0){
    proteinGroupIDs <- unique(unlist(lapply(proteinGroupIDs,
                                            function(x){strsplit(x, split = ";")})))
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetProteinGroups",
    columnNames = columns,
    filtering = paste(c(" WHERE ProteinGroupID IN (",
                        paste(c("'",
                                paste(proteinGroupIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

#' Retrieve the ProteinGroupID's of proteins via their UniqueSequenceID's
#'
#' @param db database access 'handle'
#' @param proteinUniqueIDs the UniqueSequenceID's for which the proteinGroupID's
#'  are to be retrieved. Usually these UniqueSequenceID's will come from a
#'  protein table
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame#'
#'  
#' @return a data.frame or a character vector (SQL)
#' 
#' @note the output of this is meant to serve as input for the
#'  \code{\link{dbGetProteinGroups}} function
#' @export
dbGetProteinGroupIDs <- function(db, proteinUniqueIDs, SQL = FALSE){
  if (is.Class(proteinUniqueIDs,"data.frame")){
    proteinUniqueIDs <- bit64::as.integer64(proteinUniqueIDs$UniqueSequenceID)
  } else {
    proteinUniqueIDs <- bit64::as.integer64(proteinUniqueIDs)
  }
  return(
    dbGetTable(
      db = db,
      tableName = "TargetProteinGroupsTargetProteins",
      columnNames = "TargetProteinGroupsProteinGroupID",
      filtering = paste(c(" WHERE TargetProteinsUniqueSequenceID IN (",
                          paste(c("'",
                                  paste(proteinUniqueIDs,collapse = "','"),
                                  "'"),
                                collapse = ""),
                          ")"),
                        collapse = ""),
      sortOrder = NA,
      SQL = SQL)
  )
}

#' Function to get proteinUniqueID's from a (set of) protein groupID's
#'  (eg from a proteinGroup tables, or dbGetProteinGroupIDs). This allows for
#'  getting all proteins (also non-master proteins) which together make up
#'  a protein group. Normally only the master protein is shown in a protein
#'  table
#'
#' @param db database access 'handle'
#' @param proteinGroupIDs the protein group(s) for which the UniqueSequenceID's
#'  should be retrieved
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame#'
#'
#' @return a data.frame or a character vector (SQL)#' 
#' 
#' @note every protein in the protein table has a ProteinGroupID & a
#'  UniqueSequenceID. The UniqueSequenceID is untique to the protein. A protein
#'  group may contain more than one protein (and thus also more than one
#'  UniqueSequenceID)
#' 
#' @export
dbGetProteinIDs <- function(db, proteinGroupIDs, SQL = FALSE){
  if (is.Class(proteinGroupIDs,"data.frame")){
    proteinGroupIDs <- proteinGroupIDs$proteinGroupID
  }
  if (sum(grepl(proteinGroupIDs, pattern = ";")) > 0){
    proteinGroupIDs <- unique(unlist(lapply(proteinGroupIDs,
                                            function(x){strsplit(x, split = ";")})))
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetProteinGroupsTargetProteins",
    columnNames = "TargetProteinsUniqueSequenceID",
    filtering = paste(c(" WHERE TargetProteinGroupsProteinGroupID IN (",
                        paste(c("'",
                                paste(proteinGroupIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

# ---- Annotation ----

#' Function to get the functional group annotation group ID's for proteins.
#'  This function does essentially the reverse of
#'  \code{\link{dbGetAnnotatedProteins}}. The output of this function can serve
#'  as the input for \code{\link{dbGetAnnotationGroups}}
#'
#' @param db database access 'handle'
#' @param UniqueSequenceIDs the UniqueSequenceID's (unique protein identifier),
#'  usually coming from protein table
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame or a character vector (SQL)
#' @export
dbGetProteinAnnotationGroupIDs <- function(db, UniqueSequenceIDs, SQL = FALSE){
  if (is.Class(UniqueSequenceIDs,"data.frame")){
    UniqueSequenceIDs <- UniqueSequenceIDs$UniqueSequenceID
  }
  return(
    dbGetTable(db = db, tableName = "AnnotationProteinGroupsTargetProteins",
               columnNames = "AnnotationProteinGroupsProteinAnnotationGroupID",
               filtering = paste(c(" WHERE TargetProteinsUniqueSequenceID IN ('",
                                   paste(UniqueSequenceIDs , collapse = "','"),
                                   "') "),
                                 collapse = ""), SQL = SQL))
}

#' Function to get the UniqueSequenceID's for proteins which are in an protein
#'  annotation group. Essentially does the reverse of
#'  \code{\link{dbGetProteinAnnotationGroupIDs}}. The output of this function
#'  can serve as the input for \code{\link{dbGetProteins}}
#'  
#' @param db database access 'handle'
#' @param ProteinAnnotationGroupIDs the protein annotation group ID's for which
#'  to get the UniqueSequenceID's
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame or a character vector (SQL)
#' @export
dbGetAnnotatedProteins <- function(db, ProteinAnnotationGroupIDs, SQL = FALSE){
  if (is.Class(ProteinAnnotationGroupIDs,"data.frame")){
    ProteinAnnotationGroupIDs <-
      ProteinAnnotationGroupIDs$ProteinAnnotationGroupID
  }
  return(
    dbGetTable(db = db, tableName = "AnnotationProteinGroupsTargetProteins",
               columnNames = "TargetProteinsUniqueSequenceID",
               filtering = paste(
               c(" WHERE AnnotationProteinGroupsProteinAnnotationGroupID IN ('",
                   paste(ProteinAnnotationGroupIDs , collapse = "','"),
                   "') "),
                 collapse = ""), SQL = SQL))
}

#' Function to get the info for (protein) annotation groups. Takes eg
#'  \code{\link{dbGetProteinAnnotationGroupIDs}} as input
#'
#' @param db database access 'handle'
#' @param ProteinAnnotationGroupIDs the protein annotation group ID's for which
#'  to get information
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame or a character vector (SQL)
#' @export
dbGetAnnotationGroups <- function(db, ProteinAnnotationGroupIDs = NA,
                                  columnNames = NA, SQL = FALSE){
  if (is.Class(ProteinAnnotationGroupIDs,"data.frame")){
 ProteinAnnotationGroupIDs <- ProteinAnnotationGroupIDs$ProteinAnnotationGroupID
  }
  return(
    dbGetTable(db = db, tableName = "AnnotationProteinGroups",
               columnNames = columnNames,
              filtering = ifelseProper(identical(ProteinAnnotationGroupIDs, NA),
                                        " ",
                                paste(c(" WHERE ProteinAnnotationGroupID IN ('",
                                        paste(ProteinAnnotationGroupIDs,
                                              collapse = "','"),
                                        "') "),
                                      collapse = "")),
               SQL = SQL))
}

#' Get Group Annotation information from the table: AnnotationProteinGroups.
#'  This can be done via the GroupAnnotationAccession or via the description of
#'  an annotation. When using the Description it's possible to use the SQL 
#'  'like'
#'
#' @param db database access 'handle'
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param GroupAnnotationAccession identification of the annotation, usually
#'  something like GO:....  (gene ontology) or pF.... (protein family). Note
#'  that when this argument is not NAm the arguments dealing with description
#'  etc are ignored 
#' @param description character vector specifying a word or sequence of word
#'  which is to be selected. If the 'like' argument is TRUE then it doesn't need
#'  to be exactly the same as the GroupAnnotationDescription field/column (in
#'  most cases the 'like' argument should be set to TRUE !)
#' @param UpperCase if set to TRUE then BOTH description and the
#'  GroupAnnotationDescription field/column are entirely put to uppercase in the
#'  SQL used for the query. Note that if both UpperCase and LowerCase are TRUE,
#'  then UpperCase is used 
#' @param LowerCase if set to TRUE then BOTH description and the
#'  GroupAnnotationDescription field/column are entirely put to lowercase in the
#'  SQL used for the query.
#' @param like if set to TRUE then the SQL 'LIKE' in stead of 'IN' is used to
#'  query the data. This only applies when the argument 'discription' is used. 
#'  This is ignored when 'GroupAnnotationAccession' is used. If like = TRUE,
#'  then using eg 'locomotion' will result in the SQL query being: WHERE ...
#'   LIKE '%locomotion%' (the %'s are coming from likePro & likePost). The
#'   resulting table will give all rows, where the description part contains
#'   'locomotion'. If like = FALSE, then only rows where the description exactly
#'   matches 'locomotion' will be selected. It's also possible to use the '_'
#'   (underscore) to make the LIKE function more or less specific. See eg
#'   \href{https://www.w3schools.com/sql/sql_like.asp}{SQL LIKE Operator} for
#'   more info
#' @param likePre default is '%'. This character vector gets put in front of the
#'  'description' argument to facilitate (partial) matching. It's better to set
#'  to '' (empty string) when creating LIKE arguments directly via the
#'  'description' argument 
#' @param likePost  default is '%'. This character vector gets added to the end
#'  of the 'description' argument to facilitate (partial) matching. It's better
#'  to set to '' (empty string) when creating LIKE arguments directly via the
#'  'description' argument
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame or a character vector (SQL)
#' @export
dbGetAnnotationGroupsFiltered <- function(
    db, columnNames = NA,
    GroupAnnotationAccession = NA, description = NA,
    UpperCase = FALSE, LowerCase = FALSE,
    like = FALSE, likePre = "%", likePost = "%",
    # in case of like needs to single string
    SQL = FALSE){
  if ((length(description) > 1) & (like)){
    stop("Cannot combine multi-element 'description' argument with like = TRUE")
  }
  if (identical(GroupAnnotationAccession,NA)){
    if (!identical(description,NA)){
      description <- paste(c(" WHERE ",
                             ifelse(UpperCase,
                                    "UPPER(GroupAnnotationDescription)",
                                    ifelse(LowerCase,
                                           "LOWER(GroupAnnotationDescription)",
                                           "GroupAnnotationDescription")),
                             ifelse(like,
                                    paste0(" LIKE '", likePre) ,
                                    " IN ('"),
                             paste(ifelseProper(UpperCase,
                                                toupper(description),
                                                ifelseProper(LowerCase,
                                                           tolower(description),
                                                           description)),
                                   collapse = "','"),
                             ifelse(like,
                                    paste0(likePost,"'"),
                                    "')")),
                           collapse = "") 
      
      return(
        dbGetTable(db = db, tableName = "AnnotationProteinGroups",
                   columnNames = columnNames,
                   filtering = description,
                   SQL = SQL))
    } else {
      return(NA)
    }
  } else {
    if (is.Class(GroupAnnotationAccession,"data.frame")){
   GroupAnnotationAccession <- GroupAnnotationAccession$GroupAnnotationAccession
    }
    return(
      dbGetTable(db = db, tableName = "AnnotationProteinGroups",
                 columnNames = columnNames,
                 filtering = paste(c(" WHERE GroupAnnotationAccession IN ('",
                                     paste(GroupAnnotationAccession,
                                           collapse = "','"),
                                     "') "),
                                   collapse = ""),
                 SQL = SQL))
  }
}

# ---- Peptides & Peptide Spectral Matches (PSM's) ----

#' get the peptideID's from (a set of) proteinGroupIDs
#' 
#' @param db database access 'handle'
#' @param proteinGroupIDs the proteinGroupIDs usually come from the
#'  TargetProtein Table. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the
#'  TargetProteinGroupsTargetPeptideGroups table or a character string
#'  specifying an SQL query
#' @note to get the proteinpeptidelink (table =
#'  "TargetProteinGroupsTargetPeptideGroups"). In goes "ProteinGroupID" from
#'  the table "TargetProteins" (Note: it's possible to use a c(,,,) to get
#'  the result for a number of proteins at the same time). The result is a
#'  list of numbers which are the "TargetProteinGroupsProteinGroupID" in the
#'  "TargetPeptideGroups" table
#' @export
dbGetPeptideIDs <- function(db, proteinGroupIDs, SQL = FALSE){
  if (is.Class(proteinGroupIDs,"data.frame")){
    # if so then assumed to be output from dbGetProteinTable
    # for speed set columnNames = "ProteinGroupID"
    proteinGroupIDs <-
      as.character(proteinGroupIDs$ProteinGroupID)
  } else {
    if (!is.character(proteinGroupIDs)){
      proteinGroupIDs <- as.character(proteinGroupIDs)
    }
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetProteinGroupsTargetPeptideGroups",
    columnNames = "TargetPeptideGroupsPeptideGroupID",
    filtering = paste(c(" WHERE TargetProteinGroupsProteinGroupID IN (",
                        paste(c("'",
                                paste(proteinGroupIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
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
        tableName = tableNames("peptides"),
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
    if (is.Class(PeptideIDs,"data.frame")){
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
        tableName = tableNames("peptides"),
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
#'  an SQL query
#' @export
dbGetPsmIDs <- function(db, PeptideGroupIDs, SQL = FALSE){
  if (is.Class(PeptideGroupIDs,"data.frame")){
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
#' @param filtering allows for " WHERE <expression>" additions to the SQL
#'  statement default = " " (no filtering). Note: always put a space (" ")
#'  before any statement. If NA then no filtering is applied. Note that
#'  filtering is only used when the argument PsmIDs is not NA
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
                          filtering = "MasterProteinAccessions IS NOT NULL",
                          SQL = FALSE){
  if (identical(PsmIDs,NA)){
    return(
      dbGetTable(
        db = db,
        tableName = tableNames("psms"),
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
            ")))) AND MasterProteinAccessions IS NOT NULL"),
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
    if (is.Class(PsmIDs,"data.frame")){
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
        tableName = tableNames("psms"),
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE PeptideID IN ",
            "(", paste(PsmIDs, collapse = ","),")",
            ifelse(is.na(filtering),
                   "",
                   paste(c(" AND ", filtering), collapse =""))),
            collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
  }
}

# ---- Consensus & QuanSpectrum tables ----

#' get the ConsensusID's from (a set of) PeptideGroupIDs
#' 
#' @param db database access 'handle'
#' @param PeptideGroupIDs the PeptideGroupIDs usually come from the
#'  TargetPeptideGroups Table. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the
#'  TargetPeptideGroupsConsensusFeatures table or a character string specifying
#'  a SQL query
#' @export
dbGetConsensusIDs <- function(db, PeptideGroupIDs, SQL = FALSE){
  if (is.Class(PeptideGroupIDs,"data.frame")){
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
        tableName = tableNames("Consensus"),
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
    )
  } else {
    if (is.Class(ConsensusIDs,"data.frame")){
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
        tableName = tableNames("Consensus"),
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE Id IN ",
            "(", paste(ConsensusIDs, collapse = ","),")"), collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL)
    )
  }
}

#' get the SpectrumID's from (a set of) PeptideIDs
#' 
#' @param db database access 'handle'
#' @param PeptideIDs the PeptideIDs usually come from the
#'  TargetPsms Table. This can be in numeric or character vector format
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame  
#' @return a data.frame containing requested data from the
#'  TargetPsmsQuanSpectrumInfo table or a character string specifying
#'  an SQL query
#' @export
dbGetQuanSpectrumIDs <- function(db, PeptideIDs, SQL = FALSE){
  if (is.Class(PeptideIDs,"data.frame")){
    PeptideIDs <- as.character(PeptideIDs$PeptideID)
  } else {
    if (!is.character(PeptideIDs)){
      PeptideIDs <- as.character(PeptideIDs)
    }
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetPsmsQuanSpectrumInfo",
    columnNames = "QuanSpectrumInfoSpectrumID",
    filtering = paste(c(" WHERE TargetPsmsPeptideID IN (",
                        paste(c("'",
                                paste(PeptideIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

#' get the QuanSpectrumInfo table belonging to the SpectrumID's
#'
#' @param db database access 'handle'
#' @param SpectrumIDs the SpectrumID's to be retrieved. This can be in numeric
#'  or character vector format OR the output from the dbGetQuanSpectrumIDs
#'  function (a data.frame with column "QuanSpectrumInfoSpectrumID")
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
#' @return a data.frame containing requested data from the QuanSpectrumInfo
#'  table or a character string specifying a SQL query
#' @export
dbGetQuanSpectrumInfoTable <- function(db,
                                       SpectrumIDs = NA,
                                       columnNames = NA,
                                       masterProtein = TRUE,
                                       sortOrder = NA,
                                       SQL = FALSE){
  if (identical(SpectrumIDs,NA)){
    dbGetTable(
      db = db,
      tableName = "QuanSpectrumInfo",
      columnNames = columnNames,
      filtering = paste(
        c(" WHERE SpectrumID IN ",
          "(SELECT QuanSpectrumInfoSpectrumID FROM TargetPsmsQuanSpectrumInfo",
          " WHERE TargetPsmsPeptideID IN ",
          "(SELECT PeptideID FROM TargetPsms",
          " WHERE PeptideID IN ",
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
          ")))) AND MasterProteinAccessions IS NOT NULL))"),
        collapse = ""),
      sortOrder = sortOrder,
      SQL = SQL)
  } else {
    if (is.Class(SpectrumIDs, "data.frame")){
      SpectrumIDs <- as.character(SpectrumIDs$QuanSpectrumInfoSpectrumID)
    } else {
      if (!is.character(SpectrumIDs)){
        SpectrumIDs <- as.character(SpectrumIDs)
      }
    }
    return(
      dbGetTable(
        db = db,
        tableName = "QuanSpectrumInfo",
        columnNames = columnNames,
        filtering = paste(
          c(" WHERE SpectrumID IN ",
            "(", paste(SpectrumIDs, collapse = ","),")"), collapse = ""),
        sortOrder = sortOrder,
        SQL = SQL
      )
    )
  }
}

# ---- Modifications ----

#' function to get the modificationSite ID's from (a set of) proteinUniqueID's
#'
#' @param db database access 'handle'
#' @param proteinUniqueIDs the protein identifier for which the modificationSite
#'  ID's are to be fetched. This is a vector of one or more integer64 (package:
#'  bit64 ) values. In protein tables this is the UniqueSequenceUD column
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame containing the requested data from the
#'  TargetProteinsModificationSites table or a character string specifying an
#'  SQL query
#'  
#' @note the data from modificationSitesUd's in the result can be used to query
#'  the ModificationSites table via \code{\link{dbGetModificationsTable}}
#' @export
dbGetModificationsSitesIDs <- function(db, proteinUniqueIDs, SQL = FALSE){
  if (is.Class(proteinUniqueIDs,"data.frame")){
    proteinUniqueIDs <- bit64::as.integer64(proteinUniqueIDs$UniqueSequenceID)
  } else {
    proteinUniqueIDs <- bit64::as.integer64(proteinUniqueIDs)
  }
  return(dbGetTable(
    db = db,
    tableName = "TargetProteinsModificationSites",
    columnNames = "ModificationSitesId",
    filtering = paste(c(" WHERE TargetProteinsUniqueSequenceID IN (",
                        paste(c("'",
                                paste(proteinUniqueIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

#' function to get data from the ModificationSides table using the
#'  modificiationSiteId's
#'
#' @param db database access 'handle'
#' @param modificatonSitesIDs the modification site identifiers to get from
#'  the ModificationSites table
#' @param columnNames allows the selection of columns to take from the table,
#'  default = NA (all columns)
#' @param sortOrder allows for sorting of the selected columns,
#'  default = NA, (no sorting). Other valid values are a single character
#'  string ("ASC" or "DESC") or a character vector of the same length as the
#'  columnNames vector containing a series of "ASC" or "DESC"
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame
#'
#' @return a data.frame or a character vector (SQL)
#'  
#' @note the easiest way to get the modificationSitesIDs is via the
#'  \code{\link{dbGetModificationsSitesIDs}} function
#' @export
dbGetModificationsTable <- function(db,
                                    modificatonSitesIDs,
                                    columnNames = NA,
                                    sortOrder = NA,
                                    SQL = FALSE){
  if (is.Class(modificatonSitesIDs,"data.frame")){
    modificatonSitesIDs <- as.character(modificatonSitesIDs$ModificationSitesId)
  } else {
    if (!is.character(modificatonSitesIDs)){
      modificatonSitesIDs <- as.character(modificatonSitesIDs)
    }
  }
  dbGetTable(
    db = db,
    tableName = "ModificationSites",
    columnNames = columnNames,
    filtering = paste(
      c(" WHERE Id IN ",
        "(", paste(modificatonSitesIDs, collapse = ","),")"), collapse = ""),
    sortOrder = sortOrder,
    SQL = SQL)
}

#' Function to get the peptideID's 'belonging' to a modification site
#'
#' @param db database access 'handle'
#' @param modificationIDs the modification site identifiers to get from
#'  the ModificationSites table. This should be the 'Id' field of a modifciation
#'  table row
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame#'
#'  
#' @return a data.frame or a character vector (SQL)
#' 
#' @export
dbGetModificationPeptideIDs <- function(db, modificationIDs, SQL = FALSE){
  if (is.Class(modificationIDs,"data.frame")){
    modificationIDs <- modificationIDs$Id
  } 
  return(dbGetTable(
    db = db,
    tableName = "TargetPeptideGroupsModificationSites",
    columnNames = "TargetPeptideGroupsPeptideGroupID",
    filtering = paste(c(" WHERE ModificationSitesId IN (",
                        paste(c("'",
                                paste(modificationIDs,collapse = "','"),
                                "'"),
                              collapse = ""),
                        ")"),
                      collapse = ""),
    sortOrder = NA,
    SQL = SQL))
}

