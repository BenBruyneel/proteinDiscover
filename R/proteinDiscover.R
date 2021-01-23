library(DBI)
library(RSQLite)
library(dplyr)
library(stringr)
library(pool)
library(gtools)
library(ggplot2)

source("../proteinDrawing/proteinDrawing.R")

# ---- raw vector (blob) related functions ----

# function that converts a raw vector into its numeric counterpart(s)
# note: numeric vectors are 9 bytes long, the first 8 are the actual
#       number, if the last byte (9) == 0 then the result is NA 
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

# function that converts a raw vector into its integer counterpart(s)
# note: integer vectors are 5 bytes long, the first 4 are the actual
#       number, if the last byte (5) == 0 then the result is NA 
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

# gives boolean columns
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

# specials are not numeric or integer, but have chunks of size 
specials <- data.frame(names = c("AspectBiologicalProcess","AspectCellularComponent","AspectMolecularFunction"),
                       size = c(2,2,2))

# forceBlob must be data.frame with 3 columns: name (columnName), what (type), minimumSize (number of values in a cell)

# function that converts a data.frame column with raw vectors into
# one or more columns containing the integer/numeric counterparts
# columnVector = the data.frame column, eg df[,1] or df$column1
# what         = "integer" or "numeric" (can only be one of these two!)
# columnName   = character string, can be new name, can be the original
#                column's name. Must be suitable name.
# Note: if a rawVector column contains multiple values per 'cell' then
#       the column will get split into an equal number of columns with
#       the name 'columnName'+"_"+number, eg column1 becomes:
#       column1_1, column1_2, etc
convertRawColumn <- function(columnVector, what, columnName, minimumSize = 1, forceBlob = NA){
  if (columnName %in% specials$names){
    specialSize <- specials[specials$names == columnName,]$size
    converted <- unlist(lapply(columnVector, function(x){convertRawSpecial(x,specialSize, minimumSize = minimumSize)}))
  } else {
    if (!identical(forceBlob,NA)){
      if (columnName %in% forceBlob$name){
        what <- forceBlob[forceBlob$name == columnName,]$what
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

# A seperate function to determine what type of raw vector 'rawVector" is.
# It can return a number (length of the bytes) or the type (numeric or integer)
# If not one of those two, then NA will be returned
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

# Some of the data in the tables that come out of the database system in
# proteome discoverer are in the <blob> raw Vector format. So some columns
# are lists of raw Vectors. To convert them, one has to know the type of
# raw vector (integer/numeric). Because some elements in the 'blob' list
# are NA, this function seeks the first non-NA element and tries to determine
# it's type. If it the column contains only NA's then the naValue will be
# returned (default = NA). This was included to force a certain type
# even if the column (temporarily?) contains no data
whichRawList <- function(blobList, naValue = NA){
  anyValidBlobs <- sum(!is.na(blobList))
  if (anyValidBlobs < 1){
    return(naValue)
  } else {
    firstBlob <- blobList[[which(!is.na(blobList))[1]]]
    theType <- whichRaw(firstBlob)
    numberValues <- ifelse(theType == "integer",length(firstBlob)/5,
                           ifelse(theType == "numeric", length(firstBlob)/9, length(firstBlob)))  # unknowns will give length of the blob itself
    return(list(type = theType, number = numberValues))
  }
}

whichRawListSpecial <- function(blobList, naValue = NA, specialName){
  anyValidBlobs <- sum(!is.na(blobList))
  if (anyValidBlobs < 1){
    return(naValue)
  } else {
    firstBlob <- blobList[[which(!is.na(blobList))[1]]]
    theType <- "special"
    numberValues <- length(firstBlob) / specials[specials$names == specialName,]$size  # unknowns will give length of the blob itself
    return(list(type = theType, number = numberValues))
  }
}

# function determines the raw types (numeric or integer) of the columns in the data.frame
# note: if a column is not a raw vector column then its type will be NA
determineRawTypes <- function(df){
  rawType = list()
  for (counter in 1:(ncol(df))){
    blobList <- df[,counter][[1]]  # as.list(df[,counter])[1] --> doesn't work when first item = NA ?!
    if (colnames(df[,counter]) %in% specials$names){ # special
      rawType[[counter]] <- whichRawListSpecial(blobList, specialName = colnames(df[,counter]))
    } else {
      rawType[[counter]] <- whichRawList(blobList)
    }
  }
  return(rawType)
}

# the tables/data.frame's coming from Proteome Discoverer have columns of the type
# raw vecotr (blob). These can be converted automatically or semi-automatically by this function
# Note: can only do integer & numeric blobs atm
# Note: some raw vector columns are actually two (or possibly more) columns in one. In those cases
# each element/cell of the column is two (or more) values. This function splits these columns into
# two seperate ones.
# df = data.frame coming from a table from a Proteome Discoverer database (.pdResult)
# If there are no raw vector columns, then this function has no use and may even trigger errors/warnings
df_transform_raws <- function(df, forceBlob = NA){
  # figure out which columns are raw vector (class 'blob')
  blobColumns <- which(unlist(lapply(unname(lapply(df, function(x){class(x)})),function(x){return("blob" %in% x)})))
  # determine the type of each blob column
  colClass <- determineRawTypes(df[blobColumns])
  # start a new data.frame w/o the blob columns
  newdf <- df[,-blobColumns]
  blobNames <- names(df[,blobColumns])
  for (counter in seq_along(blobColumns)){
    if (!identical(colClass[counter],NA)){
      newColumns <- convertRawColumn(df[,blobColumns[counter]][[1]],what = colClass[[counter]]$type, columnName = blobNames[counter], minimumSize = colClass[[counter]]$number, forceBlob = forceBlob)
      newdf <- bind_cols(newdf, newColumns)
    } else {
      warningMessage <- paste(c("Warning: cannot automatically convert column '",blobNames[counter],"' "), collapse = "")
      warning(warningMessage)
      newdf <- bind_cols(newdf, df[,blobColumns[counter]])
    }
  }
  return(newdf)
}

# ---- database (related) functions ----

# proteome discoverer .pdResult files are in the SQLite format
# functions in this section deal with the database. Some of the
# functions work with database (query) commands, which is the preferred
# way.
# note: for future developments when buidling a shiny app it's recommended
# to use the 'pool' & 'dbplyr' packages to streamline things

# function to open a databse
# note: if no file with the name 'fileName' exists, then it will be created
# (but obviously it will be empty, so most further commands will fail)
db_open <- function(fileName){
  return(dbPool(
    drv = RSQLite::SQLite(),
    dbname = fileName
    )
  )
  #return(dbConnect(RSQLite::SQLite(), fileName))
}

# function to close a database
db_close <- function(db){
  poolClose(db)
}
# ---- database functions, all need an active database (db) (connection open!) ----

# function to get a table from the database, tabl = tablename (character)
# note: what comes out is NOT a data.frame, but a dynamic 'connection' of
# sorts that only works as long as the database is 'open'. To make it into
# a (static) database: either use as.data.frame() or collect()
db_colNames <- function(db,tabl){
  return(dbListFields(db,tabl))
}

# perform a get query command on a table
# tabl = tablename (character)
# this command is essentially the same as: 'SELECT qcmd FROM tabl' 
db_getQ <- function(db, qcmd , tabl){
  return(dbGetQuery(db,paste(c(qcmd," ",tabl), collapse = "")))
}

# find the number of rows in a table of an (open) database
# tabl = tablename (character)
db_nrow <- function(db,tabl){
  return(dbGetQuery(db,paste(c("SELECT COUNT(*) FROM ",tabl), collapse = ""))[1,1])
}

# gather info on all tables in the database
# functionn returns a table with the following columns:
# tablename, number of rows in each table and the fields in each table
db_tbl_def <- function(db){
  datazz <- dbListTables(mydb)
  dftbl <- data.frame(name = datazz,
                      nrows = unlist(lapply(datazz, function(x){db_nrow(mydb,x)})),
                      colnames = NA)
  for (counter in 1: length(datazz)){
    dftbl$colnames[counter] <- list(db_colNames(mydb,datazz[counter]))
  }
  return(dftbl)
}

# function returns all tables in the database in the tbl() format
# these are NOT data.frame's themselves. as.data.frame()/collect()
# is needed for that
db_get_tbls <- function(db){
  temp <- lapply(dbListTables(db), function(x){tbl(db, x)})
  names(temp) <- dbListTables(db)
  return(temp)
}

# function to get a single table from a database in tbl() format
# Use as.data.frame()/collect() to convert to a regular R data.frame
db_table <- function(db, tableName){
  return(tbl(db,tableName))
}

# ---- Proteome Discoverer specific code ----

# Function to transform the table target proteins into a proper data.frame
# note: the argument dbtble must be the in tbl() format. Under normal circumstances
# it is the "TargetProteins" table in the .pdResult file
# using IsMasterProtein = TRUE will give the table that proteome discoverer automatically
# displays (if filters are default)
db_get_proteins <- function(dbtbl, MasterProtein = TRUE){
  if (MasterProtein){
    return(df_transform_raws(dbtbl %>% filter(IsMasterProtein == 0) %>% collect()))
  } else {
    return(df_transform_raws(dbtbl %>% collect())) 
  }
}

# function to get the table "TargetProteins" from the database db and to transform it to
# data.frame properly
db_get_proteinTable <- function(db, MasterProtein = TRUE, sorting = NA, descending = TRUE){
  if (!identical(sorting, NA)){
    if (descending){
      return(db_get_proteins(tbl(db,"TargetProteins"), MasterProtein = MasterProtein) %>% arrange(desc(!!sym(sorting))))
    } else {
      return(db_get_proteins(tbl(db,"TargetProteins"), MasterProtein = MasterProtein) %>% arrange(!!sym(sorting)))
    }
  } else {
    return(db_get_proteins(tbl(db,"TargetProteins"), MasterProtein = MasterProtein))
  }
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

# leave out for now
# # to get the peptidepsmlink ("TargetPsmsTargetProteins")
# # in goes "TargetProteinsUniqueSequenceID" from the table "TargetProteins" (Note: it's possible to use a c(,,,)
# # result ia a list of numbers which are the "PeptideID" in the "TargetPsms" table
# db_get_psmIDsProtein <- function(db, ProteinUniqueSequenceIDs){
#   return((tbl(db, "TargetPsmsTargetProteins") %>% filter(TargetProteinsUniqueSequenceID %in% ProteinUniqueSequenceIDs) %>% collect())$TargetPsmsPeptideID)
# }

proteinIDTypes <- function(db){
  return((tbl(db,"ProteinIdentificationGroups") %>% collect())$GroupName)
}

MSfileInfo <- function(db, type = "XcaliburRawfile"){
  if (nchar(type) == 0){
    return(tbl(db,"WorkFlowInputFiles") %>% collect())
  } else {
    return(tbl(db,"WorkFlowInputFiles") %>% collect() %>% filter(FileType %in% type))
  }
}

MSfile <- function(db, clean = TRUE, type = "XcaliburRawfile"){
  if (!clean){
    return(MSfileInfo(db, type = type)$FileName)
  } else {
    return(gsub(MSfileInfo(db, type = type)$FileName,pattern = "\\\\", replacement ="/"))
  }
}

MSfile.DateTime <- function(db, type = "XcaliburRawfile"){
  return(MSfileInfo(db, type = type)$CreationDate)
}

MSfile.Date <- function(db, type = "XcaliburRawfile"){
  return(unlist(lapply(MSfileInfo(db, type = type)$CreationDate, function(x){strsplit(x, split = " ")[[1]][1]})))
}

MSfile.Time <- function(db, type = "XcaliburRawfile"){
  return(unlist(lapply(MSfileInfo(db, type = type)$CreationDate, function(x){strsplit(x, split = " ")[[1]][2]})))
}

MSfile.rtRange <- function(db, type = "XcaliburRawfile"){
  templ <- (strsplit(MSfileInfo(db, type = type)$RetentionTimeRange, split = " - "))
  tempdf <- data.frame(start = as.numeric(), end = as.numeric())
  for (counter in 1:(length(templ))){
    tempdf <- bind_rows(tempdf, data.frame(start = as.numeric(templ[[counter]][1]), end = as.numeric(templ[[counter]][2])))
  }
  return(tempdf)
}

MSfile.dbInfo <- function(db, type = "XcaliburRawfile"){
  return(MSfileInfo(db, type = type) %>%
           select(FileID,StudyFileID,FileName, CreationDate, RetentionTimeRange) %>%
           rename(StudyID = StudyFileID,Date = CreationDate, rtRange = RetentionTimeRange))
}

pd.searchInfo <- function(db){
  return(tbl(db,"WorkflowMessages") %>% collect())
}

# in seconds
pd.searchTime <- function(db){
  temptbl <- searchInfo(db)
  (temptbl$Time[nrow(temptbl)] - temptbl$Time[1])/10000000
}

# generates a (sort of) chromatogram using data in the LcmsPeaks:
# ApexRt (retentiontime at the apex of peak found (probably Minora))
# PeakHeight or PeakArea (Minora?)
pd.Chromatogram <- function(db, useArea = FALSE){
  tempdf <- tbl(db,"LcmsPeaks") %>% collect()
  if (!useArea){
    g <- ggplot(tempdf,aes(x = ApexRT, y = PeakHeight))
  } else {
    g <- ggplot(tempdf,aes(x = ApexRT, y = PeakArea))
  }
  g <- g + geom_line()
  g <- g + ggtitle("LCMS") + xlab("Rt (min)") + ylab("Intensity")
  g <- g + theme_classic()
  return(g)
}

pd.Chromatogram.limits <- function(db, useArea = FALSE, xlim = NA, ylim = NA){
  g <- pd.Chromatogram(db, useArea = useArea)
  if (identical(xlim,NA)){
    g <- g + scale_x_continuous(expand = c(0,0))
  } else {
    g <- g + scale_x_continuous(limits = xlim, expand = c(0,0))
  }
  if (identical(ylim,NA)){
    g <- g + scale_y_continuous(expand = c(0,0))
  } else {
    g <- g + scale_y_continuous(limits = ylim, expand = c(0,0))
  }
  return(g)
}

# for translation of the value (1..5) in the psmTable to words (like in Proteome Discoverer)
# note: <...> --> means not encountered/undefined/inference
PsmAmbiguity <- c("<Not Considered>","Rejected","<Ambiguous>","Selected","Unambiguous")

# for translation of the value (0..4) in the proteinTable to words (like in Proteome Discoverer)
# note: <...> --> means not encountered/undefined/inference
# note: range of IsMasterProtein = 0 - 4
IsMasterProtein <- c("Master Protein","Master Protein Candidate","<Undefined>","Master Protein Rejected","Master Protein Rejected")

# ---- proteinDrawing ----

# function to create a 4-column table with modification data from a peptide table
# the columns are: sequence, modification, aminoacid, position (within the peptide)
# there are occasions where the exact position cannot be determined, eg in case
# of unclear position of the modification (multiple possiblities). In those cases
# the position of the modification will become NA and the amino acid column
# may contain something like "K/C". Also this function will generate a warning
# message. These can be suppressed via suppressWarnings()
# This function should always be used together with adjustModification(), as
# that function will attempt to 'correct' some NA's and such. This may not work
# properly (it's difficult to foresee all things...)
translateModification <- function(peptideTable, duplicates = FALSE){
  moddf <- data.frame(sequence = as.character(), modification = as.character(), aa = list(), pos = list())
  allsplit <- strsplit(peptideTable$Modifications,split = " ")
  # remove all "[","[" and ";"
  allsplit <- lapply(allsplit,function(x){str_replace_all(x,pattern = "\\[", replace = "")})
  allsplit <- lapply(allsplit,function(x){str_replace_all(x,pattern = "\\]", replace = "")})
  allsplit <- lapply(allsplit,function(x){str_replace_all(x,pattern = ";", replace = "")})
  for (counter in 1:length(allsplit)){
    counter2 <- 1
    while (counter2 <= length(allsplit[[counter]])){
      # first piece = modification, the rest are position pieces
      currentMod <- allsplit[[counter]][counter2]
      currentaa <- as.character()
      currentpos <- as.integer()
      counter2 <- counter2 + 1
      while ((!grepl(chr(215),allsplit[[counter]][counter2])) & (counter2 <= length(allsplit[[counter]]))) { # 215 is ascii of '×' which is NOT x
        currentaa <- append(currentaa,str_extract(allsplit[[counter]][counter2],"\\D*"))
        currentpos <- append(currentpos,as.integer(str_extract(allsplit[[counter]][counter2],"\\d*$")))
        counter2 <- counter2 + 1
      }
      moddf <- bind_rows(moddf, data.frame(sequence = peptideTable$Sequence[counter],modification = currentMod, aa = currentaa, pos = currentpos, stringsAsFactors = FALSE))
    }
  }
  moddf$modification <- gsub(pattern = paste("\\d+",chr(215), sep = ""),replace = "",moddf$modification)
  if (!duplicates){
    moddf <- moddf %>% distinct(.keep_all = TRUE)
  }
  return(moddf)
}

# separate function to allow for completely manual adjustments/flexibility
# This function attempts to 'correct' some minor things in the tables generated
# by translateModification(). Eg N-terminal modifications' position.
# It's intedend to be extended/rewritten as more 'solvable' problems are
# found. Currently it's abilities are limited.
adjustModification <- function(translatedTable){
  nas <- which(is.na(translatedTable$pos))
  for (counter in 1:(length(nas))){
    if (translatedTable$aa[nas[counter]] == "N-Term"){
      translatedTable$pos[nas[counter]] <- 1
    } else {
      if ((translatedTable$aa[nas[counter]] == "C-Term")){  # note: not seen yet!
        translatedTable$pos[nas[counter]] <- nchar(translatedTable$sequence[nas[counter]])
      } else {
        countaas <- str_count(translatedTable$sequence[nas[counter]], pattern = translatedTable$aa[nas[counter]])
        if (countaas == 1){
          str_locate(translatedTable$sequence[nas[counter]],pattern = translatedTable$aa[nas[counter]])[1,1]
        } else {
          warningMessage <- paste(c("Cannot determine position of modification: ", translatedTable$modification[nas[counter]],
                                    " (",translatedTable$aa[nas[counter]],") in row: ",toString(nas[counter])), collapse = "")
          warning(warningMessage)
        }
      }
    }
  }
  return(translatedTable)
}

# function that creates an experimentTable compatible with the experimentChain class in proteinDrawing
# db = proteome discoverer result database (.pdResult, must be connected obviously)
# foundPeptides = table generated by db_get_peptideTable()
# foundProtein = a number! should be the same proteinGroupID as used to generate the peptideTable
# description = character string, usually the name of the protein
# modificationTable = table, generated via translateModification() %>% adjustModification() (adjustModification() is optional)
# ExperimentOrder = number of the protein chain when drawing, see proteinDrawing.R
# seekPSM = when tesitng it was found that proteome discoverer sometimes generates peptide sequences, that are NOT in the
#           protein, in it's tables. This is due to difficulties like eg distinguishing between 'ILP' and 'LLP' sequences
#           if seekPSM == TRUE, then this function will attempt to search the PSM table of the specific protein/peptide
#           table for an alternative. Note if a 'wrong' sequence is found a message will be generated and when an
#           alternative is found, an additional message is generated. These can be suppressed via suppressMessages()
# extraColumns = here extra columns from the peptideTable can be defined to be added to the experimentTable, eg 'Abundances'
generateExperimentTable <- function(db, foundPeptides, foundProtein, description, modificationTable, ExperimentOrder = 1, seekPSM = TRUE, extraColumns = NA){
  # first create COVER part
  # get sequence of the protein
  foundProtein <- (db_get_proteins(db_table(db,"TargetProteins")) %>% filter(ProteinGroupIDs == foundProtein))$Sequence[1]
  df <- data.frame()
  for (counter in 1:(nrow(foundPeptides))){
    positions <- str_locate_all(foundProtein, pattern = foundPeptides$Sequence[counter])
    if (nrow(positions[[1]]) == 0){
      messageText <- paste(c("The following sequence is not in protein sequence: ",foundPeptides$Sequence[counter]),collapse = "")
      message(messageText)
      if (seekPSM){
        psmList <- db_get_psmIDs(db,foundPeptides$PeptideGroupID[counter])
        psmList <- db_get_psmTable(db,psmList)
        for (counterx in 1:(nrow(psmList))){
          positions <- str_locate_all(foundProtein, pattern = psmList$Sequence[counterx])
          if (nrow(positions[[1]]) != 0){
            foundPeptides$Sequence[counter] <- psmList$Sequence[counterx]  # replace sequence in peptide table
            messageText <- paste(c("The peptide was replaced by: ",foundPeptides$Sequence[counter]),collapse = "")
            message(messageText)
            break                                                          # exit loop since a match is found
          }
        }
      }
    }
    # check again if there's something in positions
    if (nrow(positions[[1]]) != 0){
      for (counter2 in 1:(nrow(positions[[1]]))){
        tempdf <- data.frame(
          type = "COVER",
          description = description,
          begin = positions[[1]][[counter2,1]],
          end = positions[[1]][[counter2,2]],
          length = nchar(foundPeptides$Sequence[counter]),
          order = ExperimentOrder,
          sequence = foundPeptides$Sequence[counter],
          mods = foundPeptides$Modifications[counter],
          expMr = foundPeptides$mzDabySearchEngine[counter],
          expCharge = foundPeptides$ChargebySearchEngine[counter],
          calcMr = foundPeptides$TheoreticalMass[counter],
          rt = foundPeptides$RTminbySearchEngine[counter],
          XCorr = foundPeptides$XCorrbySearchEngine[counter],
          stringsAsFactors = FALSE # for compatibility with R versions < 4
        )
        if (!identical(extraColumns,NA)) {
          tempdf <- bind_cols(tempdf, foundPeptides[counter,] %>% select(all_of(extraColumns)))
        }
        df <- bind_rows(df,tempdf)
      }
    }
  }
  df <- df %>% arrange(begin)
  tempdf <- data.frame(
    type = "CHAIN",
    description = description,
    begin = 1,
    end = nchar(foundProtein),
    length = nchar(foundProtein),
    order = ExperimentOrder,
    sequence = foundProtein,
    mods = NA,
    expMr = NA,
    expCharge = NA,
    calcMr = NA,
    rt = NA,
    XCorr = NA,
    stringsAsFactors = FALSE # for compatibility with R versions < 4
  )
  df <- bind_rows(tempdf,df)
  # add modificationTable
  tempdf <- data.frame()
  for (counterM in 1:(nrow(modificationTable))){
    # search again for the position of the sequences in modification table
    positions <- str_locate_all(foundProtein, pattern = modificationTable$sequence[counterM])
    if (nrow(positions[[1]]) == 0){
      # no message in this case, assumption is made that this is already covered in the "coverage"
      if (seekPSM){
        psmList <- db_get_psmIDs(db,foundPeptides$PeptideGroupID[counter])
        psmList <- db_get_psmTable(db,psmList)
        for (counterx in 1:(nrow(psmList))){
          positions <- str_locate_all(foundProtein, pattern = psmList$Sequence[counterx])
          if (nrow(positions[[1]]) != 0){
            foundPeptides$Sequence[counter] <- psmList$Sequence[counterx]  # replace sequence in peptide table
            # no message in this case, see before
            break                                                          # exit loop since a match is found
          }
        }
      }
    }
    # check again if there's something in positions
    if (nrow(positions[[1]]) != 0){
      for (counter2 in 1:(nrow(positions[[1]]))){
        tempdf <- bind_rows(tempdf,
                            data.frame(
                              type = "MOD_RES",
                              description = modificationTable$modification[counterM],
                              begin = positions[[1]][[counter2,1]]+modificationTable$pos[counterM]-1,
                              end = positions[[1]][[counter2,1]]+modificationTable$pos[counterM]-1,
                              length = 1,
                              order = ExperimentOrder,
                              sequence = modificationTable$aa[counterM],
                              mods = NA,
                              expMr = NA,
                              expCharge = NA,
                              calcMr = NA,
                              rt = NA,
                              XCorr = NA,
                              stringsAsFactors = FALSE # for compatibility with R versions < 4
                            )
        )
      }
    }
  }
  tempdf <- tempdf %>% distinct(.keep_all = TRUE) %>% arrange(begin)
  df <- bind_rows(df, tempdf)
  return(df)
}