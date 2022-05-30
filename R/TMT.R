# ---- general ----

#' helper function to calculate a row-wise function (like mean, median etc)
#'  across a data.frame
#'  
#' @param data the data.frame. Note that all rows and columns are used, so
#'  selection, filtering, etc should be done beforehand
#' @param setNAZero default = NA, when NA this is ignored. Otherwise all cells
#'  containing NA will be set to the value of setNAZero. When removeNAs = TRUE,
#'  this parameter is ignored
#' @param removeNAs default = FALSE, if TRUE all rows containing NA's will be
#'  removed via na.omit()
#' @param keepData if TRUE, then the original data is returned also
#' @param calcName name of the column with the calculated values in it
#' @param calcFunc function to be applied row-wise across the data.frame
#' @param ... serves to pass on "extra" arguments on to the calcFunc function,
#'  eg na.rm = TRUE in case of calcFunc = mean
#'
#' @returns a data.frame with the calculated values as the only column or with
#'  the calculated values as a mew column
#' @export
calcData <- function(data,
                     setNAZero = NA, removeNAs = FALSE, keepData = FALSE,
                     calcName = "median", calcFunc = stats::median, ...){
  if (removeNAs){
    data <- data %>% stats::na.omit()
  } else {
    if (!is.na(setNAZero)){
      data[is.na(data)] <- setNAZero 
    }
  }
  data$calculated <- NA
  for (counter in 1:nrow(data)){
    data$calculated[counter] <- calcFunc(unlist(data %>%
                                                  dplyr::slice(counter) %>%
                                   dplyr::select(-dplyr::all_of("calculated"))),
                                         ...)
  }
  if (keepData) {
    colnames(data)[(colnames(data) == "calculated")] <- calcName
  } else {
    data <- data %>%
      dplyr::select("calculated")
    colnames(data) <- calcName
  }
  return(data)
}

#' helper function to generate the a data.frame of proteins info for other
#'  functions
#'  
#' @returns a data.frame with three columns: short (character vector),
#'  Accession (character vector, uniprot "style") and knockout (logical)
#'  
#' @note even though it's called knockOutProteins, 2 of the proteins are
#'  not knock out proteins.
#' 
#' @export
knockOutProteins <- function(){
  return(data.frame(
    short = c("His4","Eno2","Met6","Ura2","PABP_YEAST"),
    Accession = c("P00815","P00925","P05694","P07259","P04147"),
    knockout = c(TRUE,FALSE,TRUE,TRUE,FALSE)) %>% dplyr::arrange("Accession"))
}

#' helper function to generate the a data.frame of TMT knockout strain (TKO)
#'  info for other functions. This function generates a data.frame based on
#'  the 11-plex TMT TKO knockout
#'  
#' @note the rows define the order of the abundance (etc) columns in the
#'  protein, peptide and psms table in a pdResult file. The order is
#'  alphabetical in protein & peptide tables, but not in the psms tables: there
#'  it is based based on the order of the isotopes
#'  
#' @note psmsChannels & isotopeChannels columns match each other
#' 
#' @returns a data.frame with four columns: all are character vectors
#' @export
tmt11Channels <- function(){
  return(data.frame(
    proteinChannels = c(rep("His4",3),
                        rep("Met6",3),
                        rep("Parental",2),
                        rep("Ura2",3)),
    peptideChannels = c(rep("His4",3),
                        rep("Met6",3),
                        rep("Parental",2),
                        rep("Ura2",3)),
    psmsChannels = c(rep("Met6",3),
                     rep("His4",3),
                     rep("Ura2",3),
                     rep("Parental",2)),
    isotopeChannels = c("126","127N","127C","128N","128C","129N",
                        "129C","130N","130C","131N","131C")))
}

#' helper function to generate the a data.frame of TMT knockout strain (TKO)
#'  info for other functions. This function generates a data.frame based on
#'  the 10-plex TMT TKO knockout (this was the original TMT-knockout-digest
#'  available)
#'  
#' @note the rows define the order of the abundance (etc) columns in the
#'  protein, peptide and psms table in a pdResult file. The order is
#'  alphabetical in protein & peptide tables, but not in the psms tables: there
#'  it is based based on the order of the isotopes
#'  
#' @note psmsChannels & isotopeChannels columns match each other
#'  
#' @returns a data.frame with four columns: all are character vectors
#' @export
tmt10Channels <- function(){
  return(data.frame(
    proteinChannels = c(rep("His4",3),
                        rep("Met6",3),
                        rep("Ura2",3),
                        rep("Wildtype",1)),
    peptideChannels = c(rep("His4",3),
                        rep("Met6",3),
                        rep("Ura2",3),
                        rep("Wildtype",1)),
    psmsChannels = c(rep("Met6",3),
                     rep("His4",3),
                     rep("Ura2",3),
                     rep("Wildtype",1)),
    isotopeChannels = c("126","127N","127C","128N","128C",
                        "129N","129C","130N","130C","131")))
}

# ---- info functions (proteins/peptides) ----

#' get protein info (without translation of columns) from a list of protein
#'  Accessions (uniprot code). Essentially this is a wrapper function for
#'  \code{\link[proteinDiscover]{dbGetTable}}
#'  
#' @param db database access 'handle'
#' @param columnNames allows the selection of columns to take from the table
#' @param proteinAccessions defines which protein(s) info will be retrieved
#'  (character vector)
#' @param sortOrder allows for sorting of the resulting data.frame by on of it's
#'  columns (default = "Accession")
#' @param SQL allows the function to return the SQL query statement in stead of
#'  a data.frame (for debugging purposes)
#'  
#' @return a data.frame containing requested data from the protein table or
#'  a character string specifying an SQL query
#' @export
getProteinInfoRaw <- function(db,
                              columnNames = c("Accession",
                                              "ProteinGroupIDs",
                                              "AbundancesNormalized",
                                              "AbundanceRatios",
                                              "AbundanceRatioPValue",
                                              "AbundanceRatioAdjPValue"),
                              proteinAccessions = knockOutProteins()$Accession,
                              sortOrder = "Accession",
                              SQL = FALSE){
  return(dbGetTable(db = db,
                    tableName = tableNames("proteins"),
                    columnNames = columnNames,
           filtering = paste(c(" WHERE IsMasterProtein = 0 AND Accession IN ('",
                                    paste(proteinAccessions, collapse = "', '"),
                                        "')"),
                                      collapse = ""),
                    sortOrder = sortOrder,
                    SQL = SQL))
}

#' get protein info (with translation of columns) from a list of protein
#'  Accessions (uniprot code). Essentially this is a wrapper function for
#'  \code{\link{getProteinInfoRaw}}
#'  
#' @param db database access 'handle'
#' @param columnNames allows the selection of columns to take from the table
#' @param proteinAccessions defines which protein(s) info will be retrieved
#'  (character vector)
#' @param sortOrder allows for sorting of the resulting data.frame by on of it's
#'  columns (default = "Accession")
#'
#' @note this function uses the default
#'  \code{\link[proteinDiscover]{getProteinInfoRaw}} function. If more control
#'  over the "translation" of raw columns is needed, then use
#'  \code{\link{getProteinInfoRaw}} and do the translation manually
#'
#' @return a data.frame containing requested data from the protein table after
#'  "translation" of the raw columns
#' @export
getProteinInfo <- function(db,
                           columnNames = c("Accession",
                                           "ProteinGroupIDs",
                                           "AbundancesNormalized",
                                           "AbundanceRatios",
                                           "AbundanceRatioPValue",
                                           "AbundanceRatioAdjPValue"),
                           proteinAccessions = knockOutProteins()$Accession,
                           sortOrder = "Accession"){
  return(getProteinInfoRaw(db = db,
                           columnNames = columnNames,
                           proteinAccessions = proteinAccessions,
                           sortOrder = sortOrder) %>%
           dfTransformRaws())
}

#' get peptide information from the peptide table from a pdResult file based on
#'  the provided proteinAccession (uniprot) codes. Raw columns are not
#'  "translated"
#'  
#' @param db database access 'handle'
#' @param columnNames allows the selection of columns to take from the table. 
#'  The columns: PeptideGroupID, Sequence, Modifications, QuanInfo are
#'  automatically included. Default column to be retrieved is
#'  AbundancesNormalized
#' @param proteinAccessions defines from which protein(s) info will be retrieved
#'  (character vector)
#'  
#' @returns a named list of data.frames (the names are the proteinAccessions)
#' @export
getPeptideInfoRaw <- function(db,
                              columnNames = "AbundancesNormalized",
                              proteinAccessions = knockOutProteins()$Accession){
  ProteinGroupIDs <- unlist(lapply(proteinAccessions,
                            function(x){
                              getProteinInfoRaw(db = db,
                                                columnNames = c("Accession","ProteinGroupIDs"),
                                                proteinAccessions = x,
                                                sortOrder = NA)$ProteinGroupIDs
                            }))
  tempResult <- lapply(ProteinGroupIDs,
                   function(x){
                     dbGetPeptideIDs(db = db, proteinGroupIDs = x) %>%
                           dbGetPeptideTable(db = db,
                                        columnNames = append(c("PeptideGroupID",
                                                                    "Sequence",
                                                                "Modifications",
                                                                    "QuanInfo"),
                                                                  columnNames),
                                             sortOrder = c("Sequence",
                                                           "Modifications"))
                       })
  names(tempResult) <- proteinAccessions
  return(tempResult)
}

#' get peptide information from the peptide table from a pdResult file based on
#'  the provided proteinAccession (uniprot) codes. Raw columns are "translated"
#'  
#' @param db database access 'handle'
#' @param columnNames allows the selection of columns to take from the table. 
#'  The columns: PeptideGroupID, Sequence, Modifications, QuanInfo are
#'  automatically included. Default column to be retrieved is
#'  AbundancesNormalized
#' @param proteinAccessions defines from which protein(s) info will be retrieved
#'  (character vector)
#' @param removeUnusedQuantInfo default = TRUE. IF TRUE then only peptide info
#'  rows with NA as QuantInfo are kept (the others contain problematic abundance
#'  info or none at all)
#'  
#' @note this function uses the default
#'  \code{\link[proteinDiscover]{getProteinInfoRaw}} function. If more control
#'  over the "translation" of raw columns is needed, then use
#'  \code{\link{getPeptideInfoRaw}} and do the translation manually
#'  
#' @returns a named list of data.frames (the names are the proteinAccessions)
#' @export
getPeptideInfo <- function(db,
                           columnNames = "AbundancesNormalized",
                           proteinAccessions = knockOutProteins()$Accession,
                           removeUnusedQuantInfo = TRUE){
  tempResult <- getPeptideInfoRaw(db = db,
                                  columnNames = columnNames,
                                  proteinAccessions = proteinAccessions)
  tempResult <- lapply(tempResult, function(x){
    if (removeUnusedQuantInfo){
      x[is.na(x$QuanInfo),] %>% dfTransformRaws()
    } else {
      x %>%
        dfTransformRaws()
    }
  })
  names(tempResult) <- proteinAccessions
  return(tempResult)
}

# ---- Interference Free Index calculations ----

#' function to calculate the IFI (interference free index) of a protein
#'  entry in the protein table of a pdResult files. Note this can only be
#'  calculated on the knockout proteins in the TKO control sample: see
#'  \code{\link{tmt10Channels}} or \code{\link{tmt11Channels}} for the eligible
#'  proteins
#'  
#' @param db database access 'handle'
#' @param selected (short) name of the selected protein
#' @param accession uniprot accession code of the selected protein. If parameter
#'  "selected" is one of the short names in \code{\link{knockOutProteins}} then
#'  doesn't need to be specified. Note that the accession does not need to be
#'  one of the accessions of the knockout proteins
#' @param columns usually this will be "Abundances". It allows the selection
#'  of the correct (raw) columns as they come out of dfTransformRaws(), eg
#'  Abunances_1, Abundances_2, etc
#' @param groups usually either tmt10Channels() or tmt11Channels: data.frame
#'  that specifies which (abundance) column belongs to which knock out group.
#'  Note that the 'selected' argument should be in groups
#' @param IFIName specifies the name to give to the calculated values, usually
#'  "IFI"
#' @param calcFunc function to be applied row-wise across the data.frame. Used
#'  in the calculation of the IFI values. Default = mean
#' @param calcName name of the column with the calculated values in it, used
#'  in the related function calcData()
#' @param na.rm default = TRUE. This specifies that NA's should be removed when
#'  using eg mean, median, etc
#'
#' @returns a data.frame with two columns: one with the (short) name of the
#'  (selected) protein  and one with the calculated values (named IFI)
#' @export
calcIFIs <- function(db,
                     selected = "His4",
                     accession = 
             knockOutProteins()$Accession[knockOutProteins()$short == selected],
                     columns = "Abundances",
                     groups = tmt11Channels(),
                     IFIName = "IFI",
                     calcFunc = mean, calcName = "mean", na.rm = TRUE){
  data = getPeptideInfo(db = db, columnNames = columns)[[accession]] %>%
    dplyr::select(dplyr::starts_with(paste0(columns,"_")))
  tempGroups <- unique(groups$peptideChannels)
  tempResult <- lapply(tempGroups, function(x){
    calcData(data %>% dplyr::select(which(groups$peptideChannels == x)),
             calcFunc = calcFunc,
             calcName = x,
             na.rm = na.rm)})
  tempResult <- dplyr::bind_cols(tempResult)
  tempResult <- dplyr::bind_cols(tempResult %>%
                                   dplyr::select(which(tempGroups == selected)),
                            calcData(tempResult[,which(tempGroups != selected)],
                                     calcFunc = calcFunc, calcName = calcName,
                                     na.rm = na.rm))
  tempResult <- 1-(tempResult[,1]/tempResult[,2])
  tempResult <- data.frame(variable = selected, value = tempResult)
  colnames(tempResult) <- c("Short", IFIName)
  return(tempResult)
}

#' Wrapper function that uses \code{\link{tmt11Channels}} to calculate the
#'  IFI's for a set of (knock out) protein channels
#'
#' @param db database access 'handle'
#' @param proteinsKnockedOut character vector that specifies the (knock out)
#'  protein channels for which the IFI's are to be calculated
#' @param accession single element character vector specifying the accession
#'  of the protein whose abundances are to be used for the IFI calculation
#' @param groups usually either tmt10Channels() or tmt11Channels: data.frame
#'  that specifies which (abundance) column belongs to which knock out group
#' @param joined defines the type of output: if TRUE then a single data.frame
#'  with all IFI's for all (knock out) proteins is generated. Otherwise a list
#'  of data.frame's is generated for all proteins separately
#'
#' @return a data.frame with two columns: one with the (short) name of the
#'  (selected) proteins  and one with the calculated values (named IFI) or a
#'  list of data.frame's with the same structure
#' @export
calcAllIFIs <- function(db,
                        proteinsKnockedOut = 
                          knockOutProteins()$short[knockOutProteins()$knockout],
                        accession = NA,
                        groups = tmt11Channels(),
                        joined = TRUE){
  if (!identical(accession, NA)){
    result <- lapply(proteinsKnockedOut,
                     function(x){calcIFIs(db = db,
                                          selected = x,
                                          accession = accession,
                                          groups = groups)})
  } else {
    result <- lapply(proteinsKnockedOut,
                     function(x){calcIFIs(db = db,
                                          selected = x,
                                          groups = groups)})
  }
  if (joined){
    return(dplyr::bind_rows(result))
  } else {
    return(result)
  }
}