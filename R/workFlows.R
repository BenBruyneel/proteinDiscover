library(stringr)
library(XML)

#' function that generates the default data.frame for the function
#' df_replace().  
#'  
#' @return a data.frame with columns:
#'  value, replacement and singleChar
#'
#' @note value is the string to be searched, replacement is what it needs
#'  to be replaced with. singleChar sets whether the replacement should only
#'  take place when dealing with single character strings. This is because
#'  single character strings sometimes 'act' different when rendering
#'  markdown documents in HTML
#'  
#' @export
replacementStrings <- function(){
  return(
    data.frame(value = c("\\+","@"),
               replacement = c("$+$"," $@$ "),
               singleChar = c(TRUE, FALSE))
  )
} 

#' function that replaces (parts of) strings in a data.frame according to
#'  a provided table of replacements
#'  
#' @param df data.frame that needs to have strings replaced. Each cell is
#'  processed with str_replace_all from the stringr package for all elements of
#'  the str_replacements data.frame
#' @param str_replacements data.frame defining the replacements, see
#'  replacementStrings for more information
#'  
#' @return the data.frame with (parts of) strings replaced if present  
#' 
#' @note this function can be called just before passing a data.frame over
#'  to eg kableExtra::kbl(). When used in HTML markdown this function
#'  sometimes generates unintended behavior, eg converting (part of) strings
#'  to email addresses when they contain an @@ sign. This functions
#'  can replace possible problematic parts with something else. This can be eg
#'  latex. For example: replace '@@' with '$@@$' will solve the email address
#'  'problem'
#' @note for obvious reasons only character vector columns are processed
#'
#' @export
df_replace <- function(df, str_replacements = replacementStrings()){
  if (identical(df,NA)){
    return(NA)
  } else {
    if (!is.Class(df,"data.frame")){
      warning("df_replace only works on a data.frame, tibble or similar")
      return(df)
    } else {
      if (identical(str_replacements,NA)){
        return(df)
      }
    }
  }
  df <- as.data.frame(df)
  for (dfColumn in 1:ncol(df)){
    if (class(df[,dfColumn]) == "character"){
      for (rowCounter in 1:nrow(df)){
        for (replaceCounter in 1:nrow(str_replacements)){
          if (nchar(df[rowCounter,dfColumn]) == 1){
            df[rowCounter,dfColumn] <-
              stringr::str_replace_all(df[rowCounter,dfColumn],
                                       pattern = 
                                         str_replacements$value[replaceCounter],
                                       replacement = 
                                         str_replacements$replacement[replaceCounter])
          } else {
            if (!str_replacements$singleChar[replaceCounter]){
              df[rowCounter,dfColumn] <-
                stringr::str_replace_all(df[rowCounter,dfColumn],
                                         pattern = 
                                           str_replacements$value[replaceCounter],
                                         replacement = 
                                           str_replacements$replacement[replaceCounter])
            }
          }
        }
      }
    }
  }
  return(df)
}

#' function to get the workflow information from a .pdResult file
#' 
#' @param db database access 'handle' pointing to a .pdResult file
#' @param workflowsTable name of the table containing the info. Default is
#'  'WorkFlows'
#' @param returnNodeData if TRUE then the node parameters are included in the
#'  returned data
#' 
#' @return either a single data.frame containing basic info on the workflows or
#'  (if returnNodeData is TRUE) a list of the data.frame with the second list
#'  element containing information on the nodes that make up the
#'   processing & the consensus workflows (in xmlToList result
#'  format). This second element (called nodeInfo) is used in additional
#'  functions to show/display the processing/consensus workflows.
#'  
#' @export
workflowInfo <- function(db, workflowsTable = "WorkFlows",
                         returnNodeData = TRUE){
  sqlString <- paste(c("SELECT * FROM ",
                       workflowsTable),
                     collapse = "")
  resultTable <- pool::dbGetQuery(db,sqlString)
  resultXMLs <- lapply(resultTable$WorkflowXML, XML::xmlToList)
  names(resultXMLs) <- resultTable$WorkflowType
  resultTable$numberOfNodes <- unlist(lapply(resultXMLs,
                                             function(x){length(x[[1]])}))
  # get template names from .attrs
  resultTable$template <- unlist(lapply(resultXMLs,
                                        function(x){x[[2]]["TemplateName"]}))
  # info to get from resultTable (+ new names)
  workflowsTableNames = list(name = "WorkflowName", 
                             description = "WorkflowDescription", 
                             startDate = "WorkflowStartDate", 
                             user = "User", 
                             software = "SoftwareVersion", 
                             pc = "MachineName", 
                             type = "WorkflowType", 
                             study = "Study")
  resultTable <- resultTable %>%
    dplyr::select(dplyr::all_of(unlist(workflowsTableNames)),
                  .data$template,
                  .data$numberOfNodes)
  if (!returnNodeData){
    return(resultTable)
  } else {
    return(list(workflowInfo = resultTable,
                nodeInfo = resultXMLs))
  }
}

#' function to display an overview table of the processing/consensus workflows
#'  in the nodeInfo coming out of the workflowInfo function
#'  
#' @param nodeInfo either the processing or consensus part of the nodeInfo
#' 
#' @note an example of it's use: (workflowInfo(db))$nodeInfo$Consensus
#'  %>% nodeTable()
#' 
#' @export
nodeTable <- function(nodeInfo){
  # info to get from nodeInfo (+ new names)
  nodeTableNames = list(node = "ProcessingNodeNumber",
                        name = "FriendlyName",
                        description = "Description",
                        category = "Category",
                        parent = "ParentProcessingNodeNumber")
  numberNodes <- length(nodeInfo[[1]])
  nodeTable <- dplyr::bind_rows(lapply(nodeInfo[[1]], function(x){x[[".attrs"]]})) %>%
    dplyr::select(dplyr::all_of(unlist(nodeTableNames)))
  return(nodeTable)
}

#' function to create a DiagrammeR string that can be used by
#'  DiagrammeR::grViz() to plot a visual representation of the workflow
#'  
#' @param nodesTable output from the nodeTable function. Columns that
#'  need to be present are node, name & parent
#' @param showBelow boolean, default = TRUE. Set to FALSE when troubleshooting.
#'  Note that if set to FALSE, the parameter returnString will be ignored
#'  It is not recommended to depend on this parameter, as it will probably be
#'  removed in a newer version of the package
#' @param returnString default = TRUE. Set to FALSE when troubleshooting. Note
#'  that the parameter showBelow makes it so that this parameter is ignored.
#'  It is not recommended to depend on this parameter, as it will probably be
#'  removed in a newer version of the package
#' @param hideDoubleParents either NA (ignored) or a data.frame specifying what
#'  tp do in case of multiple parents. The data.frame should have the columns
#'  name and parent. The parent column should specify
#'  which parent to use ('first' or 'last') for connections
#'  
#' @return character vector that can be passed on to DiagrammeR::grViz()
#'  
#' @note during development it was noticed that some elements (nodes in the
#'  diagram) have more than one parent which is not seen in the proteome
#'  discoverer software of Thermo Scientific. The default data.frame 'corrects'
#'  known multiple parent nodes. If the parameter hideDoubleParents is
#'  set to NA, then the double parent connections are drawn.
#' 
#' @note an example of it's use:
#' (workflowInfo(db))$nodeInfo$Consensus %>%
#'  nodeTable() %>%
#'   createDiagrammeRString() %>%
#'    grViz()
#' 
#' @export
createDiagrammeRString <- function(nodesTable, showBelow = TRUE,
                                   returnString = TRUE,
                                   hideDoubleParents = data.frame(
                                     name = c("Precursor Ions Quantifier",
                                              "Feature Mapper",
                                              "Reporter Ions Quantifier",
                                              "Protein Marker",
                                              "Peptide in Protein Annotation",
                                              "Modification Sites",
                                              "Peptide Isoform Grouper"),
                                     parent = c("last","first","last","first",
                                                "first","last","first"))){
  nodesString <- ""  # node definitions to be placed in diagram
  nodesString2 <- "" # node connections
  nodesString3 <- "" # node labels
  nodesString4 <- "" # node definitions to be placed below diagram
  nodesString4c <- "" # node connections below diagram (fake, but)
  for (counter in 1:nrow(nodesTable)){
    newNodeString <- paste(c("node [label = '@@",
                             as.character(counter),
                             "'] n",
                             as.character(nodesTable$node[counter]),
                             " \n"), collapse = "")
    if ((nodesTable$parent[counter] == "-1") &
        (nodesTable$node[counter] != "0")){
      # not connected to a parent and not top of diagram
      nodesString4 <- paste0(nodesString4, newNodeString)
      if (nchar(nodesString4c) == 0){
        nodesString4c <- paste0("n", as.character(nodesTable$node[counter]))
      } else {
        nodesString4c <- paste(c(nodesString4c,
                                 " -> n",
                                 as.character(nodesTable$node[counter])),
                               collapse = "")
      }
    } else {
      nodesString <- paste0(nodesString, newNodeString)
    }
    if (nodesTable$parent[counter] != "-1"){ # has connection if != -1
      newNodesString <- ""
      parents <- strsplit(nodesTable$parent[counter],";")[[1]]
      if (!identical(hideDoubleParents,NA)){
        hideIt <-  hideDoubleParents %>%
          dplyr::filter(.data$name == nodesTable$name[counter])
        if (nrow(hideIt) > 0){
          if (hideIt$parent[1] == "first"){
            parents <- parents[1]
          } else {
            if (hideIt$parent[1] == "last"){
              parents <- parents[length(parents)]
            } else {
              parents <- parents[as.integer(hideIt$parent[1])]
            }
          }
        }
      }
      for (counter2 in 1:length(parents)){
        newNodeString <- paste(c("n",
                                 as.character(parents[counter2]),
                                 " -> n",
                                 as.character(nodesTable$node[counter]),
                                 "\n"),
                               collapse = "")
        nodesString2 <- paste0(nodesString2, newNodeString)
      }
    }
    newNodeString <- paste(c("[",
                             as.character(counter),
                             "]: '",
                             as.character(nodesTable$node[counter]),
                             ": ",
                             nodesTable$name[counter],
                             "' \n"),collapse = "")
    nodesString3 <- paste0(nodesString3, newNodeString)
    
  }
  if ((nchar(nodesString4) == 0) | ((nchar(nodesString4) != 0) & !showBelow)){
    result <- list(defineString = "digraph thegraph { \n node [fontname = Helvetica, shape = rectangle, fixedsize = true, width = 3]\n",
                   nodeDefinitions = nodesString,
                   nodeConnections = nodesString2,
                   endDefineString = "}\n",
                   nodeLabels =  nodesString3)
  } else {
    if (grepl(nodesString4c, pattern = "->")){
      # set alpha for connections to 0
      nodesString4c <- paste0(nodesString4c, " [color = '#00000000']")
    } else {
      # otherwise wipe connections (not needed)
      nodesString4c <- ""
    }
    result <- list(defineString = "digraph thegraph { \n rankdir = TB\n node [fontname = Helvetica, shape = rectangle, fixedsize = true, width = 3]\n",
                   clusterTopStart = "subgraph clustertop { peripheries = 0 \n",
                   nodeDefinitions = nodesString,
                   nodeConnections = nodesString2,
                   clusterTopEnd = "\n}\n",
                   clusterBelowStart = "subgraph clusterbelow { peripheries = 0 \n",
                   nodeBelowDefinitions = nodesString4,
                   nodeBelowConnections = nodesString4c,
                   clusterBelowEnd = "\n}\n",
                   endDefineString = "}\n",
                   nodeLabels =  nodesString3)
  }
  if (returnString){
    return(paste0(result))
  } else {
    return(result)
  }
}

#' internal function used by the nodes function
#' 
#' @param nodeParameter element from a (xmlToList type) workflow 
#' 
#' @return data.frame with parameters (or attributes)
#' 
#' @noRd
nodeParameters <- function(nodeParameter){
  if (".attrs"  %in% names(nodeParameter)){
    return(nodeParameter[[".attrs"]])
  } else {
    return(nodeParameter)
  }
}

#' internal function used by the nodes function
#' 
#' @param nodeParameter a (xmlToList type) workflow
#' 
#' @return character vector with the names of all nodes in the specified
#'  workflow
#'  
#' @noRd
nodeNames <- function(Workflow){
  if ("WorkflowTree" %in% names(Workflow)){
    unname(unlist(lapply(Workflow[["WorkflowTree"]], function(x){x[[".attrs"]]["FriendlyName"]})))
  } else {
    unname(unlist(lapply(Workflow, function(x){x[[".attrs"]]["FriendlyName"]})))
  }
}

#' function that takes a (xmlToList type) workflow and returns a list of 
#' 
#' @param Workflow a (xmlToList type) workflow
#' @param showHidden if TRUE then rows with hidden = TRUE are included (default:
#'  false)
#' @param showAdvanced if TRUE then rows with advanced = TRUE are included
#'  (default: TRUE)
#' @param showConfiguration if TRUE then rows with configuration = TRUE are
#'  included (default: FALSE)
#' 
#' @return a list of named data.frame objects containing all the parameters/
#'  settings in the nodes of the workflow
#' 
#' @note an example of it's use: (workflowInfo(db))$nodeInfo$Consensus %>%
#'  nodes()
#' 
#' @export
nodes <- function(Workflow,
                  showHidden = FALSE, showAdvanced = TRUE,
                  showConfiguration = FALSE){
  result <- lapply(Workflow$WorkflowTree, function(x){
    resultx <- dplyr::bind_rows(lapply(x$ProcessingNodeParameters, nodeParameters))
    if (ncol(resultx) == 0){
      return(resultx)
    } else {
      return(resultx %>%
           dplyr::select(.data$FriendlyName, .data$Category, .data$IsAdvanced,
                         .data$IsHidden, .data$IsConfig, .data$DisplayValue) %>%
               dplyr::rename(name = .data$FriendlyName,
                             category = .data$Category,
                             advanced = .data$IsAdvanced,
                             hidden = .data$IsHidden,
                             configuration = .data$IsConfig,
                             value = .data$DisplayValue))
    }
  })
  names(result) <- nodeNames(Workflow)
  if (!showHidden){
    result <- lapply(result, function(x){
      if (ncol(x)>0){
        return(x %>%
                 dplyr::filter(.data$hidden != "True") %>%
                 dplyr::select(-.data$hidden))
      } else {
        return(x)
      }
    })
  }
  if (!showConfiguration){
    result <- lapply(result, function(x){
      if (ncol(x)>0){
        return(x %>%
                 dplyr::filter(.data$configuration != "True") %>%
                 dplyr::select(-.data$configuration))
      } else {
        return(x)
      }
    })
  }
  if (!showAdvanced){
    result <- lapply(result, function(x){
      if (ncol(x)>0){
        return(x %>%
                 dplyr::filter(.data$advanced != "True") %>%
                 dplyr::select(-.data$advanced))
      } else {
        return(x)
      }
    })
  }
  return(result)
}