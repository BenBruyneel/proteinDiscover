purrr::walk(list.files(path = "R", pattern = "*.R",
                       full.names = TRUE), source)

library(pool)
library(dbAccess)
library(dplyr)

tkodb <- db_open("/home/ben/Documents/Thermo/Hela2/210105-data07.pdResult")

tkoProt <- dbGetProteinTable(db = tkodb)
# 
# blobLength <- function(blobList){
#   return(mean(unlist(lapply((lapply(blobList,length)),
#                             function(x){ifelse(x==0,NA,x)})),
#               na.rm = TRUE))
# }
# 
# db_getBlobs <- function(db, tableName = tableName()){
#   return(db_columnInfo(db = db, tableName = tableName) %>%
#            filter(type == "blob"))
# }
# 
# db_getBlobs(db = tkodb, tableName = tableName()) %>%
#   filter(!grepl(name, pattern = "Aspect"))
# 
# determineBlobLengths <- function(columnInfo, theTable){
#   columnInfo$length <- unlist(lapply(1:nrow(columnInfo),function(x){blobLength(theTable[,columnInfo$name[x]])}))
#   return(columnInfo)
# }

determineBlobLengths(blobDF = db_getBlobs(db = tkodb, tableName = tableName()),
                     theTable = tkoProt)
determineBlobLengths(blobDF = getBlobs(tkoProt),
                     theTable = tkoProt)

determineBlobLengths(blobDF = db_getBlobs(db = tkodb, tableName = tableName()), theTable = tkoProt)

frcblb <- determineBlobLengths(blobDF = db_getBlobs(db = tkodb, tableName = tableName()),
                               theTable = tkoProt)
frcblb

# determineBlobTypeRaw <- function(blobLength){
#   return(
#     unlist(
#       lapply(blobLength, function(x){
#         switch(toString(x),
#                "5"="integer",
#                "9"="numeric",
#                NA)
#       }
#       )
#     )
#   )
# }

determineBlobTypeRaw(9)
determineBlobTypeRaw(36)

# # minimumNumber is first checked only when preferMinimumNumber = TRUE
# determineBlobType <- function(blobLength, minimumNumber,
#                               numberOfGroups = minimumNumber, preferNumberOfGroups = TRUE,
#                               ratioNumberOfGroups = numberOfGroups - 1, useRatios = TRUE){
#   if (blobLength %in% c(5,9)){
#     return(data.frame(what = determineBlobTypeRaw(blobLength), minimumSize = 1))
#   }
#   if (preferNumberOfGroups){
#     if ((blobLength %% numberOfGroups) == 0){
#       return(data.frame(what = determineBlobTypeRaw(blobLength = blobLength %/% numberOfGroups),
#                         minimumSize = numberOfGroups))
#     } else {
#       if ((blobLength %% minimumNumber) != 0){
#         if (useRatios){
#           if ((blobLength %% ratioNumberOfGroups) == 0){
#             return(data.frame(what = determineBlobTypeRaw(blobLength = blobLength %/% ratioNumberOfGroups),
#                               minimumSize = ratioNumberOfGroups))
#           }
#         }
#         warning("Possibly invalid blob type")
#       }
#       return(data.frame(what = determineBlobTypeRaw(blobLength = blobLength %/% minimumNumber),
#                         minimumSize = minimumNumber))
#     }
#   } else {
#     if ((blobLength %% minimumNumber) == 0){
#       return(data.frame(what = determineBlobTypeRaw(blobLength = blobLength %/% minimumNumber),
#                         minimumSize = minimumNumber))
#     } else {
#       if ((blobLength %% numberOfGroups) != 0){  # unsure if this is needed/good
#         if (useRatios){
#           if ((blobLength %% ratioNumberOfGroups) == 0){
#             return(data.frame(what = determineBlobTypeRaw(blobLength = blobLength %/% ratioNumberOfGroups),
#                               minimumSize = ratioNumberOfGroups))
#           }
#         }
#         warning("Possibly invalid blob type!")
#       }
#       return(data.frame(what = determineBlobTypeRaw(blobLength = blobLength %/% numberOfGroups),
#                         minimumSize = numberOfGroups))
#     }
#   } 
# }


determineBlobType(5, numberOfGroups = 4, minimumNumber = 11)
determineBlobType(9, numberOfGroups = 4, minimumNumber = 11)
determineBlobType(20, numberOfGroups = 4, minimumNumber = 11)
determineBlobType(99, numberOfGroups = 4, minimumNumber = 11)
determineBlobType(55, minimumNumber = 11, numberOfGroups = 4, ratioNumberOfGroups = 3)
determineBlobType(21, numberOfGroups = 3, minimumNumber = 11)
determineBlobType(18, minimumNumber = 11, numberOfGroups = 4, ratioNumberOfGroups = 2)

determineBlobType(5, numberOfGroups = 4, minimumNumber = 11, preferNumberOfGroups = FALSE)
determineBlobType(9, numberOfGroups = 4, minimumNumber = 11, preferNumberOfGroups = FALSE)
determineBlobType(20, numberOfGroups = 4, minimumNumber = 11, preferNumberOfGroups = FALSE)
determineBlobType(99, numberOfGroups = 4, minimumNumber = 11, preferNumberOfGroups = FALSE)
determineBlobType(55, numberOfGroups = 4, minimumNumber = 11, preferNumberOfGroups = FALSE)
determineBlobType(21, numberOfGroups = 4, minimumNumber = 11, preferNumberOfGroups = FALSE)

blobEstimateTypes <- function(blobLengths, minimumNumber,
                              numberOfGroups = minimumNumber,
                              ratioNumberOfGroups = numberOfGroups - 1){
  return(
    bind_rows(
      lapply(blobLengths, function(x){determineBlobType(blobLength = x,
                                                        numberOfGroups = numberOfGroups,
                                                        minimumNumber = minimumNumber,
                                                        ratioNumberOfGroups = ratioNumberOfGroups)})
    )
  )
}






blobEstimateTypes(c(5,12,9,18,27), numberOfGroups = 4, minimumNumber = 5)

blobEstimateTypes(frcblb$length, numberOfGroups = 4, minimumNumber = 11)
blobEstimateTypes((frcblb %>% filter(!grepl(name, pattern = "Aspect")))$length, numberOfGroups = 4, minimumNumber = 11)

blobEstimateTypes(blobLengths = frcblb$length, numberOfGroups = 4, minimumNumber = 11)

bind_cols(frcblb,blobEstimateTypes(blobLengths = frcblb$length, numberOfGroups = 4, minimumNumber = 11))

frcblb2 <- frcblb
frcblb2[12:17,]$length <- c(99,55,99,36,20,5)
bind_cols(frcblb2,blobEstimateTypes(blobLengths = frcblb2$length, numberOfGroups = 4, minimumNumber = 11))

determineBlobLengths(blobDF = db_getBlobs(db = tkodb, tableName = "TargetProteins") %>%
                       filter(!grepl(name, pattern = "Aspect")), theTable = tkoProt)

determineBlobTypes <- function(blobDF, theTable, minimumNumber,
                               numberOfGroups = minimumNumber,
                               ratioNumberOfGroups = numberOfGroups - 1,
                               specials = TRUE){
  blobDF <- determineBlobLengths(blobDF = blobDF, theTable = theTable)
  blobDF <- bind_cols(blobDF, blobEstimateTypes(blobLengths = blobDF$length,
                                                minimumNumber = minimumNumber,
                                                numberOfGroups = numberOfGroups,
                                                ratioNumberOfGroups = ratioNumberOfGroups))
  return(blobDF)
}

# blobDF = db_getBlobs(db = tkodb, tableName = "TargetProteins"), 
determineBlobTypes(theTable = tkoProt,
                   minimumNumber = 1, numberOfGroups = 1,ratioNumberOfGroups = 1)



tkol <- determineBlobTypes(theTable = tkoProt,
                           minimumNumber = 1, numberOfGroups = 1,ratioNumberOfGroups = 1)

tkoProtTranslated <- tkoProt %>% dfTransformRaws(blobDF = tkol)
tkoProtTranslated2 <- tkoProt %>% dfTransformRaws()
tkoProtTranslated3 <- tkoProt %>% dfTransformRaws(specials = FALSE)

determineRawTypes(df = tkoProt[,colnames(tkoProt) %in% getBlobs(tkoProt)$name])
