
test_that("dbGetTable works",{
  testingDB <- openTest(testfile = 1)
  dftest <- dbGetTable(db = testingDB,
                       tablename = "TargetProteins", SQL = T)
  expect_equal(dftest, "SELECT * FROM TargetProteins")
  dftest <- dbGetTable(db = testingDB,
                       tablename = "TargetProteins", SQL = F)
  expect_equal(nrow(dftest), 20)
  expect_equal(ncol(dftest), 64)
  # column names that have 'Abundance' in them
  expect_equal(length(colnames(dftest)[grepl(colnames(dftest),
                                             pattern = "Abundance")]), 13)
  # Abundances column is blob type
  expect_equal(sum(grepl(class(dftest$Abundances), pattern = "blob")), 1)
  # presence of knockout protein P00815
  expect_equal(sum(grepl(dftest$Accession, pattern = "P00815")), 1)
  closeTest(testingDB)
})

test_that("proteinIDTypes works", {
  testingDB1 <- openTest(testfile = 1)
  testingDB2 <- openTest(testfile = 2)
  expect_equal(proteinIDTypes(testingDB1)$GroupName,
               "Sequest HT")
  expect_equal(proteinIDTypes(testingDB2)$GroupName,
               c("Sequest HT","MSPepSearch"))
  closeTest(testingDB1)
  closeTest(testingDB2)
})

test_that("MSfileInfo works", {
  testingDB <- openTest(testfile = 1)
  testdf <- MSfileInfo(db = testingDB)
  expect_equal(testdf$RetentionTimeRange, "0.00 - 100.00")
  expect_equal(testdf$InstrumentDescription, "Q Exactive HF-X Orbitrap")
  expect_equal(grepl(testdf$FileName, pattern = "_top20_APDoff.raw"), TRUE)
  closeTest(testingDB)
})

test_that("SearchInfo works", {
  testingDB <- openTest(testfile = 1)
  testdf <- SearchInfo(db = testingDB)
  expect_equal(dim(testdf), c(324, 7))
  closeTest(testingDB)
})

test_that("totalSearchTime works", {
  testingDB <- openTest(testfile = 1)
  teststr <- format(totalSearchTime(db = testingDB), digits = 10)
  expect_equal(teststr, "2033.737044")
  closeTest(testingDB)
})

test_that("getAcquisitionDateTime & getAcquisitionDate work", {
  testingDB <- openTest(testfile = 1)
  teststr <- as.character(getAcquistionDateTime(db = testingDB))
  expect_equal(teststr, "2017-11-23 13:54:36")
  teststr2 <- as.character(getAcquistionDate(db = testingDB))
  expect_equal(teststr2, "2017-11-23")
  closeTest(testingDB)
})

test_that("dbGetProteinTable works",{
  testingDB <- openTest(testfile = 1)
  testdf <- dbGetProteinTable(db = testingDB)
  expect_equal(dim(testdf), c(17, 64))
  testdf <- dbGetProteinTable(db = testingDB, masterProtein = F)
  expect_equal(dim(testdf), c(20, 64))
  closeTest(testingDB)
})

test_that("dbGetProteinTable works",{
  testingDB <- openTest(testfile = 1)
  testdf <- dbGetProteinTable(db = testingDB)
  testSQL <- dbGetProteins(testingDB,
                           testdf$UniqueSequenceID[1:2],
                           columns = "Accession", SQL = T)
  expect_equal(testSQL,
               "SELECT Accession FROM TargetProteins WHERE UniqueSequenceID IN ('-9183054829930716487','-8166449411917027120')")
  testdf <- dbGetProteins(testingDB,
                          testdf$UniqueSequenceID[1:2],
                          columns = "Accession",
                          SQL = F)
  expect_equal(testdf$Accession, c("P06634","P39935"))
  closeTest(testingDB)
})

test_that("dbGetProteinFiltered works", {
  testingDB <- openTest(testfile = 1)
  testSQL <- dbGetProteinFiltered(testingDB, columns = "Accession",
                                  filtering = " Accession = 'P06634'", SQL = T)
  expect_equal(testSQL,
               "SELECT Accession FROM TargetProteins WHERE    Accession = 'P06634'")
  testSQL <- dbGetProteinFiltered(tesingDB, columns = "Accession",
                                  filtering = "(Accession = 'P06634' OR Accession = 'P39935')",
                                  masterProtein = TRUE,
                                  SQL = T)
  expect_equal(testSQL,
               "SELECT Accession FROM TargetProteins WHERE  IsMasterProtein = 0 AND (Accession = 'P06634' OR Accession = 'P39935')")
  testdf <- dbGetProteinFiltered(testingDB, columns = "Accession",
                                 filtering = " Accession = 'P06634'", SQL = F)
  expect_equal(testdf$Accession, "P06634")
  testdf <- dbGetProteinFiltered(testingDB,
                                 filtering = "(Accession = 'P06634' OR Accession = 'P39935')",
                                 masterProtein = TRUE,
                                 SQL = F)
  expect_equal(dim(testdf), c( 2, 64))
  expect_equal(testdf$Accession, c("P06634","P39935"))
  closeTest(testingDB)
})

test_that("dbGetProteinUniqueSequenceUDs works",{
  testingDB <- openTest(testfile = 1)
  testSQL <- dbGetProteinUniqueSequenceIDs(testingDB,
                                           accession = c("P06634","P39935"),
                                           SQL = T)
  expect_equal(testSQL,
               "SELECT UniqueSequenceID FROM TargetProteins WHERE   Accession IN ('P06634','P39935')")
  testdf <- dbGetProteinUniqueSequenceIDs(testingDB,
                                          accession = c("P06634","P39935"),
                                          SQL = F)
  expect_equal(as.character(testdf$UniqueSequenceID[1]), "-9183054829930716487")
  expect_equal(as.character(testdf$UniqueSequenceID[2]), "-8166449411917027120")
  closeTest(testingDB)
})

test_that("dbGetProteinGroups works",{
  testingDB <- openTest(testfile = 1)
  testSQL <- dbGetProteinGroups(db = testingDB,
                                proteinGroupIDs = c("1655","592","1131",
                                                    "1181","1428"),
                                SQL = T)
  expect_equal(testSQL,
               "SELECT * FROM TargetProteinGroups WHERE ProteinGroupID IN ('1655','592','1131','1181','1428')")
  testdf <- dbGetProteinGroups(db = testingDB,
                               proteinGroupIDs = c("1655","592","1131",
                                                   "1181","1428"),
                               SQL = F)
  expect_equal(dim(testdf),c(5, 10))
  expect_equal(sum(testdf$PeptideGroupCount), 160)
  closeTest(testingDB)
})

test_that("dbGetProteinGroupIDs works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetProteinGroupIDs(db = testingDB,
                                     proteinUniqueIDs = c("9147706995934957525",
                                                          "-2768548653852576336",
                                                          "3069497284891092343",
                                                          "3649982194632465886",
                                                          "6611310638582806557" ),
                                     SQL = F)[,1]
  expect_equal(testResult, c(592, 1131, 1181, 1428, 1655))
  closeTest(testingDB)
})

test_that("dbGetProteinIDs works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetProteinIDs(db = testingDB,
                                proteinGroupIDs = c(592, 1131, 1181, 1428, 1655),
                                SQL = F)[,1]
  testResult <- as.character(testResult)
  expect_equal(testResult, c("-2768548653852576336",
                             "3069497284891092343",
                             "3649982194632465886",
                             "6611310638582806557",
                             "9147706995934957525"))
  closeTest(testingDB)
})

test_that("dbGetProteinAnnotationGroupIDs works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetProteinAnnotationGroupIDs(
    db = testingDB,
    uniqueSequenceIDs = c("9147706995934957525",
                         "-2768548653852576336",
                         "3069497284891092343",
                         "3649982194632465886",
                         "6611310638582806557" ),
    SQL = F)[,1]
  expect_equal(sum(testResult), 642169)
  closeTest(testingDB)
})

test_that("dbGetAnnotatedProteins works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetAnnotatedProteins(db = testingDB,
                                       proteinAnnotationGroupIDs = c(264,
                                                                     14,
                                                                     1016,
                                                                     14475),
                                SQL = F)[,1]
  testResult <- as.character(testResult)
  expect_equal(testResult, c("-9183054829930716487",
                             "-8166449411917027120",
                             "-5974306765773997050",
                             "-2051302152036108770",
                             "-1452154832651305000",
                             "-635082257874195660",
                             "3069497284891092343",
                             "4093177069699858524"))
  closeTest(testingDB)
})

test_that("dbGetAnnotationGroups works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetAnnotationGroups(db = testingDB,
                                       proteinAnnotationGroupIDs = c(264,
                                                                     14,
                                                                     1016,
                                                                     14475),
                                       SQL = F)
  expect_equal(dim(testResult),c(3, 9))
  closeTest(testingDB)
})

test_that("dbGetAnnotationGroupsFiltered works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetAnnotationGroupsFiltered(
    db = testingDB,
    groupAnnotationAccession = "Pf03856",
    SQL = F)
  expect_equal(testResult$ProteinAnnotationGroupID[1], 264)
  expect_equal(testResult$GroupAnnotationAccession[1], "Pf03856")
  testResult <- dbGetAnnotationGroupsFiltered(
    db = testingDB,
    description = "Oxi", like = TRUE,
    SQL = F)
  expect_equal(testResult$ProteinAnnotationGroupID[1:2], c(422, 423))
  expect_equal(testResult$GroupAnnotationAccession[2:3], c("Pf02913","GO:0016898"))
  expect_equal(testResult$GroupAnnotationDescription[10], "oxidoreduction coenzyme metabolic process")
  closeTest(testingDB)
})

test_that("dbGetPeptideIDs works",{
  testingDB <- openTest(testfile = 1)
  testResult <- (dbGetProteinTable(db = testingDB)[1,] %>%
                   dbGetPeptideIDs(db = testingDB))[,1]
  expect_equal(sum(testResult), 138202)
  closeTest(testingDB)
})

test_that("dbGetPeptideTable works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPeptideTable(db = testingDB,
                                  peptideIDs = c(117, 222))
  expect_equal(testResult$Sequence, c("ACVVYGGSPIGNQLR",
                                      "AEIAIFGVPEDPNFQSSGINFDNYDDIPVDASGK"))
  closeTest(testingDB)
})

test_that("dbGetPPsmIDs works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPsmIDs(db = testingDB, peptideGroupIDs = 222)[,1]
  expect_equal(testResult, c(145617, 145484))
  closeTest(testingDB)
})

test_that("dbGetPsmTable works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPsmTable(db = testingDB,
                                  psmIDs = c(145617, 145484))
  expect_equal(testResult$ModifiedSequence,
               rep("[K].aEIAIFGVPEDPNFQSSGINFDNYDDIPVDASGk.[D]", 2))
  closeTest(testingDB)
})

test_that("dbGetConsensusIDs works",{
  testingDB <- openTest(testfile = 2)
  testResult <- dbGetPeptideTable(db = testingDB,
                                  columns = "PeptideGroupID")
  testResult <- dbGetConsensusIDs(db = testingDB,
                                  peptideGroupIDs = testResult$PeptideGroupID[1:10])[,1]
  expect_equal(testResult, c(98998, 140321, 99103, 34100, 128968, 108655,
                             132644, 5748, 69090, 54585, 140583))
  closeTest(testingDB)
})

test_that("dbGetConsensusTable works",{
  testingDB <- openTest(testfile = 2)
  testResult <- dbGetConsensusTable(db = testingDB,
                                    consensusIDs = c(98998, 140321, 99103,
                                                     34100, 128968, 108655,
                                                     132644, 5748, 69090,
                                                     54585, 140583))
  expect_equal(dim(testResult), c(11, 20))
  expect_equal(sum(testResult$ChargeState), 23)
  closeTest(testingDB)
})

test_that("dbGetQuanSpectrumIDs works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPsmTable(db = testingDB,
                              columns = "PeptideID")
  testResult <- dbGetQuanSpectrumIDs(db = testingDB,
                                     peptideIDs = testResult$PeptideID[1:10])[,1]
  expect_equal(testResult, c(28735, 28899, 28965, 29026,
                             29332, 29420, 29513, 29521, 29605, 29668))
  closeTest(testingDB)
})

test_that("dbGetQuanSpectrumInfoTable works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetQuanSpectrumInfoTable(
    db = testingDB,
    spectrumIDs = c(28735, 28899, 28965, 29026,
                    29332, 29420, 29513, 29521, 29605, 29668))
  expect_equal(dim(testResult), c(10, 33))
  expect_equal(sum(testResult$Charge), 22)
  closeTest(testingDB)
})


test_that("dbGetModificationsSitesID works",{
  testingDB <- openTest(testfile = 2)
  testResult <- dbGetProteinTable(db = testingDB,
                                  columns = "UniqueSequenceID")
  testResult <- dbGetModificationsSitesIDs(db = testingDB,
                                           proteinUniqueIDs = testResult$UniqueSequenceID[1:10])[,1]
  expect_equal(testResult, c(288, 289, 290, 291, 3069, 3070, 3071, 3072, 3073,
                             3074, 3075, 3076, 3077, 7584, 7715, 7716, 7717,
                             7718))
  closeTest(testingDB)
})

test_that("dbGetModificationsTable works",{
  testingDB <- openTest(testfile = 2)
  testResult <- dbGetModificationsTable(
    db = testingDB,
    modificatonSitesIDs = c(288, 289, 290, 291, 3069, 3070, 3071, 3072, 3073,
                            3074, 3075, 3076, 3077, 7584, 7715, 7716, 7717,
                            7718))
  expect_equal(dim(testResult), c(18, 21))
  expect_equal(sum(testResult$PositionInPeptide), 169)
  closeTest(testingDB)
})

test_that("dbGetModificationPeptideIDs works",{
  testingDB <- openTest(testfile = 2)
  testResult <- dbGetModificationPeptideIDs(
    db = testingDB,
    modificationIDs = c(288, 289, 290, 291, 3069, 3070, 3071, 3072, 3073,
                        3074, 3075, 3076, 3077, 7584, 7715, 7716, 7717,
                        7718))[,1]
  expect_equal(sum(testResult), 796078)
  closeTest(testingDB)
})

test_that("dbGetMSnSpectrumInfo works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPsmTable(db = testingDB,
                              columns = "PeptideID")
  testResult <- dbGetMSnSpectrumInfo(db = testingDB,
                                     peptideID = testResult)
  expect_equal(testResult$ScanNumbers[1:10],
               c("3705","3943","4044","4142","4594","4712","4836","4846","4961","5061"))
  expect_equal(nrow(testResult), 563)
  closeTest(testingDB)
})

test_that("dbGetMSnSpectrumInfo works",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPsmTable(db = testingDB,
                              columns = "PeptideID")
  testResult <- dbGetMassSpectrumItems(db = testingDB,
                                       peptideID = testResult)
  expect_equal(nrow(testResult), 563)
  closeTest(testingDB)
})

test_that("Spectrum routines work",{
  testingDB <- openTest(testfile = 1)
  testResult <- dbGetPsmTable(db = testingDB,
                              columns = "PeptideID")
  testResult <- dbGetMassSpectrumItems(db = testingDB,
                                       peptideID = testResult)
  testSp1 <- transformSpectrumRaw(testResult$Spectrum[[1]])
  testSp2 <- transformSpectrumRaw(testResult$Spectrum[[100]])
  expect_equal(names(testSp1), c('Header','ScanEvent','PrecursorInfo','PeakCentroids','ProfilePoints'))
  expect_equal(names(testSp2), c('Header','ScanEvent','PrecursorInfo','PeakCentroids','ProfilePoints'))
  expect_equal(spectrum.header(testSp1)$DataType, "Centroid")
  expect_equal(spectrum.scanEvent(testSp2)$ResolutionAtMass200, "45000")
  expect_equal(nrow(spectrum.centroid(testSp1)),85)
  expect_equal(colnames(spectrum.centroid((testSp1))), c('mz','intensity','charge','resolution','signalToNoise'))
  expect_equal(nrow(spectrum.centroid((testSp2))), 172)
  expect_equal(nrow(spectrum.precursor.centroid(testSp1)), 21)
  expect_equal(colnames(spectrum.precursor.centroid(testSp1)), c('mz','intensity','charge','resolution','signalToNoise'))
  expect_equal(spectrum.precursor.additionalInfo(testSp2)$Charge, "2")
  closeTest(testingDB)
})
