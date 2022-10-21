
test_that("convertRawNumeric works", {
  testVector <- as.raw(
    as.hexmode(c('00', '00', 'c0', 'cc', 'cc', '01', 'b0', '40', '01',
                 '00', '00', '40', '98', '99', '12', 'b0', '40', '01',
                 '00', '00', '80', '9a', '19', '8b', 'b0', '40', '01',
                 '00', '00', '40', '68', '66', 'e7', 'ad', '40', '01',
                 '00', '00', 'c0', '9d', '99', '17', 'ae', '40', '01',
                 '00', '00', '80', '33', '33', '3a', 'af', '40', '01',
                 '00', '00', '00', '00', '00', '00', '00', '00', '00',
                 '00', '00', '00', '00', '00', '00', '00', '00', '01',
                 'ff', 'ff', 'bf', '99', '99', '8b', 'af', '40', '01',
                 '00', '00', '00', '32', 'b3', '09', 'b0', '40', '01')))
  resultVector <- convertRawNumeric(testVector)
  expect_equal(length(resultVector), 10)
  expect_equal(format(resultVector[[1]], digits = 1, nsmall = 4), "4097.8000")
  expect_equal(resultVector[[7]], NA)
  expect_equal(format(resultVector[[8]], digits = 1, nsmall = 4), "0.0000")
  expect_equal(format(resultVector[[10]], digits = 1, nsmall = 4), "4105.7000")
})

test_that("convertRawInteger works", {
  testVector <- as.raw(
    as.hexmode(c('04', '00', '00', '00', '01',
                 '04', '00', '00', '00', '01',
                 '04', '00', '00', '00', '01',
                 '04', '00', '00', '00', '01',
                 '00', '00', '00', '00', '00',
                 '04', '00', '00', '00', '01',
                 '04', '00', '00', '00', '01',
                 '00', '00', '00', '00', '01',
                 '04', '00', '00', '00', '01', 
                 '04', '00', '00', '00', '01')))
  resultVector <- convertRawInteger(testVector)
  expect_equal(length(resultVector), 10)
  expect_equal(resultVector[[1]], 4)
  expect_equal(resultVector[[5]], NA)
  expect_equal(resultVector[[8]], 0)
  expect_equal(resultVector[[10]], 4)
})

test_that("convertRawSpecial & ConvertRawSpecialToString work", {
  testVector <- as.raw(
    as.hexmode(c('00', '00', '00', '00', '00', '00',
                 '00', '00', '00', '00', '01', '01',
                 '00', '00', '00', '00', '00', '00',
                 '00', '00', '00', '00', '00', '00',
                 '00', '00', '01', '01', '00', '00',
                 '00', '00', '00', '00', '00', '00')))
  resultVector <- convertRawSpecial(testVector, specialSize = 18)
  expect_equal(length(resultVector), 18)
  expect_equal(resultVector[[14]], TRUE)
  expect_equal(resultVector[[15]], NA)
  resultVector <- convertRawSpecialToString(resultVector)
  expect_equal(resultVector, "000001000000010000")
})

test_that("getBlobs works",{
  testingDB <- openTest(testfile = 1)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- getBlobs(testTable)
  expect_equal(dim(testResult), c(26, 2))
  expect_equal(testResult$name[3:5], c("NumberOfPeptidesbySearchEngine",
                                       "ScoreSequestHT",
                                       "ProteinFDRConfidence"))
  expect_equal(testResult$name[12:15], c("Abundances",
                                         "AbundancesCounts",
                                         "AbundancesNormalized",
                                         "AbundancesScaled"))
  closeTest(testDB = testingDB)
  testingDB <- openTest(testfile = 2)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- getBlobs(testTable)
  expect_equal(dim(testResult), c(15, 2))
  expect_equal(testResult$name[3:5], c("NumberOfPeptidesbySearchEngine",
                                       "ScoreSequestHT",
                                       "ScoreMSPepSearch"))
  expect_equal(testResult$name[12:15], c("Abundances",
                                         "AbundancesCounts",
                                         "FoundinSamples",
                                         "FoundinSampleGroups"))
  closeTest(testDB = testingDB)
  
})

test_that("determineBlobLengths works", {
  testingDB <- openTest(testfile = 1)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- determineBlobLengths(theTable = testTable)
  expect_equal(testResult$name[12], "Abundances")
  expect_equal(testResult$length[12], 90)
  expect_equal(sum(testResult$length), 795)
  closeTest(testDB = testingDB)
  testingDB <- openTest(testfile = 2)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- determineBlobLengths(theTable = testTable)
  expect_equal(testResult$name[12], "Abundances")
  expect_equal(testResult$length[12], 9)
  expect_equal(sum(testResult$length), 175)
  closeTest(testDB = testingDB)
})

test_that("determineBlobTypes works", {
  testingDB <- openTest(testfile = 1)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- determineBlobTypes(theTable = testTable)
  expect_equal(dim(testResult), c(26, 5))
  expect_equal(sum(testResult$minimumSize, na.rm = TRUE), 141)
  expect_equal(testResult$what[10], "special")
  expect_equal(testResult$what[11], as.character(NA))
  expect_equal(testResult$what[12], "numeric")
  expect_equal(testResult$what[13], "integer")
  closeTest(testDB = testingDB)
})

test_that("dfTransformRaws works", {
  testingDB <- openTest(testfile = 1)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- dfTransformRaws(testTable)
  expect_equal(dim(testResult), c(17, 133))
  expect_equal(colnames(testResult)[50:59], c("Abundances_1",
                                              "Abundances_2",
                                              "Abundances_3",
                                              "Abundances_4",
                                              "Abundances_5",
                                              "Abundances_6",
                                              "Abundances_7",
                                              "Abundances_8",
                                              "Abundances_9",
                                              "Abundances_10"))
  expect_equal(format(sum(testResult$Abundances_2), digits = 5, nsmall = 4),
               "84496.0000")
  closeTest(testDB = testingDB)
  testingDB <- openTest(testfile = 2)
  testTable <- dbGetProteinTable(testingDB)
  testResult <- dfTransformRaws(testTable)
  expect_equal(dim(testResult), c(4, 59))
  expect_equal(colnames(testResult)[56:57], c("Abundances_1",
                                              "AbundancesCounts_1"))
  expect_equal(format(sum(testResult$Abundances_1), digits = 5, nsmall = 4),
               "153310764.1562")
  closeTest(testDB = testingDB)
})