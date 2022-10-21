test_that("intoTable works",{
  testVector <- c("1","strain","", "Factor")
  names(testVector) <- c("Version", "Name", "Value", "Kind")
  testResult <- intoTable(testVector)
  expect_equal(dim(testResult), c(1, 4))
  expect_equal(unique(unlist(unname(lapply(intoTable(testVector), class)))),
               "character")
  testVector2 <- testVector
  testVector2["Name"] <- "Test"
  testResult <- intoTable(list(testVector, testVector2))
  expect_equal(dim(testResult), c(2, 4))
  expect_equal(colnames(testResult), c("Version", "Name", "Value", "Kind"))
})

test_that("analysisDefinition works", {
  testingDB <- openTest(testfile = 1)
  testResult <- analysisDefinition(db = testingDB)
  expect_equal(class(testResult), "list")
  expect_equal(length(testResult), 5)
  expect_equal(names(testResult), c('StudyDefinition',
                                    'AnalysisInputFileSets',
                                    'StudyAnalysisExtensionSettings',
                                    'StudyAnalysisExtensions',
                                    '.attrs'))
  closeTest(testDB = testingDB)
})

test_that("studyDefinitionFactors works", {
  testingDB <- openTest(testfile = 1)
  testResult <- analysisDefinition(db = testingDB)
  testResult <- studyDefinitionFactors(testResult)
  expect_equal(dim(testResult), c(1, 5))
  expect_equal(dim(testResult$factors[[1]]), c(4, 3))
  expect_equal(testResult$factors[[1]]$value, c("met6", "his4",
                                                "ura2", "wildtype"))
  closeTest(testDB = testingDB)
})

test_that("studyDefinitionQuanMethods works", {
  testingDB <- openTest(testfile = 1)
  testResult <- analysisDefinition(db = testingDB)
  testResult <- studyDefinitionQuanMethods(testResult)
  expect_equal(length(testResult), 2)
  expect_equal(testResult[[1]]$Name, "TMT 10plex")
  expect_equal(testResult[[2]]$Name[1:3], c("126", "127N", "127C"))
  closeTest(testDB = testingDB)
})

test_that("studyDefinitionSamples works", {
  testingDB <- openTest(testfile = 1)
  testResult <- analysisDefinition(db = testingDB)
  testResult <- studyDefinitionSamples(testResult)
  expect_equal(dim(testResult), c(10, 6))
  expect_equal(testResult$Name[9:10], c("TKOTMT10plex_50min_120k_45k_86ms_top20_APDoff - [130C]",
                                        "TKOTMT10plex_50min_120k_45k_86ms_top20_APDoff - [131]"))
  closeTest(testDB = testingDB)
})

test_that("studyDefinitionExtensions works", {
  testingDB <- openTest(testfile = 1)
  testResult <- analysisDefinition(db = testingDB)
  testResult <- studyDefinitionExtensions(testResult)
  expect_equal(length(testResult), 2)
  expect_equal(names(testResult), c("Method", "CorrectionFactors"))
  expect_equal(testResult$Method$selected, "MassTags")
  expect_equal(format(sum(as.numeric(testResult$CorrectionFactors$ReporterIon)),
                      digits = 10), "1286.329531")
  closeTest(testDB = testingDB)
})

test_that("studyDefinitionExtensionSettings works", {
  testingDB <- openTest(testfile = 1)
  testResult <- analysisDefinition(db = testingDB)
  testResult <- studyDefinitionExtensionSettings(testResult)
  expect_equal(length(testResult), 4)
  expect_equal(names(testResult), c("StudyVariablesForGrouping",
                                    "StudyVariablesForSorting",
                                    "QuanRatios",
                                    "XML"))
  expect_equal(dim(testResult$StudyVariablesForGrouping), c(1, 4))
  expect_equal(dim(testResult$StudyVariablesForSorting), c(1, 5))
  expect_equal(names(testResult$QuanRatios), c("his4/wildtype",
                                               "met6/wildtype",
                                               "ura2/wildtype"))
  expect_equal(dim(testResult$QuanRatios[[1]]$NumeratorSamples), c(3, 7))
  expect_equal(dim(testResult$QuanRatios[[1]]$DenominatorSamples), c(1, 7))
  expect_equal(length(testResult$XML), 4)
  closeTest(testDB = testingDB)
})