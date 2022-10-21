test_that("calcData works",{
  testResult <- calcData(mtcars)
  expect_equal(as.integer(testResult$median[1:4]), 
               c(4, 4, 4, 3))
  expect_equal(format(sum(testResult$median), digits = 6, nsmall = 4),
         "135.6440")
  testResult <- calcData(mtcars, calcName = "Mean", calcFunc = mean)
  expect_equal(format(testResult$Mean[1:4], digits = 6, nsmall = 4),
               c("29.9073", "29.9814", "23.5982", "38.7395"))
  expect_equal(format(sum(testResult$Mean), digits = 6, nsmall = 4),
               "1267.4729")
})

test_that("getProteinInfo(Raw) works", {
  testingDB <- openTest()
  testResult <- getProteinInfoRaw(db = testingDB)
  expect_equal(dim(testResult), c(5, 6))
  expect_equal(testResult$Accession |> sort(),
               knockOutProteins()$Accession |> sort())
  testResult <- getProteinInfo(db = testingDB)
  expect_equal(dim(testResult), c(5, 21))
  expect_equal(testResult$Accession |> sort(),
               knockOutProteins()$Accession |> sort())
  expect_equal(format(testResult$AbundancesNormalized_1,
                      digits = 6, nsmall = 4),
               c(" 1265.9276", " 5795.8521", " 7790.0660",
                 "11087.3421", "10720.0537"))
  closeTest(testingDB)
})

test_that("getPeptideInfo(Raw) works", {
  testingDB <- openTest()
  testResult <- getPeptideInfoRaw(db = testingDB)
  expect_equal(length(testResult), 5)
  expect_equal(dim(testResult[[5]]), c(26, 5))
  expect_equal(names(testResult) |> sort(),
               knockOutProteins()$Accession |> sort())
  testResult <- getPeptideInfo(db = testingDB)
  expect_equal(length(testResult), 5)
  expect_equal(dim(testResult[[5]]), c(26, 14))
  expect_equal(format(sum(testResult[[5]]$AbundancesNormalized_1),
                      digits = 6, nsmall = 4),
               "7790.0660")
  closeTest(testingDB)
})

test_that("calcIFI & calcAllIFIs work", {
  testingDB <- openTest()
  testResult <- calcIFIs(db = testingDB, groups = tmt10Channels())
  expect_equal(dim(testResult), c(22, 2))
  expect_equal(testResult$Short[1], "His4")
  expect_equal(format(mean(testResult$IFI), digits = 2, nsmall = 4), "0.6496")
  testResult <- calcAllIFIs(db = testingDB, groups = tmt10Channels())
  expect_equal(dim(testResult), c(102, 2))
  expect_equal(unique(testResult$Short), c("His4", "Met6", "Ura2"))
  expect_equal(format(mean(testResult[testResult$Short == "His4", ]$IFI), digits = 2, nsmall = 4), "0.6496")
  expect_equal(format(mean(testResult$IFI), digits = 2, nsmall = 4), "0.6753")
  closeTest(testingDB)
})