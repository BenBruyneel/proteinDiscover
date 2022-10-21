test_that("df_replace works", {
  testdf <- data.frame(column1 = c("1\\+", "email@email.com", "Test"),
                       column2 = c(1,2,3),
                       column3 = c("Testing", "Testing +", "Test @ 200"))
  testdf <- df_replace(testdf)
  expect_equal(testdf$column1, c("1\\+", "email $@$ email.com", "Test"))
  expect_equal(testdf$column3, c("Testing", "Testing +",  "Test  $@$  200"))
})

test_that("workflowInfo works", {
  testingDB <- openTest(testfile = 1)
  testResult <- workflowInfo(db = testingDB)
  expect_equal(class(testResult), "list")
  expect_equal(length(testResult), 2)
  expect_equal(names(testResult), c("workflowInfo",
                                    "nodeInfo"))
  expect_equal(dim(testResult[[1]]), c(2, 10))
  expect_equal(names(testResult[[2]]), c("Consensus", "Processing"))
  closeTest(testDB = testingDB)
})

test_that("nodeTable works", {
  testingDB <- openTest(testfile = 1)
  testWF <- workflowInfo(db = testingDB)
  testConsensus <- nodeTable(testWF$nodeInfo$Consensus)
  expect_equal(dim(testConsensus), c(14, 5))
  expect_equal(testConsensus$parent[8], "5;0")
  expect_equal(testConsensus$name[1], "MSF Files")
  testProcessing <- nodeTable(testWF$nodeInfo$Processing)
  expect_equal(dim(testProcessing), c(5, 5))
  expect_equal(testProcessing$parent[4], "2")
  expect_equal(testProcessing$name[1], "Spectrum Files RC")
  closeTest(testDB = testingDB)
})

test_that("createDiagrammeRString works", {
  testingDB <- openTest(testfile = 1)
  testWF <- workflowInfo(db = testingDB)
  testConsensus <- nodeTable(testWF$nodeInfo$Consensus)
  testVector <- createDiagrammeRString(testConsensus)
  expect_equal(class(testVector), "character")
  expect_equal(length(testVector), 11)
  expect_equal(substr(testVector[[1]],1,17), "digraph thegraph ")
  expect_equal(substr(testVector[[3]],1,25), "node [label = '@@1'] n0 \n")
  closeTest(testDB = testingDB)
})

test_that("nodes works", {
  testingDB <- openTest(testfile = 1)
  testWF <- workflowInfo(db = testingDB)
  testConsensus <- allNodesTable(nodes(testWF$nodeInfo$Consensus))
  expect_equal(dim(testConsensus), c(118, 5))
  expect_equal(colnames(testConsensus), c("node",
                                          "category",
                                          "name",
                                          "advanced",
                                          "value"))
  expect_equal(testConsensus$value[118],
               "Only unique peptides based on protein groups")
  testProcessing <- allNodesTable(nodes(testWF$nodeInfo$Processing))
  expect_equal(dim(testProcessing), c(107, 5))
  expect_equal(colnames(testProcessing), c("node",
                                           "category",
                                           "name",
                                           "advanced",
                                           "value"))
  expect_equal(testProcessing$value[107], "1000")
  closeTest(testDB = testingDB)
})