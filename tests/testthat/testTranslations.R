test_that("isMasterProtein & pQuanInfo work",{
  expect_equal(isMasterProtein(c(1:5, 1:2)), c("Master Protein Candidate",
                                               "<Undefined>",
                                               "Master Protein Rejected",
                                               "Master Protein Rejected", 
                                               NA,
                                               "Master Protein Candidate",
                                               "<Undefined>"))
  expect_equal(pQuanInfo(c(1,3,4, 4:7)), c( "No Quan Values",
                                            "Redundant",
                                            "No Proteins",
                                            "No Proteins",
                                            "Filtered Out",
                                            "No Quan Labels",
                                            "Inconsistently Labelled"))
})
