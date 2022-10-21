
openTest <- function(testfile = 1){
  if (testfile == 1){
    return(dbOpen(test_path("fixtures", "test1.pdResult")))
  } else {
    return(dbOpen(test_path("fixtures", "test2.pdResult")))
  }
}

closeTest <- function(testDB){
  dbClose(testDB)
}
