# ---- database related functions ----

# Note: proteome discoverer .pdResult files are in the SQLite format
# functions in this section deal with the database. Some of the
# functions work with database (query) commands, which is the preferred
# way. A lot of these functions are merely wrappers to make programming
# a little easier
#
# Note: for future developments when buidling a shiny app it's recommended
# to use the 'pool' & 'dbplyr' packages to streamline things

#' Wrapper around pool::dbPool(): opens a databse
#'
#' @param  fileName  a character vector specifying the name and location
#'                   of the database
#' @param  drv defines database connection type, default = RSQLite::SQLite()
#' @param  ... to pass on additional parameters to pool::dbPool, exmples are
#'             host = "shiny-demo.csa7qlmguqrf.us-east-1.rds.amazonaws.com"
#'             username = "guest"
#'             password = "guest"
#' @return database access 'handle'
#' @note if no file with the name 'fileName' exists, then it will be created
#' (but obviously it will be empty, so most further commands will fail)
#' @note if fileName == ":memory:" the database will be an in-memory database
#' @export
db_open <- function(fileName, drv = RSQLite::SQLite(), ...){
  return(pool::dbPool(drv = drv,
                      dbname = fileName,
                      ...))
}

#' Wrapper around pool::pooClose():  closes an open database
#' (normally opened earlier via eg db_open())
#'
#'
#' @param db   database access 'handle' to be closed
#'
#' @export
db_close <- function(db){
  pool::poolClose(db)
}

#' find the number of rows in a table of an (open) database
#'
#' @param db    database access 'handle'
#' @param tabl  name of the table
#' @return number of rows in the table
#'
#' @export
db_nrow <- function(db,tabl){
  return(pool::dbGetQuery(db,paste(c("SELECT COUNT(*) FROM ",tabl),
                                   collapse = ""))[1,1])
}

#' get the classes of the different columns in a database table
#'
#' @param db    database access 'handle'
#' @param tabl  name of the table
#' @return a data.frame with columns name & type (= class)
#'
#' @export
db_columnInfo <- function(db, tabl){
  # get first row of the table to determine names & classes
  temp <- pool::dbGetQuery(db1, "SELECT * FROM TargetProteins LIMIT 1")
  # determine classes, note that in case of class lists
  # (as with type OOP descendants) this statement only take the first name
  # of the list
  temp <- map_chr(temp,~class(.x)[1])
  return(data.frame(name = names(temp), type = unname(temp)))
}

#' gathers info on all tables in the database
#'
#' @param db    database access 'handle'
#' @return a data.frame with the following columns: name, nrows, info
#' every entry in the column info contains a data.frame with the names
#' of the columns in the table and their type (see also db_columnInfo)
#' @export
db_tbl_def <- function(db){
  datazz <- pool::dbListTables(db)
  dftbl <- data.frame(name = datazz,
                      nrows = unlist(lapply(datazz,
                                            function(x){db_nrow(db,x)})),
                      info = NA)
  for (counter in 1: length(datazz)){
    dftbl$info[counter] <- list(db_columnInfo(db,datazz[counter]))
  }
  return(dftbl)
}
