

fName <- "210126-nonHela-Combined-all.pdResult"
fDir <- "/home/ben/Documents/Thermo/Hela2"

dbName <- paste(c(fDir,"/",fName), collapse = "")

db1 <- db_open(dbName)
dbListTables(db1)
dbListFields(db1,"TargetProteins")
lapply(tbl(db1, "TargetProteins"), class)
gl <- glimpse(tbl(db1,"TargetProteins"), width= 0)


tbl(db1,"TargetProteins") %>% lapply(type_sum) %>% unlist

gl <- dbGetQuery(db1,"SELECT count(*) FROM TargetProteins")

gl <- dbGetQuery(db1, "SELECT * FROM TargetProteins LIMIT 1")

gl2 <- map_chr(gl,~class(.x)[1])

data.frame(name = names(gl2), type = unname(gl2))

dbGetQuery(db1,"SELECT * FROM pragma_info('TargetProteins')")

