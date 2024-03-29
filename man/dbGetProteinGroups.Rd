% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteinDatabase.R
\name{dbGetProteinGroups}
\alias{dbGetProteinGroups}
\title{Gets the ProteinGroup information from the TargetProteinGroups table}
\usage{
dbGetProteinGroups(db, proteinGroupIDs, columns = NA, SQL = FALSE)
}
\arguments{
\item{db}{database access 'handle'}

\item{proteinGroupIDs}{specifies which protein groups to get, these values
can come from eg the protein table}

\item{columns}{character vector, specifies which columns to retrieve}

\item{SQL}{allows the function to return the SQL query statement in stead of
a data.frame#'}
}
\value{
a data.frame or a character vector (SQL)
}
\description{
Gets the ProteinGroup information from the TargetProteinGroups table
}
