% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tablecompet.R
\name{tablecompet}
\alias{tablecompet}
\title{tablecompet}
\usage{
tablecompet(data)
}
\arguments{
\item{data}{Data frame consisting survival events such, 'Id', 'OS', 'PFS', 'Death', 'Prog'.}
}
\description{
'tablecompet' function from 'progcomprisk' package enables user to input a
data and calculate the proportions of types of status in column 'Status'.
}
\details{
'tablecompet(data) adds a new column 'Status' to the input
dataset and assigns status 0, 1, 2 to each cell. It assigns 0 when progression is
1, death is 0 or progression 0 death is 0, assigns 1 when progression is 1 and
death is 1, whereas assigns 2 when progression is 0 and death is 1.

Concern should be given to rename the columns of the input data before using 'progcomprisk' package, as follows.

1) Subject Id column should be named as 'Id'.

2) OS column must be named as 'OS' .

3) Death status/event column should be named as 'Death'.

4) PFS column should be named as 'PFS'.

5) Progression event column should be named as 'Prog'.
}
\examples{
tablecompet(user_data)

}
