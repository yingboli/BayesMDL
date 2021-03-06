% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{loc2dirac}
\alias{loc2dirac}
\title{Convert a Location Vector to a 0-1 Indicator Vector}
\usage{
loc2dirac(loc, N)
}
\arguments{
\item{loc}{An integer vector, length does not exceed \code{N}.}

\item{N}{The length of the resulting indicator vector.}
}
\value{
A 0-1 indicator vector of length \code{N}.
}
\description{
Convert a location vector to a 0-1 indicator vector, such that the
  corresponding entries in these locations are 1.
}
\keyword{internal}
