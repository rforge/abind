\name{rm.varray}
\alias{rm.varray}
\title{
  Remove a varray and its component objects
}
\description{
  Remove a varray and its component objects
}
\usage{
rm.varray(x, list = NULL)
}
\arguments{
  \item{x}{The name of a \code{varray} object to remove, either quoted
    or unquoted.
}
\item{list}{
  A character vector of names of varray objects to be removed.
}
}
\details{
  Removes all the component objects of a varray, then removes the varray
  object itself.
}
\value{
  Returns the vector of names of objects as an invisible object.
}
\author{
  Tony Plate \email{tplate@acm.org}
}
\keyword{ misc }
