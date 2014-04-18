\name{is.virtual.array}
\alias{is.virtual.array}
\alias{is.virtual.array.varray}
\alias{is.virtual.array.default}
\title{
  Indicate whether an array-like object is native or virtual
}
\description{
  Indicate whether an array-like object is native or virtual.  A native
  array object is a primitive vector object with a \code{dim} attribute.
  A virtual array object is some object that behaves like an array.
  This is a generic function with methods defined for class \code{varray}.
}
\usage{
is.virtual.array(x)
\method{is.virtual.array}{varray}(x)
\method{is.virtual.array}{default}(x)
}
\arguments{
  \item{x}{
    An array-like object
}
}
\value{
  A logical value, \code{TRUE} or \code{FALSE}
  }
\author{
    Tony Plate \email{tplate@acm.org}
}
\keyword{ misc }
