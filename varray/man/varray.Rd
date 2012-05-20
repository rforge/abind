\name{varray}
\alias{varray}
\alias{as.array.varray}
\title{
  Construct a virtual array
}
\description{
  Construct a virtual array by binding together supplied arrays or
  matrixes.  Stores only a by-name reference to the sub-arrays; does
  not store their actual data.
}
\usage{
varray(\dots, along = 1, dimorder = NULL, env.name = TRUE, envir = NULL)
\method{as.array}{varray}(x, \dots)
}
\arguments{
  \item{\dots}{
    Either a single  argument being a character vector naming the component sub-arrays, or
    multiple arguments, each of which is a single character item or
    unquoted name of an object.
    Note that extra arguments are ignored for \code{as.array} (they are
    only present because the generic has them).
}
\item{along}{
  The dimension along which to bind (from the user's point of view).
}
\item{dimorder}{
  The order in which dimensions are stored in the sub-arrays.  The
  reverse of this permutation is applied to extract data from the sub-arrays.
}
\item{env.name}{
  \code{TRUE} or \code{FALSE}, specifiying whether to record the name of
  the environment in which the object is found.  If \code{FALSE},
  objects are searched for in the global environment.
}
\item{envir}{
  Where to find the sub-arrays (can be an environment name, or an
  environment, as long as the environment can be recovered by \code{as.environment(environmentName(envir))}).
}
\item{x}{
  A \code{varray} object.
}
}
\details{
  Sub-arrays are stored by reference by name.  At the time of object
  extraction the sub-arrays must be retrieved.  This creates the
  possibility for retrieving an unintended object of the same name.
  To minimize this possibility, the environments which are searched to
  find the sub-arrays are limited as follows:
  \itemize{
    \item If \code{env.name} is \code{FALSE} the global environment and
    the rest of the search path
    \item If \code{env.name} is \code{TRUE} the supplied environment
    \code{envir}, or, if that was \code{NULL}, the environment in which
    the object was found at the time the \code{varray} was constructed.
  }
  Note that this makes it currently impossible for a \code{varray} to
  refer to objects that reside in an unnamed environment or one not on
  the search list of environments.
}
\value{
  \code{varray()} returns an S3 object of class \code{varray}.
  \code{as.array()} returns a standard R array (which is a matrix when
  there are just two dimensions).
}
\author{
  Tony Plate \email{tplate@acm.org}
}
\seealso{
  \code{\link{[.varray}}
  \code{\link{dim.varray}}
  \code{\link{dimnames.varray}}
  \code{\link{length.varray}}
  \code{\link{mode.varray}}
  \code{\link{storage.mode.varray}}
}
\examples{
a <- array(1:6, dim=c(2,3), dimnames=list(letters[1:2],letters[23:25]))
b <- array(7:15, dim=c(3,3), dimnames=list(letters[3:5],letters[24:26]))
x <- varray(a,b)
x1 <- varray('a','b') # equivalent
x2 <- varray(c('a','b')) # equivalent
as.array(x)
x
x[c('a'),c('x','z'),drop=FALSE]
x[c('d','b','c'),c('y','z'),drop=FALSE]
}
\keyword{ misc }
