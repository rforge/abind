\name{varray}
\alias{varray}
\alias{rebuild.varray}
\alias{as.array.varray}
\alias{as.matrix.varray}
\title{
  Construct a virtual array
}
\description{
  Construct a virtual array by binding together supplied arrays or
  matrixes.  Stores only a by-name reference to the sub-arrays; does
  not store their actual data.
}
\usage{
varray(\dots, along = 1, dimorder = NULL, env.name = FALSE, envir = NULL,
       naidxok = NA, dimnames = NULL, comp.name = NULL,
       keep.ordered = TRUE, umode = NULL, fill = NULL)
\method{as.array}{varray}(x, \dots)
\method{as.matrix}{varray}(x, \dots)
rebuild.varray(x)
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
\item{naidxok}{
  Logical value indicating whether the component objects can handle
  \code{NA}'s in indices.  Notably, \code{ff} objects cannot.  The
  default value is \code{FALSE} for objects that inherit from class
  \code{ff} and \code{TRUE} for other objects.
}
\item{dimnames}{
  Optional dimnames to use instead of the union of dimnames of the
  components.
  Must be a list of character vectors.  Null components are ignored.
  E.g., supplying \code{dimnames=list(NULL, c('a','b','c'))} would tell
  varray to use the usual rownames (i.e., the union of the rownames of
  the components) but to use \code{c('a','b','c')} as column names.  If
  some components have other column names, those other columns will be
  inaccessible through the varray object.
}
\item{x}{
  A \code{varray} object.
}
\item{comp.name}{
  A character string that records how to construct the name of a
  component.  Used in \code{update.ts.varray()}.
}
\item{keep.ordered}{
  Should dimnames be kept ordered when new elements are added?
}
\item{umode}{
  The storage mode of the component objects.
}
\item{fill}{
  The value which is returned instead of explicit or implicit NAs.
}
}
\details{
  Component arrays are stored by reference (by name.)  At the time of subset
  extraction the component arrays will be retrieved.  This creates the
  possibility for retrieving an unintended object of the same name.
  To minimize this possibility, the environments searched to
  find the component arrays are limited to the following:
  \itemize{
    \item If \code{env.name} is \code{FALSE}; the global environment and
    the rest of the search path
    \item If \code{env.name} is \code{TRUE}; the supplied environment
    \code{envir}, or, if that was \code{NULL}, the environment in which
    the component was found at the time the \code{varray} was constructed.
  }
  Note that this makes it currently impossible for a \code{varray} to
  refer to objects that reside in an unnamed environment or one not on
  the search list of environments.

  All the data associated with a single element of the binding dimension
  (the \code{along} dimension) must be stored in a single one of the
  component objects.

  \code{rebuild.varray()} rebuilds a varray object and can be used when
  the underlying objects have changed.

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
b <- b[-2,]
x <- rebuild.varray(x)
}
\keyword{ misc }
