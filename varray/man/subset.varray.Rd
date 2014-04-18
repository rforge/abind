\name{subset.varray}
\alias{print.varray}
\alias{[.varray}
\alias{dim.varray}
\alias{dimnames.varray}
\alias{mode.varray}
\alias{storage.mode.varray}
\alias{vstorage.mode}
\alias{vstorage.mode.varray}
\alias{[<-.varray}
\alias{dim<-.varray}
\alias{dimnames<-.varray}
\alias{length<-.varray}
\alias{mode<-.varray}
\alias{storage.mode<-.varray}
\alias{vstorage.mode<-}
\alias{vstorage.mode<-.varray}
\title{
  Functions for extracting data and getting info about varrays
}
\description{
  Functions for extracting data and getting info about varrays.
  Replacement functions are also defined, but as \code{varray} objects
  are read-only, all the replacement functions stop with an error.
}
\usage{
\method{print}{varray}(x, ...)
\method{[}{varray}(x, ..., drop=TRUE)
\method{dim}{varray}(x)
\method{dimnames}{varray}(x)
\method{mode}{varray}(x)
\method{storage.mode}{varray}(x)
\method{vstorage.mode}{varray}(x, deep)
\method{vstorage.mode}{varray}(x) <- value
}
\arguments{
\item{x}{A \code{varray} object}
\item{\dots}{indexing arguments, treated as for ordinary array indexing}
\item{drop}{should dimensions with extent 1 in the result be dropped?}
\item{deep}{should individual components be inspected to determine the result?}
\item{value}{new values to use}
}
\details{
The functions \code{mode()} and \code{storage.mode} work by extracting a
sample from each sub-array, and concatenating the samples together to
work out the mode of the combined objects.  E.g., if one sub-array was
logical, and another numeric, the mode of the combined object would be
numeric.

The extraction functions \code{vstorage.mode} is similar to
\code{storage.mode}, but the replacement function \code{vstorage.mode<-}
will change the storage mode of all component objects.

}
\value{
The \code{[} method for \code{varray} returns an ordinary array.
}
\seealso{
\code{\link{varray}}
}
\examples{
a <- array(1:6, dim=c(2,3), dimnames=list(letters[1:2],letters[23:25]))
b <- array(7:15, dim=c(3,3), dimnames=list(letters[3:5],letters[24:26]))
x <- varray(a,b)
x
x[c('a'),c('x','z'),drop=FALSE]
x[c('d','b','c'),c('y','z'),drop=FALSE]
}
\keyword{misc}
