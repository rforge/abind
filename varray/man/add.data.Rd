\name{add.data.Matrix}
\alias{add.data.Matrix}
\alias{add.data.matrix}
\alias{update.Matrix}
\alias{update.matrix}
\title{
  Update the contents of a Matrix or matrix object
}
\description{
  Update the contents of a Matrix or matrix object, adding new dimension
  indices if necessary.
}
\usage{
add.data.Matrix(object, data, need.dimnames = list(NULL, NULL), keep.ordered = TRUE, ...)
add.data.matrix(object, data, need.dimnames = list(NULL, NULL), keep.ordered = TRUE, ...)
update.Matrix(object, data, need.dimnames = list(NULL, NULL), keep.ordered = TRUE, ...)
update.matrix(object, data, need.dimnames = list(NULL, NULL), keep.ordered = TRUE, ...)
}
\arguments{
  \item{object}{
    Character data specifying the name of the object to update.  Must
    have dimnames.
}
\item{data}{
  New data to incorporate in the object.  Should have the same number of
  dimensions as the object being updated.  Must
  have dimnames.
}
  \item{need.dimnames}{
Dimension names that should be included in the updated object.
}
\item{keep.ordered}{
  Logical.  Specifies which dimensions should be kept ordered.  Can be a
  single element or a vector with length equal to the number of
  dimensions of \code{object}.
}
\item{\dots}{
  Not used, but needed because \code{add.data()} could be a generic.
}
}
\details{
  Contents of \code{data} are inserted into \code{object}, expanding the
  dimensions of \code{object} if necessary.

  \code{update.matrix} is an obselete name for \code{add.data.matrix}.

  \code{update.Matrix} is an obselete name for \code{add.data.Matrix}.

  The \code{x} parameter is supplied by name rather than by value
  because it is modified.

}
\value{
  Invoked for its side effect, which is altering the value of the
  variable named by \code{object}.
}
\author{
Tony Plate \email{tplate@acm.org}
}
\note{
  Not really closely related to \code{varray} objects, but supplied here
  as a useful analogue to \code{\link{add.tsdata.varray}}.
}
\seealso{
  \code{\link{add.tsdata.varray}}
}
\examples{
x <- cbind(A=c(a=1))
add.data.matrix('x', cbind(B=c(b=2)))
x
}
\keyword{ misc }
