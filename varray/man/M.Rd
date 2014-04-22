\name{M}
\alias{M}
\title{
  Shorthand for as.matrix or as.array
}
\description{
  Shorthand for as.matrix or as.array.
}
\usage{
M(x)
}
\arguments{
  \item{x}{
}
}
\details{
  If \code{x} is an atomic object and is not a virtual array, returns
  \code{x}, otherwise returns \code{as.matrix(x)} if \code{x} has two
  dimensions or \code{as.array(x)} otherwise.
}
\value{
  An atomic array version of \code{x}
}
\author{
  Tony Plate \email{tplate@acm.org}
}
\seealso{
  \code{link{as.array}}
  \code{link{as.matrix}}
}
\examples{
x <- varray(cbind(A=c(a=1)), cbind(B=c(b=2)))
x
M(x)
as.array(x, fill=0)
as.matrix(x)
M(x) + 10
}
\keyword{ misc }
