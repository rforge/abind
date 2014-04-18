\name{add.tsdata.varray}
\alias{add.tsdata.varray}
\title{
  Update the contents of a varray object that stores time series data
}
\description{
  Update the contents of a varray object that stores time series data
}
\usage{
add.tsdata.varray(object,
                 data,
                 comp.name = va$comp.name,
                 dateblock = "\%Y",
%                 dates.by = "bizdays",
%                 holidays = "NYSEC",
%                 vmode = "single",
                 along = va$along,
                 dimorder = va$dimorder,
                 env.name = va$env.name,
                 envir = NULL,
                 naidxok = va$naidxok,
                 keep.ordered = va$keep.ordered,
                 umode = NULL,
                 store.env.name = FALSE,
                 fill = NA,
                 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    Character data specifying the name of the varray object to update.  Must
    have dimnames.
}
  \item{data}{
  New data to incorporate in the object.  Should have the same number of
  dimensions as the object being updated.  Must
  have dimnames.
}
\item{comp.name}{
  Base of name to use for the component objects in the varray.
}
  \item{dateblock}{
%%     ~~Describe \code{dateblock} here~~
}
%%  \item{dates.by}{}
%%  \item{holidays}{}
%%  \item{vmode}{}
  \item{along}{
%%     ~~Describe \code{along} here~~
}
  \item{dimorder}{
%%     ~~Describe \code{dimorder} here~~
}
  \item{env.name}{
%%     ~~Describe \code{env.name} here~~
}
  \item{envir}{
%%     ~~Describe \code{envir} here~~
}
\item{naidxok}{
  Set this attribute on the varray.  Specifies whether the component objects can handle \code{NA} indices.
}
  \item{keep.ordered}{
  Logical.  Specifies which dimensions should be kept ordered.  Can be a
  single element or a vector with length equal to the number of
  dimensions of \code{object}.
}
  \item{umode}{
%%     ~~Describe \code{umode} here~~
}
  \item{store.env.name}{
%%     ~~Describe \code{store.env.name} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{\dots}{
  Not used, but needed because \code{update()} is a generic.
}
}
\details{
  Contents of \code{data} are inserted into \code{object}, expanding the
  dimensions of \code{object} and creating new component objects if necessary.
}
\value{
  Invoked for its side effect, which is altering the value of the
  variable named by \code{object}.
}
\author{
Tony Plate \email{tplate@acm.org}
}

\seealso{
  \code{\link{add.data.Matrix}}, \code{\link{add.data.matrix}}
}
\examples{
x <- array(1:45, dim=c(15, 3), dimnames=list(dateSeq('2001-04-19', len=15, by='months'), c('A','B','C')))
y <- array(101:121, dim=c(7, 3), dimnames=list(dateSeq(rownames(x)[nrow(x)], len=7, by='months'), c('B','C','D')))
z <- array(1:72, dim=c(24, 3), dimnames=list(dateSeq('2001-04-19', len=24, by='months'), c('A','B','C')))
as.array(varray(x, y, along=1))
add.tsdata.varray('v', z, comp.name='vvv.\%Y', along=1)
as.array(v)
add.tsdata.varray('v', y)
as.array(v)
}
\keyword{ misc }
