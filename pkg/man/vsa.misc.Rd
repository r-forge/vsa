\name{vsa.misc}
\alias{vsa.misc}
\alias{simval}
\alias{scalar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions and objects for VSA similarity values. }
\description{

Functions to convert between \code{simval}s and ordinary scalar values.

}
\usage{
simval(x)
scalar(e1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric value (i.e., a vector of length 1 ) }
  \item{e1}{ a \code{simval} }
}
\details{

A \code{simval} ("Similarity value") is a numeric vector of length 1
with the S3 class \code{simval}.  Various functions, e.g., \code{dot()}
and \code{cosine()} return \code{simval} objects.  A object with a class
is used instead of an ordinary scalar so as to catch some possible
operator precedence errors such as writing \code{x * y \%.\% z} with the
(false) belief that it will be parsed as \code{(x * y) \%.\% z} when in
fact it will be parsed as \code{x * (y \%.\% z)}.  Having \code{\%.\%}
return a \code{simval} object allows the expression \code{x * y \%.\% z}
to result in an error, because it is illegal to multiple a vector by a
\code{simval}.  If \code{x * (y \%.\% z)} really was intended, then it
must be written as \code{x * scalar(y \%.\% z)}.

}

\value{
  \code{simval()} returns a vector with the S3 class \code{simval}.

  \code{scalar()} strips off the S3 class \code{simval} (i.e., it returns a numeric vector with no S3 class).

}
\author{ Tony Plate \email{tplate@acm.org}}

\seealso{ \code{\link{vsa}} }
\examples{
x <- newVec()
y <- newVec()
z <- newVec()
(x * y) \%.\% z
x * scalar(y \%.\% z)
}
\keyword{ misc }% at least one, from doc/KEYWORDS
