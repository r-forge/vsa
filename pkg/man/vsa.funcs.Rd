\name{vsa.funcs}
\alias{vsa.funcs}
\alias{vsaprod}
\alias{appinv}
\alias{add}
\alias{addnorm}
\alias{dot}
\alias{\%.\%}
\alias{equiv}
\alias{\%==\%}
\alias{cosine}
\alias{\%cos\%}
\alias{norm}
\alias{mag}
\alias{vsapower}
\alias{vsascale}
\alias{!.vsa}
\alias{*.vsa}
\alias{/.vsa}
\alias{+.vsa}
\alias{-.vsa}
\alias{^.vsa}
\alias{\%.\%.vsa}
\title{Generic algebraic operations and functions for VSA vectors }
\description{
  Generic functions for working with VSA vectors.  Functions for all
standard algebraic operations on vectors: multiplication, addition,
superposition, scaling, powers (including inverses), approximate
inverse, similarity.
}
\usage{
vsaprod(e1, e2, method = c("fft", "outer"))
appinv(e1)
\special{!e1}
% \method{!}{vsa}(e1)
add(e1, ...)
norm(e1)
addnorm(...)
dot(e1, e2)
e1 \%.\% e2
equiv(e1, e2, tol)
e1 \%==\% e2
cosine(e1, e2, mag1=NULL, mag2=NULL)
e1 \%cos\% e2
mag(e1, actual=NULL)
vsapower(e1, e2)
vsascale(e1, e2)
\special{e1 * e2}
\special{e1 + e2}
\special{e1 / e2}
\special{e1 - e2}
\special{e1 ^ e2}
% \method{*}{vsa}(e1, e2)
% \method{+}{vsa}(e1, e2)
% \method{/}{vsa}(e1, e2)
% \method{-}{vsa}(e1, e2)
% \method{^}{vsa}(e1, e2)
}
\arguments{
  \item{e1}{ a \code{vsa} vector -- a subclass of \code{vsa} }
  \item{e2}{ a \code{vsa} vector for \code{vsaprod}, \code{dot}, \code{equiv}, \code{cosine}, and \code{add}, or a scalar (a numeric vector of length 1) for \code{vsapower} and \code{vsascale} }
  \item{method}{ computational method for convolution (\code{vsaprod}) }
  \item{tol}{ numerical tolerance for \code{equiv} }
  \item{mag1,mag2}{ precomputed magnitudes for \code{e1} and \code{e2} }
  \item{actual}{work with data values from this object (in this mode of
    invocation \code{e1} is a template holding the class)}
  \item{...}{ addtional \code{vsa} vectors for \code{add} or  \code{addnorm} }
}
\details{

All of these functions take as arguments one or more VSA vectors, and
possibily a scalar.  A VSA vector is a subclass of \code{vsa} (e.g.,
\code{realhrr}).  The result is a VSA vector, or a scalar with class
\code{simval}.  To catch some errors due to misunderstanding of R's
precedence, a "simval" scalar cannot be supplied as the argument of
one of these functions (to do, remove the \code{simval} class using
the function \code{scalar}).

  \itemize{
  \item \code{mag} computes the magnitude of a vector (scalar value)
  \item \code{norm} normalizes a vector so that its magnitude is 1 (vector value)
  \item \code{vsaprod} and its synonym \code{*} compute the vector product of two vectors (vector value)
  \item \code{appinv} and its synonym \code{!} computes an approximate inverse of a vector.  The expression \code{x / y} where \code{x} and \code{y} are vectors computes \code{x * appinv(y)}.  (vector value)
  \item \code{vsascale} and its synonym \code{*} multiple a vector by a scalar (vector value)
  \item \code{add} computes the vector superposition of a collection of
  vectors (supplied as arguments).  Need not be associative (if it
  always was, a two-argument addition function would suffice).  \code{x
  + y} is a synonym for \code{add(x, y)}. \code{x - y} is a synonym for
\code{add(x, vsascale(y, -1))}.  (vector value)
  \item \code{addnorm} computes the vector superposition of a collection of
  vectors (supplied as arguments) and normalizes the result.
  \item \code{vsapower} and its synonym \code{^} compute the power of a vector. Negative powers are allowed (vector value)
  \item \code{equiv} compares two vectors for approximate relative equality (to a quite strict tolerance) (logical value)
  \item \code{cosine} and its synonym \code{\%cos\%} compute the
  normalized similarity of two vectors and return a scalar value with class \code{simval}
  \item \code{dot} and its synonym \code{\%.\%} compute the unnormalized
  similarity of two vectors and return a scalar value with class \code{simval}
  }
}

\value{
  \code{vsaprod}, \code{appinv}, \code{add}, \code{norm}, \code{vsapower}, \code{vsascale} return a \code{vsa} object of the same class as their argument.

  \code{mag}, \code{dot} and \code{cosine} return a scalar value with class \code{simval}
}
\references{ \url{http://www.d-reps.org} }
\author{ Tony Plate \email{tplate@acm.org}}

\seealso{\code{\link{vsa}}, \code{\link{vsa.mem}}, for operator precedence: \code{\link{Syntax}}}
\examples{
library(vsa)
a <- newVec()
a
b <- newVec()
b
dot(a, b)
a %.% b
elts(a)[1:10]
x <- a * b
(!a)
(!a * x) %.% b # not what is expected, because "*" has higher precedence than "!"
((!a) * x) %.% b
(appinv(a) * x) %.% b
(x / a) %.% b
c <- newVec()
d <- newVec()
e <- newVec()
f <- newVec()
x <- norm(a * b + c * d + e * f)
mem <- vsamem(a, b, c, d, e, f)
cosmem(mem, a)
dotmem(mem, x / b)
cosmem(mem, x / b)
cosmem(mem, x / c)
cosmem(mem, x / f)
bestmatch(mem, x / b, cos=TRUE)
bestmatch(mem, x / b, cos=FALSE)
}
\keyword{misc}% at least one, from doc/KEYWORDS

