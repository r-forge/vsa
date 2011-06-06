\name{realhrr}
\alias{realhrr}
\alias{newVec.realhrr}
\alias{vsaprod.realhrr}
\alias{vsapower.realhrr}
\alias{appinv.realhrr}
\alias{dot.realhrr}
\alias{equiv.realhrr}
\alias{cosine.realhrr}
\alias{norm.realhrr}
\alias{mag.realhrr}
\alias{add.realhrr}
\alias{vsascale.realhrr}
\title{Implementation of Holographic Reduced Representations}
\description{

Specific functions for creating and working with Holographic Reduced
Representations.  Supplies functions for all standard algebraic VSA
operations on vectors: multiplication, addition, superposition, scaling,
powers (including inverses), approximate inverse, and similarity.  These
work in functional and operator form (e.g., \code{add(x, y)} and \code{x
* y} both call the \code{add} method for class \code{realhrr}.

}
\usage{
% Just document the specific methods for realhrr - all the operators
% will use these methods.
\method{newVec}{realhrr}(what = c("rand", "I", "1", "0", "NA"), len = options("vsalen")[[1]],
        elts = NULL, cnorm = getOption("vsacnorm", TRUE), opnorm = getOption("vsaopnorm", FALSE),
        vsatype = getOption("vsatype"))
\method{vsaprod}{realhrr}(e1, e2, method = c("fft", "outer"))
\method{appinv}{realhrr}(e1)
%\special{!e1}
% \method{!}{vsa}(e1)
\method{add}{realhrr}(e1, ...)
\method{norm}{realhrr}(e1)
\method{dot}{realhrr}(e1, e2)
\method{equiv}{realhrr}(e1, e2, tol)
\method{cosine}{realhrr}(e1, e2, mag1, mag2)
\method{mag}{realhrr}(e1, actual)
\method{vsapower}{realhrr}(e1, e2)
\method{vsascale}{realhrr}(e1, e2)
\special{addnorm(...)}
%\special{e1 \%==\% e2} % Rd can't parse this no matter how many
% backslashes I use -- tried between 0 and 7
%\special{e1 \%cos\% e2}
%\special{e1 \%.\% e2}
%\special{e1 * e2}
%\special{e1 + e2}
%\special{e1 / e2}
%\special{e1 - e2}
%\special{e1 ^ e2}
%\method{*}{realhrr}(e1, e2)
%\method{+}{realhrr}(e1, e2)
%\method{/}{realhrr}(e1, e2)
%\method{-}{realhrr}(e1, e2)
%\method{^}{realhrr}(e1, e2)
}
\arguments{
  \item{what}{ the type of vector wanted: one of \code{rand},
    \code{I}, \code{1}, \code{0}, \code{NA}:
    \itemize{\item \code{rand} is a vector filled with random samples
      from a Gaussian distribution with mean 0 and standard deviation \code{1/len}
      \item \code{I} and \code{1} mean the same thing and return the
identity vector (i.e., the vector satisfying \code{x * I = x} for all
non-zero \code{x})
\item \code{0} is the zero vector, i.e., the vector
satisfying \code{x * 0 = 0} and \code{x + 0 = x} for all vectors
\code{x}
\item \code{NA} is a vector filled with \code{NA}
values.} \code{what} is ignored if \code{elts} is supplied.}

  \item{len}{ the number of elements in the vector }

  \item{elts}{ the raw elements for the vector }
  \item{cnorm}{ a logical value (\code{TRUE} or \code{FALSE}) - controls
  whether the new vectors is created to have a magnitude of 1}
  \item{opnorm}{ a logical value (\code{TRUE} or \code{FALSE}) - controls
  whether the result of a vsa operation is normalized to have a magnitude of 1}
  \item{vsatype}{ ignored for \code{newVec.realhrr} (as with all methods of \code{newVec}) }
  \item{e1}{ a \code{vsa} vector -- a subclass of \code{vsa} }
  \item{e2}{ a \code{vsa} vector for \code{vsaprod}, \code{dot}, \code{equiv}, \code{cosine}, and \code{add}, or a scalar (a numeric vector of length 1) for \code{vsapower} and \code{vsascale} }
  \item{mag1,mag2}{ precomputed magnitudes for \code{e1} and \code{e2} }
  \item{actual}{work with data values from this object (in this mode of
    invocation \code{e1} is a template holding the class)}
  \item{method}{ computational method for convolution (\code{vsaprod}) }
  \item{tol}{ numerical tolerance for \code{equiv} }
  \item{...}{ addtional \code{vsa} vectors for \code{add}}
}
\details{

These functions provide implementations of VSA operations for
real-valued Holographic Reduced Representations (in the "spatial"
domain).  All of these functions take as arguments one or more
\code{realhrr} VSA vectors, and possibily a scalar.  The result is a
\code{realhrr} VSA vector, or a scalar with class \code{simval}.  To
catch some errors due to misunderstanding of R's precedence, a \code{simval}
scalar cannot be supplied as the argument of one of these functions (to
do, remove the \code{simval} class using the function \code{scalar}).

  \itemize{
  \item \code{mag} computes the magnitude of a vector (scalar value) as
  the L2 norm (\code{sum(e1^2)})
  \item \code{norm} normalizes a vector so that its magnitude is 1 (vector value)
  \item \code{vsaprod} and its synonym \code{*} compute the circular
  (wrapped) convolution of two vectors (vector value)
  \item \code{appinv} and its synonym \code{!} computes an approximate
  inverse of a vector, which for \code{realhrr} is
  \code{e1[c(1,seq(len(e1), 2)]}.  The expression \code{x / y} where \code{x} and \code{y} are vectors computes \code{x * appinv(y)}.  (vector value)
  \item \code{vsascale} and its synonym \code{*} multiple a vector
  elementwise by a
  scalar (vector value).
  \item \code{add} computes the vector superposition of a collection of
  vectors (supplied as arguments).  For \code{realhrr}s this is simply
  the elementwise sum of its arguments.  \code{x
  + y} is a synonym for \code{add(x, y)}. \code{x - y} is a synonym for
\code{add(x, vsascale(y, -1))}.  (vector value)
  \item \code{addnorm} computes the vector superposition of a collection of
  vectors (supplied as arguments) and normalizes the result.
  \item \code{vsapower} and its synonym \code{^} compute the power of a vector. Negative powers are allowed (vector value)
  \item \code{equiv} compares two vectors for approximate relative equality (to a quite strict tolerance) (logical value)
  \item \code{cosine} and its synonym \code{\%cos\%} compute the
  normalized similarity of two vectors (normalized dot product) (scalar
  value with class \code{simval}).  For \code{realhrr}s, the normalized
  similarity is the cosine of the angle
  between the two vectors.
  \item \code{dot} and its synonym \code{\%.\%} compute the
  un-normalized similarity of two vectors (scalar value with class
  \code{simval}).    For \code{realhrr}s, the unnormalized
  similarity is the dot-product of the two vectors.
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
(!a * x) %.% b # probably not what is expected, because "*" has higher precedence than "!"
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

