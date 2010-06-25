\name{vsa.create}
\alias{vsa.create}
\alias{newVec}
\alias{elts}
\alias{print.vsa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create and print VSA objects}
\description{
  Functions for creating, extracting data from, and printing VSA objects.
}
\usage{
newVec(what = c("rand", "I", "1", "0", "NA"), len = options("vsalen")[[1]],
        elts = NULL, norm = options("vsanorm")[[1]], vsatype = options("vsatype")[[1]])
elts(e1)
\method{print}{vsa}(x, \dots, values = FALSE)
}
\arguments{

  \item{what}{ the type of vector wanted: one of \code{rand},
    \code{I}, \code{1}, \code{0}, \code{NA}:
    \itemize{\item \code{rand} is a random
      vector
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

  \item{vsatype}{ the subclass of vsa: must be a character string for
  which a function \code{newVec.<vsatype>} exists}

  \item{norm}{ a logical value (\code{TRUE} or \code{FALSE}) - controlls
  whether the new vectors is created to have a magnitude of 1}

  \item{e1}{a vsa vector}

  \item{x}{a vsa vector}

  \item{\dots}{additional arguments for \code{print}}

  \item{values}{ a logical value indicating whether or not to print the actual values elements of the vsa vector}
}
\details{
  VSA objects can be of different types or classes.  The \code{vsa}
  package currently provides an implementation of real-valued HRR
  vectors.  Other possible types include binary vectors,
  frequency-domain HRRs, and more.
  
  The length and type of the vector created by \code{newVec()} is
  determined by the current settings for \code{options("vsatype")},
  \code{options("vsanorm")}, and \code{options("vsalen")}.
}

\value{
  \code{newVec} generates a new vsa vector.

  \code{elts} returns the elements of a vsa vector as the most suitable R atomic class (usually double, but other possibilities are integer, raw (bytes), or complex).
}

\author{ Tony Plate \email{tplate@acm.org} }

\seealso{ \code{\link{vsa.funcs}},
\code{\link{vsa.misc}},
\code{\link{vsa.mem}}
}
\examples{
a <- newVec()
b <- newVec()
c <- newVec()
d <- newVec()
c <- a * b
}
\keyword{misc}% at least one, from doc/KEYWORDS
