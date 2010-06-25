\name{vsa.mem}
\alias{vsa.mem}
\alias{vsamem}
\alias{addmem}
\alias{delmem}
\alias{bestmatch}
\alias{cleanup}
\alias{cosmem}
\alias{dotmem}
\alias{conformable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to work with VSA memories. }
\description{

  A VSA memory is a collection (list) of VSA vectors, which are usually
named.  These functions create and modify VSA memories, and compare
VSA vectors against the vectors in a memory.

}
\usage{
vsamem(..., call = match.call())
addmem(..., call = match.call())
delmem(mem, items)
dotmem(mem, x, cos=FALSE)
cosmem(mem, x)
cleanup(mem, x, cos=TRUE, threshold=NA)
bestmatch(mem, x, cos=TRUE, n=1, num=FALSE)
conformable(x=NULL, list, stop.on.error=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ A collection of conformable VSA vectors (i.e., of same
    type and length) to be stored in a VSA memory, or VSA memory objects
    containing such conformable VSA vectors.}
  \item{call}{ The construction call to be recorded with the memory }
  \item{mem}{ A VSA memory object }
  \item{items}{ A character vector containing names of the items to
    delete from the moery}
  \item{x}{ A VSA vector }
  \item{list}{ A list of VSA vectors (can be a VSA memory object) }
  \item{cos}{ \code{TRUE} or \code{FALSE} indicating whether or not to
    use a cosine similarity measure (which is scaled).  The alternative
    is the unscaled dot-product similarity measure.}
  \item{threshold}{ A scalar value }
  \item{n}{ An integer }
  \item{num}{ \code{TRUE} or \code{FALSE}.  If \code{TRUE}, 
    return numeric indices of the vectors in memory.  If \code{FALSE}
    return names of the vectors in memory. }
  \item{stop.on.error}{ Should the function stop with an error if it finds a problem, or return
    \code{FALSE} }
}
\details{
  The functions \code{vsamem()}, \code{addmem()}, and \code{delmem()}
  return the created or modified VSA memory as their value.

}
\value{
  The functions \code{vsamem()}, \code{addmem()}, and \code{delmem()}
  return a \code{vsamem} object.

  The functions \code{dotmem()} and \code{cosmem()} return a
  \code{simval} object.

  The function \code{cleanup()} returns a VSA vector.

  The function \code{conformable()} returns \code{TRUE} if \code{x} and all objects in \code{list}
  are of the same VSA type and have the same length.  Otherwise it return \code{FALSE} or stops with
  an error depending on the value of \code{stop.on.error}.

  The function \code{bestmatch()} returns a numeric
  containing up to the specified number of best matches.  If
  \code{num=FALSE}, the numeric values are the similarity scores, and
  the names are the names of the objects from the \code{mem}.  If
  \code{num=TRUE}, the numeric values are indices of the positions of
  the best matches in \code{mem}, and an attribute \code{scores} is
  attached to it.
}
\author{ Tony Plate \email{tplate@acm.org} }

\seealso{ \code{\link{vsa}} }
\examples{
options(vsatype="realhrr")
options(vsanorm=TRUE)
options(vsalen=1024)
set.seed(0)

circle <- newVec()
square <- newVec()
triangle <- newVec()
star <- newVec()
small <- newVec()
large <- newVec()

large_circle <- addnorm(circle, large)
large_square <- addnorm(square, large)
large_triangle <- addnorm(triangle, large)
large_star <- addnorm(star, large)
small_circle <- addnorm(circle, small)
small_square <- addnorm(square, small)
small_triangle <- addnorm(triangle, small)
small_star <- addnorm(star, small)

sizeMem <- vsamem(small, large)
shapeMem <- vsamem(circle, square, triangle, star, large_circle, large_square, large_triangle,
large_star, small_circle, small_square, small_triangle, small_star)
bestmatch(shapeMem, large_circle, n=3)
bestmatch(sizeMem, large_circle, n=3)

allMem <- addmem(sizeMem, shapeMem)
bestmatch(allMem, large_circle, n=3)
}
\keyword{ misc }% at least one, from doc/KEYWORDS
