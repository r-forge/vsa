\name{vsa}
\alias{VSA}
\alias{vsa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vector Symbolic Architectures}
\description{
  
  Vector Symbolic Architectures (VSAs) are ways of representing complex
  concepts in distributed representations using associative memory
  techniques.  These representations can be used models for cognitive
  science and artificial intelligence.  This package contains data
  structures and functions for implementing VSA schemes.  It currently
  provides an implementation of one VSA scheme: "Holographic Reduced
  Representations".
  
}
  
\details{
  VSAs involve two broad classes of objects:
  \enumerate{
    \item Objects for representing entities, roles, fillers,
    relationships, etc.  These are all represented as uniformly-sized
    one-dimensional arrays (vectors).  The documentation refers to these
    as "VSA objects" or "VSA vectors."
    \item Objects that contain collections of VSA vectors.  These are
    called VSA memories.  The one type of VSA memory implemented here is
    a list of VSA vectors.  The principal functionality provided by a
    VSA memory is returning the closest, or k closest matches to a probe vector.
    }
    
  VSA vectors can be of different types or classes.  The \code{vsa}
  package currently provides an implementation of real-valued
  "Holographic Reduced Representations" (HRR) vectors.  Other possible
  types include binary vectors, frequency-domain HRRs, and more.
  
  The code here uses classes and methods.  This allows different VSA
  schemes to provide methods for common operations.  For example, if
  \code{X} and \code{Y} are VSA objects (vectors), then \code{X * Y}
  will call an appropriate binding method (product), while \code{X + Y}
  will call an appropriate superposition method (addition).  However,
  \code{3 * X} will scale \code{X} by a factor of three.
  
  There are several groups of functions for working with VSAs:
  \itemize{
    \item functions for creating and printing VSA vectors.
    \item binding, unbinding, superposition, and scaling functions
    \item functions for computing the similarity of two VSA vectors.
    \item functions for creating VSA memories, and finding closest matches.
    }
}

\references{
\url{http://d-reps.org/vsapubs} \url{http://www.d-reps.org/tplate}
  
Tony Plate (2003) "Holographic Reduced Representations", CSLI Lecture
Notes Number 150, CSLI Publications, Stanford, CA, ISBN 1575864304.
}
\author{ Tony Plate \email{tplate@acm.org} }

\seealso{ \code{\link{vsa.create}}, \code{\link{vsa.funcs}},
\code{\link{vsa.misc}}, \code{\link{vsa.mem}},
\code{\link{vsamem}}
}
\examples{
# Set the random seed
set.seed(0)
# Generate some HRRs (default is 1024 elements)
a <- newVec()
b <- newVec()
c <- newVec()
d <- newVec()
e <- newVec()
f <- newVec()
# Store these in a memory
mem <- vsamem(a, b, c, d, e, f)
# Find the 3 best matchs in the memory for the vector 'a' (using cosine score)
bestmatch(mem, a, n=3)
# Convolution of 'a' and 'b', prints length, magnitude & cosine with the identity vector
a * b
# Create 'x' as the superposition of 'c' with the binding of 'a' and 'b'
x <- a * b + c
# Find the 3 best matches for x in the memory (expect to see 'c' as the best match)
bestmatch(mem, x, n=3)
# Find what's bound to 'a' in 'x' (!a is the approximate inverse of x, synonym x / a)
bestmatch(mem, x * !a, n=3)
# Find what's bound to 'b' in 'x' (!a is the approximate inverse of x, synonym x / a)
bestmatch(mem, x * !b, n=3)

# Be careful with operator precedence: !a * x parses as !(a * x), not as (!a) * x
bestmatch(mem, !a * x, n=3)
bestmatch(mem, (!a) * x, n=3)
rm(a, b, c, d, e, f, x, mem)

# This code reproduces the experiments described
# in Chapter 6 of the book "Holographic Reduced
# Representations" (the geometric shapes example).
# It reproduces the results in Table 16 on p179 of
# that book. (Note that results in the table vary
# by the random seed used, but the results generated
# here should usually be with in ranges in the
# square brackets.)

options(vsatype="realhrr")
options(vsanorm=TRUE)
options(vsalen=1024)
set.seed(0)

# weights for different components of HRRs
w_rpl <- 1 # relation-predicate label
w_rpb <- 1 # relation-predicate binding
w_cow <- 1 # object weight in chunk
w_cpw <- 1 # relation-predicate weight in chunk

# objects.ca
# Experiment from Markman, Gentner, & Wisniewski, draft
# Stimuli from Figure 2, discussed p12

horizontal <- newVec()
leftR <- newVec()
rightR <- newVec()

vertical <- newVec()
aboveR <- newVec()
belowR <- newVec()

circle <- newVec()
square <- newVec()
triangle <- newVec()
star <- newVec()
small <- newVec()
large <- newVec()

# shapes.ca
large_circle <- addnorm(circle, large)
large_square <- addnorm(square, large)
large_triangle <- addnorm(triangle, large)
large_star <- addnorm(star, large)
small_circle <- addnorm(circle, small)
small_square <- addnorm(square, small)
small_triangle <- addnorm(triangle, small)
small_star <- addnorm(star, small)

# lcalq : large circle above large square
lcalq_rf <- w_rpb * aboveR * large_circle + w_rpb * belowR * large_square
lcalq <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_circle, large_square),
                 w_cpw * lcalq_rf)

# lcalt : large circle above large triangle
lcalt_rf <- w_rpb * aboveR * large_circle + w_rpb * belowR * large_triangle
lcalt <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_circle, large_triangle),
                 w_cpw * lcalt_rf)

# lcnlt : large circle nextto large triangle
lcnlt_rf <- w_rpb * leftR * large_circle + w_rpb * rightR * large_triangle
lcnlt <- addnorm(w_rpl * horizontal,
                 w_cow * addnorm(large_circle, large_triangle),
                 w_cpw * lcnlt_rf)

# lqalc : large square above large circle
lqalc_rf <- w_rpb * aboveR * large_square + w_rpb * belowR * large_circle
lqalc <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_square, large_circle),
                 w_cpw * lqalc_rf)

# lqnls : large square nextto large star
lqnls_rf <- w_rpb * leftR * large_square + w_rpb * rightR * large_star
lqnls <- addnorm(w_rpl * horizontal,
                 w_cow * addnorm(large_square, large_star),
                 w_cpw * lqnls_rf)

# lsalq : large star above large square
lsalq_rf <- w_rpb * aboveR * large_star + w_rpb * belowR * large_square
lsalq <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_star, large_square),
                 w_cpw * lsalq_rf)

# lsalt : large star above large triangle
lsalt_rf <- w_rpb * aboveR * large_star + w_rpb * belowR * large_triangle
lsalt <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_star, large_triangle),
                 w_cpw * lsalt_rf)

# lsasq : large star above small square
lsasq_rf <- w_rpb * aboveR * large_star + w_rpb * belowR * small_square
lsasq <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_star, small_square),
                 w_cpw * lsasq_rf)

# lsnlq : large star nextto large square
lsnlq_rf <- w_rpb * leftR * large_star + w_rpb * rightR * large_square
lsnlq <- addnorm(w_rpl * horizontal,
                 w_cow * addnorm(large_star, large_square),
                 w_cpw * lsnlq_rf)

# ltalc : large triangle above large circle
ltalc_rf <- w_rpb * aboveR * large_triangle + w_rpb * belowR * large_circle
ltalc <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_triangle, large_circle),
                 w_cpw * ltalc_rf)

# ltals : large triangle above large star
ltals_rf <- w_rpb * aboveR * large_triangle + w_rpb * belowR * large_star
ltals <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_triangle, large_star),
                 w_cpw * ltals_rf)

# ltasc : large triangle above small circle
ltasc_rf <- w_rpb * aboveR * large_triangle + w_rpb * belowR * small_circle
ltasc <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(large_triangle, small_circle),
                 w_cpw * ltasc_rf)

# ltnlc : large triangle nextto large circle
ltnlc_rf <- w_rpb * leftR * large_triangle + w_rpb * rightR * large_circle
ltnlc <- addnorm(w_rpl * horizontal,
                 w_cow * addnorm(large_triangle, large_circle),
                 w_cpw * ltnlc_rf)

# ssalq : small star above large square
ssalq_rf <- w_rpb * aboveR * small_star + w_rpb * belowR * large_square
ssalq <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(small_star, large_square),
                 w_cpw * ssalq_rf)

# stalc : small triangle above large circle
stalc_rf <- w_rpb * aboveR * small_triangle + w_rpb * belowR * large_circle
stalc <- addnorm(w_rpl * vertical,
                 w_cow * addnorm(small_triangle, large_circle),
                 w_cpw * stalc_rf)

relMem <- vsamem(vertical, horizontal)
roleMem <- vsamem(leftR, rightR, aboveR, belowR)
sizeMem <- vsamem(small, large)
shapeMem <- vsamem(circle, square, triangle, star,
                   large_circle, large_square, large_triangle, large_star,
                   small_circle, small_square, small_triangle, small_star)

arrMem <- vsamem(ltnlc, lqnls, lcalt, lcnlt, lsalq, lsnlq,
                 ltalc, lqalc, lcalq, ltals, lsalt, stalc,
                 ltasc, ssalq, lsasq)

# res1.ca

options(digits=3)
makeTable16 <- function() {
    x <- ltalc
    cat("Comparing with ltalc", "\n")
    C1A <- x %.% ltnlc
    C1B <- x %.% lqnls
    cat("1  ltnlc: ", format(C1A, wid=5),
        "    lqnls: ", format(C1B, wid=5), "\n")
    C2aA <- x %.% lcalt
    C2aB <- x %.% lcnlt
    cat("2a lcalt: ", format(C2aA, wid=5),
        "    lcnlt: ", format(C2aB, wid=5), "\n")
    C2bA <- x %.% lsalq
    C2bB <- x %.% lsnlq
    cat("2b lsalq: ", format(C2bA, wid=5),
        "    lsnlq: ", format(C2bB, wid=5), "\n")
    C3A <- x %.% ltalc
    C3B <- x %.% lcalt
    cat("3  ltalc: ", format(C3A, wid=5),
        "    lcalt: ", format(C3B, wid=5), "\n")
    C4aA <- x %.% lqalc
    C4aB <- x %.% lcalq
    cat("4a lqalc: ", format(C4aA, wid=5),
        "    lcalq: ", format(C4aB, wid=5), "\n")
    C4bA <- x %.% ltals
    C4bB <- x %.% lsalt
    cat("4b ltals: ", format(C4bA, wid=5),
        "    lsalt: ", format(C4bB, wid=5), "\n")
    x <- stalc
    cat("Comparing with stalc", "\n")
    C5A <- x %.% stalc
    C5B <- x %.% ltasc
    cat("5  stalc: ", format(C5A, wid=5),
        "    ltasc: ", format(C5B, wid=5), "\n")
    C6A <- x %.% ssalq
    C6B <- x %.% lsasq
    cat("6  ssalq: ", format(C6A, wid=5),
        "    lsasq: ", format(C6B, wid=5), "\n")
}
makeTable16()

# most similar arrangments to "large triangle above large circle"
bestmatch(arrMem, ltalc, n=3)
# the relation in "large triangle above large circle"
bestmatch(relMem, ltalc, n=3)
# the object in role 'above' in "large triangle above large circle"
bestmatch(shapeMem, ltalc * appinv(aboveR), n=3)
# size of object 'above' in "large triangle above large circle"
bestmatch(sizeMem, ltalc * appinv(aboveR), n=3)
# the object in role 'below' in "large triangle above large circle"
bestmatch(shapeMem, ltalc * appinv(belowR), n=3)
# what size is the object 'below' in "large triangle above large circle"
bestmatch(sizeMem, ltalc * appinv(belowR), n=3)
}
\keyword{misc}% at least one, from doc/KEYWORDS
