# Contains all the methods that must be define for a new type of VSA
# These are all defined as the methods for class "vsa" and give errors
# when they are invoked.

# To define a new subclass of vsa, say "XVSAX", make a copy of this
# file, change the name of each function that ends in ".vsa" to ".XVSAX",
# and give appropriate definitions.  Each of these functions in this
# file has a stop() in it.  After the stop() is code implementing the
# function for "realhrr" -- replace this code by your own code.

# All operators for vsa's are defined in terms of generic functions, so
# operators will work for XVSAX once all the methods in this file are
# defined for XVSAX.

# "newVec" is the one function that uses classed objects in a slightly
# trickly way -- the default method randVec() creates a character object
# with the class of vector we want, and then calls newVec again, which
# will now dispatch to the appropriate newVec method for that class.

newVec.vsa <- function(what=c("rand", "I", "1", "0", "NA"),
                        len=options("vsalen")[[1]],
                        elts=NULL,
                        norm=options("vsanorm")[[1]],
 		        vsatype=options("vsatype")[[1]]) {
    stop("cannot work with 'vsa' vectors -- need to work with a subclass of 'vsa'")
    if (is.null(norm))
        norm <- FALSE
    what <- match.arg(what)
    elts.supplied <- !is.null(elts)
    if (is.null(elts))
	elts <- if (what=="rand") rnorm(len) else double(len)
    if (what=="I" || what=="1" || what=="0" || what=="NA") {
        if (length(elts)!=len)
            elts <- rep(elts, len=len)
        if (what=="NA") {
            elts[] <- NA
        } else {
            elts[] <- 0
            if (what!="0")
                elts[1] <- 1
        }
    }
    if (!is.numeric(elts) || is.array(elts))
        stop("elts must be numeric vector")
    res <- structure(as.vector(elts), class=c("realhrr", "vsa"))
    if (norm && what=="rand")
        res <- norm(res)
    res
}

used.vsa.method <- function(func, cl) paste("need '", func, "' method for '", cl, "'", sep="")

vsaprod.vsa <- function(e1, e2, method=c("fft", "outer"))
{
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    stop(used.vsa.method("vsaprod", class(e1)[1]))
    if (is.null(method) || match.arg(method)=="fft") {
        # fast method using FFT's
        res <- Re(fft(fft(unclass(e1)) * fft(unclass(e2)), inv=T)) / length(e1)
    } else {
        # very slow method!
        x <- outer(e1, e2)
        cell <- (row(x) + col(x) - 2) %% length(e1) + 1
        res <- unlist(tapply(x, cell, sum, simplify=F), use.names=F)
    }
    class(res) <- class(e1)
    return(res)
}

vsapower.vsa <- function(e1, e2) {
    if (is(e2, "vsa") || is(e2, "simval") || !is.numeric(e2) || length(e2)!=1)
	stop("e2 must be a scalar")
    stop(used.vsa.method("vsapower", class(e1)[1]))
    if (e2==1) {
        return(e1)
    } else if (e2==0) {
        e1[] <- 0
        e1[1] <- 1
        return(e1)
    } else {
        res <- fft(unclass(e1))
        e1[] <- Re(fft(res ^ e2, inv=T) / length(e1))
        return(e1)
    }
}

appinv.vsa <- function(e1) {
    stop(used.vsa.method("appinv", class(e1)[1]))
    res <- e1
    res[seq(2,length(e1))] <- res[seq(length(e1),2,by=-1)]
    return(res)
}

dot.vsa <- function(e1, e2) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    stop(used.vsa.method("dot", class(e1)[1]))
    structure(sum(unclass(e1) * unclass(e2)), class="simval")
}

equiv.vsa <- function(e1, e2, tol=1e-6) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    stop(used.vsa.method("equiv", class(e1)[1]))
    e1 <- elts(e1)
    e2 <- elts(e2)
    return(sum((e1-e2)^2) < tol * max(sum(e1^2), sum(e2^2)))
}

cosine.vsa <- function(e1, e2) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    stop(used.vsa.method("cosine", class(e1)[1]))
    e1 <- unclass(e1)
    e2 <- unclass(e2)
    structure(sum(e1 * e2)/sqrt(sum(e1^2)*sum(e2^2)), class="simval")
}

norm.vsa <- function(e1) {
    stop(used.vsa.method("norm", class(e1)[1]))
    return(e1 / sqrt(sum(unclass(e1)^2)))
}

mag.vsa <- function(e1) {
    stop(used.vsa.method("mag", class(e1)[1]))
    return(sqrt(sum(unclass(e1)^2)))
}

add.vsa <- function(e1, ...) {
    stop(used.vsa.method("add", class(e1)[1]))
    res <- elts(e1)
    for (e2 in list(...)) {
        if (class(e1)[1]!=class(e2)[1] || length(e1)!=length(e2))
            stop("all vsa's in a superposition must have same class and length")
        res <- res + elts(e2)
    }
    e1[] <- res
    e1
}

vsascale.vsa <- function(e1, e2) {
    if (is(e2, "vsa") || is(e2, "simval") || !is.numeric(e2) || length(e2)!=1)
	stop("e2 must be a scalar")
    stop(used.vsa.method("vsascale", class(e1)[1]))
    e1[] <- elts(e1) * as.vector(e2)
    e1
}

obselete.conformable.vsa <- function(vsa, vsalist, stop.on.error=T)
{
    stop(used.vsa.method("conformable", class(vsa)[1]))
    if (length(vsalist)==0)
        return(TRUE)
    if (any(sapply(vsalist, class)!=class(vsa)))
        if (stop.on.error)
            stop("all vsas must have the same class")
        else
            return(FALSE)
    if (any(sapply(vsalist, length)!=length(vsa)))
        if (stop.on.error)
            stop("all vsas must have the same length")
        else
            return(FALSE)
    return(TRUE)
}
