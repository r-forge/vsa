# setOldClass(c("realhrr", "vsa")) # the class for vectors

newVec.realhrr <- function(what=c("rand", "I", "1", "0", "NA"),
                           len=options("vsalen")[[1]],
                           elts=NULL,
                           norm=options("vsanorm")[[1]],
                           vsatype) {
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

vsaprod.realhrr <- function(e1, e2, method=c("fft", "outer"))
{
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
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

vsapower.realhrr <- function(e1, e2) {
    if (is(e2, "vsa") || is(e2, "simval") || !is.numeric(e2) || length(e2)!=1)
	stop("e2 must be a scalar")
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

appinv.realhrr <- function(e1) {
    res <- e1
    res[seq(2,length(e1))] <- res[seq(length(e1),2,by=-1)]
    return(res)
}

dot.realhrr <- function(e1, e2) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    structure(sum(unclass(e1) * unclass(e2)), class="simval")
}

equiv.realhrr <- function(e1, e2, tol=1e-6) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    e1 <- elts(e1)
    e2 <- elts(e2)
    return(sum((e1-e2)^2) < tol * max(sum(e1^2), sum(e2^2)))
}

cosine.realhrr <- function(e1, e2) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    e1 <- unclass(e1)
    e2 <- unclass(e2)
    structure(sum(e1 * e2)/sqrt(sum(e1^2)*sum(e2^2)), class="simval")
}

norm.realhrr <- function(e1) {
    return(e1 / sqrt(sum(unclass(e1)^2)))
}

mag.realhrr <- function(e1) {
    return(sqrt(sum(unclass(e1)^2)))
}

add.realhrr <- function(e1, ...) {
    res <- elts(e1)
    for (e2 in list(...)) {
        if (class(e1)[1]!=class(e2)[1] || length(e1)!=length(e2))
            stop("all vsa's in a superposition must have same class and length")
        res <- res + elts(e2)
    }
    e1[] <- res
    e1
}

vsascale.realhrr <- function(e1, e2) {
    if (is(e2, "vsa") || is(e2, "simval") || !is.numeric(e2) || length(e2)!=1)
	stop("e2 must be a scalar")
    e1[] <- elts(e1) * as.vector(e2)
    e1
}
