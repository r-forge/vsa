# setOldClass(c("realhrr", "vsa")) # the class for vectors

newVec.realhrr <- function(what=c("rand", "I", "1", "0", "NA"),
                           len=options("vsalen")[[1]],
                           elts=NULL,
                           cnorm=getOption("vsacnorm", TRUE),
                           opnorm=getOption("vsaopnorm", FALSE),
                           vsatype=getOption("vsatype")) {
    if (is.null(cnorm))
        cnorm <- FALSE
    what <- match.arg(what)
    elts.supplied <- !is.null(elts)
    if (is.null(elts))
        elts <- if (what=="rand") rnorm(len)*(1/sqrt(len)) else double(len)
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
    if (cnorm && what=="rand")
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

cosine.realhrr <- function(e1, e2, mag1=NULL, mag2=NULL) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    e1 <- unclass(e1)
    e2 <- unclass(e2)
    if (is.null(mag1))
        mag1 <- sqrt(sum(e1^2))
    if (is.null(mag2))
        mag2 <- sqrt(sum(e2^2))
    structure(sum(e1 * e2)/(mag1 * mag2), class="simval")
}

norm.realhrr <- function(e1) {
    return(e1 / sqrt(sum(unclass(e1)^2)))
}

mag.realhrr <- function(e1, actual=NULL) {
    # supply actual in the cases where e1 is a template/class-holder
    if (is.null(actual))
        return(sqrt(sum(unclass(e1)^2)))
    else
        return(sqrt(sum(unclass(actual)^2)))
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

# This cumbersomely-named method-method will get called by
# dotmem(mem, x), cosmem(mem, x), bestmatch(mem, x), cleanup(mem, x)
# via a method call to dotmem.vsamat() which calls the generic
# dotmem.vsamat.compute().
# The default will work, but a more efficent version can be supplied that
# dispatches off the vsa subclass (i.e., the type of the vsa vector).
# The method can safely assume the columns of mem conform with x.
dotmem.vsamat.compute.realhrr <-
    function(x, mem, ..., cos=FALSE,
             method=c("fast", "R", "crossprod", "matprod", "dotmem", "dotmempp", "dotmempp1", "dotmemmp1", "dotmemmp2"),
             cores=FALSE, usenames=TRUE) {
    method <- match.arg(method)
    if (!is(x, "vsa"))
        stop("x must be a vsa object")
    if (!is(mem, "vsamem"))
        stop("x must be a vsamem object")
    if (storage.mode(mem)!="double")
        stop("storage.mode(mem)!='double'")
    if (storage.mode(x)!="double")
        stop("storage.mode(x)!='double'")
    if (length(x) != nrow(mem))
        stop("length(x) != nrow(mem)")
    xmag <- mag(x)
    if (xmag==0)
        xmag <- 1
    if (method=="fast") {
        res <- numeric(ncol(mem))
        .C("crossprod_skipna", mem, nrow(mem), ncol(mem), x, length(x), 1L, res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    } else if (method %in% c("R", "crossprod")) {
        res <- drop(crossprod(x, unclass(mem)))
    } else if (method %in% c("matprod")) {
        res <- numeric(ncol(mem))
        .C("matprod_skipna", t(mem), ncol(mem), nrow(mem), x, length(x), 1L, res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    } else if (method %in% c("dotmem")) {
        res <- numeric(ncol(mem))
        .C("realhrr_dotmem", x, length(x), mem, ncol(mem), res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    } else if (method %in% c("dotmempp")) {
        res <- numeric(ncol(mem))
        .C("realhrr_dotmempp", x, length(x), mem, ncol(mem), res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    } else if (method %in% c("dotmempp1")) {
        res <- numeric(ncol(mem))
        .C("realhrr_dotmempp1", x, length(x), mem, ncol(mem), res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    } else if (method %in% c("dotmemmp1")) {
        res <- numeric(ncol(mem))
        coreused <- integer(8)
        .C("realhrr_dotmemmp1", x, length(x), mem, ncol(mem), res, coreused, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
        if (cores)
            cat("coreused:", paste(coreused, collapse=" "), "\n")
    } else if (method %in% c("dotmemmp2")) {
        res <- numeric(ncol(mem))
        coreused <- integer(8)
        .C("realhrr_dotmemmp2", x, length(x), mem, ncol(mem), res, coreused, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
        if (cores)
            cat("coreused:", paste(coreused, collapse=" "), "\n")
    } else {
        stop("unknown method:", method)
    }
    if (usenames && is.null(names(res)))
        names(res) <- colnames(mem)
    if (cos) {
        memmag <- attr(mem, "mag")
        if (length(memmag) != ncol(mem))
            memmag <- sqrt(colSums(unclass(mem)^2, na.rm=T))
        if (any(i <- memmag==0))
            memmag[i] <- 1
        res <- res / (xmag * memmag)
    }
    res
}

