# ihrr32: spatial HRRs implemented with 32-bit integers
# Normalized real HRRs have an sd of 1/sqrt(n)
# A HRR[1024] has an sd of 1/32.  So, all values
# should be in +/- 8/32.
# A HRR[64] has an sd of 1/8, so all values should
# be in +/- 8/8
# So, the interpretation of 32 bit int x can be as
# x/(2^32) = x/M.  This will produce numbers in the
# right range for any HHR with >= 64 bits (even 32,
# but that's too short to be useful anyway)
# Thus:
#   z/M = x/M * y/M => z = (x * y) / M.
#   z/M = x/M + y/M => z = x + y
#
# Generate random numbers approx normally distributed
# as ((sum of 6 uniform[0-1]) - 3)*sqrt(2)*M
# sd(colSums(matrix(runif(6000), nrow=6))-3)*sqrt(2)
#
# Can rescale to
#   ((sum of 6 uniform[0-16]) - 48)*sqrt(2)*(M/16)

newVec.ihrr32 <- function(what=c("rand", "I", "1", "0", "NA"),
                           len=options("vsalen")[[1]],
                           elts=NULL,
                           cnorm=options("vsacnorm")[[1]],
                           opnorm=options("vsaopnorm")[[1]],
                           vsatype) {
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
    res <- structure(as.vector(elts), class=c("ihrr32", "vsa"))
    if (cnorm && what=="rand")
        res <- norm(res)
    res
}

vsaprod.ihrr32 <- function(e1, e2, method=c("fft", "outer"))
{
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    if (is.null(method) || match.arg(method)=="fft") {
        # fast method using FFT's
        res <- Re(fft(fft(unclass(e1)) * fft(unclass(e2)), inverse=TRUE)) / length(e1)
    } else {
        # very slow method!
        x <- outer(e1, e2)
        cell <- (row(x) + col(x) - 2) %% length(e1) + 1
        res <- unlist(tapply(x, cell, sum, simplify=F), use.names=F)
    }
    class(res) <- class(e1)
    return(res)
}

vsapower.ihrr32 <- function(e1, e2) {
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
        e1[] <- Re(fft(res ^ e2, inverse=TRUE) / length(e1))
        return(e1)
    }
}

appinv.ihrr32 <- function(e1) {
    res <- e1
    res[seq(2,length(e1))] <- res[seq(length(e1),2,by=-1)]
    return(res)
}

dot.ihrr32 <- function(e1, e2) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    structure(sum(unclass(e1) * unclass(e2)), class="simval")
}

equiv.ihrr32 <- function(e1, e2, tol=1e-6) {
    if (class(e1)[1]!=class(e2)[1])
        stop("e1 and e2 must have same class")
    if (length(e1) != length(e2))
        stop("e1 and e2 must have the same length")
    e1 <- elts(e1)
    e2 <- elts(e2)
    return(sum((e1-e2)^2) < tol * max(sum(e1^2), sum(e2^2)))
}

cosine.ihrr32 <- function(e1, e2, mag1=NULL, mag2=NULL) {
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

norm.ihrr32 <- function(e1) {
    return(e1 / sqrt(sum(unclass(e1)^2)))
}

mag.ihrr32 <- function(e1, actual=NULL) {
    # supply actual in the cases where e1 is a template/class-holder
    if (is.null(actual))
        return(sqrt(sum(unclass(e1)^2)))
    else
        return(sqrt(sum(unclass(actual)^2)))
}

add.ihrr32 <- function(e1, ...) {
    res <- elts(e1)
    for (e2 in list(...)) {
        if (class(e1)[1]!=class(e2)[1] || length(e1)!=length(e2))
            stop("all vsa's in a superposition must have same class and length")
        res <- res + elts(e2)
    }
    e1[] <- res
    e1
}

vsascale.ihrr32 <- function(e1, e2) {
    if (is(e2, "vsa") || is(e2, "simval") || !is.numeric(e2) || length(e2)!=1)
        stop("e2 must be a scalar")
    e1[] <- elts(e1) * as.vector(e2)
    e1
}

# The default will work, but a more efficent version can be supplied that
# dispatches off the vsa subclass (i.e., the type of the vsa vector).
# The method can safely assume the columns of mem conform with x.
dotmem.vsamat.compute.ihrr32 <- function(x, mem, ..., cos=FALSE, method=c("fast", "R"), usenames=TRUE) {
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
        if (usenames)
            names(res) <- colnames(mem)
    } else {
        res <- drop(crossprod(x, unclass(mem)))
    }
    if (cos) {
        memmag <- attr(mem, "mag")
        if (length(memmag) != ncol(mem))
            memmag <- sqrt(colSums(unclass(mem)^2, na.rm=TRUE))
        if (any(i <- memmag==0))
            memmag[i] <- 1
        res <- res / (xmag * memmag)
    }
    res
}

ihrr32.dotmem <- function(x, mem, ..., cos=FALSE, method=c("crossprod", "matprod")) {
    method <- match.arg(method)
    if (!is(x, "vsa"))
        stop("x must be a vsa object")
    if (!is(mem, "vsamem"))
        stop("x must be a vsamem object")
    xmag <- mag(x)
    if (xmag==0)
        xmag <- 1
    if (storage.mode(mem)!="double")
        stop("storage.mode(mem)!='double'")
    if (storage.mode(x)!="double")
        stop("storage.mode(x)!='double'")
    res <- numeric(ncol(mem))
    if (method=="matprod")
        .C("matprod_skipna", t(mem), ncol(mem), nrow(mem), x, length(x), 1L, res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    else
        .C("crossprod_skipna", mem, nrow(mem), ncol(mem), x, length(x), 1L, res, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    names(res) <- colnames(mem)
    if (cos) {
        memmag <- attr(mem, "mag")
        if (length(memmag) != ncol(mem))
            memmag <- sqrt(colSums(unclass(mem)^2, na.rm=TRUE))
        if (any(i <- memmag==0))
            memmag[i] <- 1
        res <- res / (xmag * memmag)
    }
    res
}
