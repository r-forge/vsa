# Structure of a realhrr ramdb:
# list with components:
#     db: External ptr
#     mag: External ptr
#     veclen: integer
#     memsize: integer
#     labels: character
# and attributes:
#     example:
#     * (no longer used) mag
# Use attributes for example and mag to be same as other vsamem classes
addmem.vsadb <- function(..., mem=NULL, memsize=NULL, fill=FALSE, labels=NULL, call=match.call(expand.dots=FALSE), datatype=NULL, lookup=NULL) {

    call$...$memsize <- NULL
    call$...$mem <- NULL
    call$...$fill <- NULL
    call$...$datatype <- NULL
    call$...$lookup <- NULL
    if (!is.null(mem) && !is.null(memsize))
        stop("cannot supply both mem and memsize")
    if (length(call$...)) {
        vsamat <- addmem.vsamat(..., labels=labels, call=call)
        example <- attr(vsamat, "example")
        vsamat.memsize <- memsize(vsamat)
        vsamat.labels <- memlabels(vsamat)
    } else {
        example <- newVec("1")
        vsamat.memsize <- 0
        vsamat.labels <- labels
    }
    if (is.null(mem)) {
        memsize <- as.integer(max(memsize, vsamat.memsize))
        if (length(vsamat.labels) < memsize)
            length(vsamat.labels) <- memsize
        mem <- create.realhrr_ramdb(memsize, example=example, fill=fill, labels=vsamat.labels, datatype=datatype, lookup=lookup)
        n <- as.integer(mem$veclen)
        for (i in as.integer(seq(length.out=vsamat.memsize))) {
            mag <- .Call(mem$setter, mem$ramdb, n, memsize, i, vsamat[,i], PACKAGE="vsa")
            .Call("realvec_set", mem$mag, i, mag, PACKAGE="vsa")
        }
    } else {
        if (fill)
            stop("can only supply fill=TRUE when creating a new memory")
        conformable(example, attr(mem, "example"))
        memsize <- as.integer(memsize(mem))
        n <- as.integer(mem$veclen)
        if (vsamat.memsize > memsize)
            stop("not enough room in mem (total ", memsize, ") for vsamat (", vsamat.memsize, ")")
        memlabels <- memlabels(mem)
        ii <- match(vsamat.labels, memlabels)
        if (sum(is.na(ii)) > sum(is.na(memlabels)))
            stop("not enough room in mem (have ", sum(is.na(memlabels)),
                 " empty slots) for vsamat (have ", sum(is.na(ii)), " new labels)")
        if (any(is.na(ii))) {
            jj <- which(is.na(memlabels))
            if (length(jj) > sum(is.na(ii)))
                length(jj) <- sum(is.na(ii))
            memlabels[jj] <- vsamat.labels[is.na(ii)]
            mem$labels <- memlabels
            ii[is.na(ii)] <- jj
        }
        for (i in as.integer(ii)) {
            mag <- .Call(mem$setter, mem$ramdb, n, memsize, i, vsamat[,i], PACKAGE="vsa")
            .Call("realvec_set", mem$mag, i, mag, PACKAGE="vsa")
        }
    }
    return(mem)
}

contents.realhrr_ramdb <- function(mem, ii=seq(length.out=memsize(mem)), ...) {
    veclen <- as.integer(mem$veclen)
    memsize <- as.integer(memsize(mem))
    x <- sapply(as.integer(ii), function(i)
                .Call(mem$getter, mem$ramdb, veclen, memsize, i, PACKAGE="vsa"))
    colnames(x) <- memlabels(mem)[ii]
    x
}

create.realhrr_ramdb <- function(memsize, example=newVec(), fill=FALSE, labels=NULL,
                                 datatype=c("float", "double", "int16", "int8", "int4", "int2",
                                 "oint16", "oint8", "oint4", "oint2", "oint1", "bit"),
                                 lookup=FALSE) {
    if (is.null(lookup)) lookup <- FALSE
    if (!inherits(example, "vsa"))
        stop("'example' must be a vsa object")
    if (!is.null(labels) && (length(labels)!=memsize || !is.character(labels)))
        stop("labels must be character vector with length = memsize")
    if (is.null(datatype)) datatype <- "float"
    datatype <- match.arg(datatype)
    dt <- switch(datatype, float="f", double="d", int16="i16", int8="i8", int4="i4", int2="i2",
                 xint16="x16", xint8="x8", xint4="x4", xint2="x2", bit=, xint1="x1")
    creator <- paste("realhrr_ramdb", dt, "_create", sep="")
    rander <- paste("realhrr_ramdb", dt, "_set_rand", sep="")
    mem <- list(ramdb=.Call(creator, length(example), as.integer(memsize), NULL, NULL, PACKAGE="vsa"))
    mem$datatype <- datatype
    # If the _getraw function doesn't exist getmemraw() won't work, but otherwise no matter
    mem$rawgetter <- paste("realhrr_ramdb", dt, "_getraw", sep="")
    mem$rawsetter <- paste("realhrr_ramdb", dt, "_setraw", sep="")
    mem$getter <- paste("realhrr_ramdb", dt, "_get", sep="")
    mem$setter <- paste("realhrr_ramdb", dt, "_set", sep="")
    mem$dotter <- paste("realhrr_ramdb", dt, "_dot", sep="")
    mem$veclen <- length(example)
    mem$memsize <- as.integer(memsize)
    cnorm <- as.logical(getOption("vsacnorm", TRUE))
    mem$mag <- .Call("realvec_create", mem$memsize, PACKAGE="vsa")
    if (fill) {
        mag <- .Call(rander, mem$ramdb, length(example), as.integer(memsize), NULL, cnorm, NULL, PACKAGE="vsa")
        .Call("realvec_set", mem$mag, seq.int(mem$memsize), mag, PACKAGE="vsa")
    }
    attr(mem, "example") <- example
    # NA's in labels show empty locations
    mem$labels <- if (!is.null(labels)) labels else rep(as.character(NA), memsize)
    if (lookup && datatype=="int8") {
        ii <- c(seq.int(0L, 127L), seq(-128L, -1L))
        mem$mlookup <- .Call("intvec_create", 256L * 256L, PACKAGE="vsa")
        .Call("intvec_set", mem$mlookup, seq.int(length.out=256^2), as.integer(outer(ii, ii, "*")), PACKAGE="vsa")
    }
    class(mem) <- c("realhrr_ramdb", "vsamem")
    mem
}

memlabels.realhrr_ramdb <- function(mem) return(mem$labels)
memsize.realhrr_ramdb <- function(mem) return(mem$memsize)
memmags.realhrr_ramdb <- function(mem) {
   .Call("realvec_get", mem$mag, seq.int(mem$memsize), PACKAGE="vsa")
}

print.realhrr_ramdb <- function(x, ...) {
    vsaclass <- "realhrr"
    len <- x$veclen
    m <- x$labels
    n <- x$memsize
    cat("RAM DB (", x$datatype, ") containing ", n, " ", vsaclass[1], "[", len, "]: ",
        if (is.null(m)) "(no labels)\n"
        else paste(paste(m[seq(1, length.out=min(length(m), 3))], collapse=", "),
                   if (length(m)>3) " ...", "\n", sep=""), sep="")
    invisible(x)
}

getmem.realhrr_ramdb <- function(mem, i) {
    example <- attr(mem, "example")
    if (is.character(i) && length(i)==1) {
        j <- match(i, memlabels(mem))
        if (is.na(j))
            stop("label '", i, "' not found")
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>memsize(mem))
            stop("index ", i, " out of range (mem has ", memsize(mem), " elements)")
    } else {
        stop("i must be single character or numeric")
    }
    memsize <- mem$memsize
    vec <- .Call(mem$getter, mem$ramdb, length(example), as.integer(memsize), as.integer(j), PACKAGE="vsa")
    example[] <- vec
    return(example)
}

getmemraw.realhrr_ramdb <- function(mem, i) {
    example <- attr(mem, "example")
    if (is.character(i) && length(i)==1) {
        j <- match(i, memlabels(mem))
        if (is.na(j))
            stop("label '", i, "' not found")
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>memsize(mem))
            stop("index ", i, " out of range (mem has ", memsize(mem), " elements)")
    } else {
        stop("i must be single character or numeric")
    }
    memsize <- mem$memsize
    vec <- .Call(mem$rawgetter, mem$ramdb, length(example), as.integer(memsize), as.integer(j), PACKAGE="vsa")
    return(vec)
}

setmemraw.realhrr_ramdb <- function(mem, i, x) {
    example <- attr(mem, "example")
    if (is.character(i) && length(i)==1) {
        j <- match(i, memlabels(mem))
        if (is.na(j))
            stop("label '", i, "' not found")
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>memsize(mem))
            stop("index ", i, " out of range (mem has ", memsize(mem), " elements)")
    } else {
        stop("i must be single character or numeric")
    }
    memsize <- mem$memsize
    vec <- .Call(mem$rawsetter, mem$ramdb, length(example), as.integer(memsize), as.integer(j), as.raw(x), PACKAGE="vsa")
    return(vec)
}

setmem.realhrr_ramdb <- function(mem, i, x, label=NULL) {
    example <- attr(mem, "example")
    conformable(example, list(x))
    if (is.character(i) && length(i)==1) {
        j <- match(i, memlabels(mem))
        if (is.na(j)) {
            j <- which(is.na(memlabels(mem)))[1]
            if (is.na(j))
                stop("label '", i, "' not found, and no empty labels to use")
            mem$labels[j] <- i
        }
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>memsize(mem))
            stop("index ", i, " out of range (mem has ", memsize(mem), " elements)")
        if (!is.null(label))
            mem$labels[j] <- label
    } else {
        stop("i must be single character or numeric")
    }
    mag <- .Call(mem$setter, mem$ramdb, length(example), as.integer(memsize(mem)), as.integer(j), x, PACKAGE="vsa")
    .Call("realvec_set", mem$mag, as.integer(j), mag, PACKAGE="vsa")
    return(mem)
}

dotmem.realhrr_ramdb <- function(mem, x, ..., cos=FALSE, debug=FALSE) {
    params <- NULL
    example <- attr(mem, "example")
    conformable(example, list(x))
    ii <- which(!is.na(mem$labels))
    memsize <- ii[length(ii)]
    if (debug)
        params <- as.numeric(c(0, debug))
    dot <- .Call(mem$dotter, mem$ramdb, length(x), as.integer(memsize), ii, x, params, mem$mlookup, PACKAGE="vsa")
    magii <- .Call("realvec_get", mem$mag, ii, PACKAGE="vsa")
    if (cos)
        dot <- dot / (mag(x) * magii)
    names(dot) <- mem$labels[ii]
    dot
}
