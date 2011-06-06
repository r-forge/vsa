addmem.vsadb <- function(..., mem=NULL, memsize=NULL, labels=NULL, call=match.call(expand.dots=FALSE)) {
    if (!is.null(mem) && !is.null(memsize))
        stop("cannot supply both mem and memsize")
    vsamat <- addmem.vsamat(..., labels=labels, call=call)
    if (is.null(mem)) {
        memsize <- as.integer(min(memsize, memsize(vsamat)))
        labels <- memlabels(vsamat)
        if (length(labels) < memsize)
            length(labels) <- memsize
        mem <- create.realhrr_ramdb(memsize, example=attr(vsamat, "example"), fill=FALSE, labels=labels)
        n <- as.integer(attr(mem, "veclen"))
        for (i in as.integer(seq(len=memsize(vsamat))))
            attr(mem, "mag")[i] <- .Call("realhrr_ramdb_set", mem, n, memsize, i, vsamat[,i], PACKAGE="vsa")
    } else {
        memsize <- as.integer(memsize(mem))
        n <- as.integer(attr(mem, "veclen"))
        if (memsize(vsamat) > memsize)
            stop("not enough room in mem (total ", memsize, ") for vsamat (", memsize(vsamat), ")")
        matlabels <- memlabels(vsamat)
        memlabels <- memlabels(mem)
        ii <- match(matlabels, memlabels)
        if (sum(is.na(ii)) > sum(is.na(memlabels)))
            stop("not enough room in mem (room for ", sum(is.na(memlabels)), ") for vsamat (have ", sum(is.na(ii)), " new)")
        if (any(is.na(ii))) {
            jj <- which(is.na(memlabels))
            if (length(jj) > sum(is.na(ii)))
                length(jj) <- sum(is.na(ii))
            memlabels[jj] <- matlabels[is.na(ii)]
            attr(mem, "labels") <- memlabels
            ii[is.na(ii)] <- jj
        }
        for (i in as.integer(ii))
            attr(mem, "mag")[i] <- .Call("realhrr_ramdb_set", mem, n, memsize, i, vsamat[,i], PACKAGE="vsa")
    }
    return(mem)
}

contents.realhrr_ramdb <- function(mem, ii=seq(len=memsize(mem))) {
    veclen <- as.integer(attr(mem, "veclen"))
    memsize <- as.integer(memsize(mem))
    x <- sapply(as.integer(ii), function(i)
                .Call("realhrr_ramdb_get", mem, veclen, memsize, i, PACKAGE="vsa"))
    colnames(x) <- memlabels(mem)[ii]
    x
}

create.realhrr_ramdb <- function(memsize, example=newVec(), fill=FALSE, labels=NULL) {
    if (!inherits(example, "vsa"))
        stop("'example' must be a vsa object")
    if (!is.null(labels) && (length(labels)!=memsize || !is.character(labels)))
        stop("labels must be character vector with length = memsize")
    mem <- .Call("realhrr_ramdb_create", length(example), as.integer(memsize), NULL, NULL, PACKAGE="vsa")
    attr(mem, "veclen") <- length(example)
    attr(mem, "memsize") <- as.integer(memsize)
    cnorm <- as.logical(getOption("vsacnorm", TRUE))
    if (fill)
        attr(mem, "mag") <- .Call("realhrr_ramdb_setrand", mem, length(example), as.integer(memsize), NULL, cnorm, NULL, PACKAGE="vsa")
    else
        attr(mem, "mag") <- numeric(memsize)
    attr(mem, "example") <- example
    # NA's in labels show empty locations
    attr(mem, "labels") <- if (!is.null(labels)) labels else rep(as.character(NA), memsize)
    class(mem) <- c("realhrr_ramdb", "vsamem")
    mem
}

memlabels.realhrr_ramdb <- function(mem) return(attr(mem, "labels"))
memsize.realhrr_ramdb <- function(mem) return(attr(mem, "memsize"))

print.realhrr_ramdb <- function(x, ...) {
    vsaclass <- "realhrr"
    len <- attr(x, "veclen")
    m <- attr(x, "labels")
    n <- attr(x, "memsize")
    cat("RAM DB (single precision) containing ", n, " ", vsaclass[1], "[", len, "]: ",
        if (is.null(m)) "(no labels)\n"
        else paste(paste(m[seq(1, len=min(length(m), 3))], collapse=", "),
                   if (length(m)>3) " ...", "\n", sep=""))
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
    memsize <- attr(mem, "memsize")
    vec <- .Call("realhrr_ramdb_get", mem, length(example), as.integer(memsize), as.integer(j), PACKAGE="vsa")
    example[] <- vec
    return(example)
}

setmem.realhrr_ramdb <- function(mem, i, x, label=NULL) {
    example <- attr(mem, "example")
    conformable(x, example)
    if (is.character(i) && length(i)==1) {
        j <- match(i, memlabels(mem))
        if (is.na(j)) {
            j <- which(is.na(memlabels(mem)))[1]
            if (is.na(j))
                stop("label '", i, "' not found, and no empty labels to use")
            attr(mem, "memlabels")[j] <- i
        }
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>memsize(mem))
            stop("index ", i, " out of range (mem has ", memsize(mem), " elements)")
        if (!is.null(label))
            attr(mem, "memlabels")[j] <- label
    } else {
        stop("i must be single character or numeric")
    }
    mag <- .Call("realhrr_ramdb_set", mem, length(example), as.integer(memsize), as.integer(j), x, PACKAGE="vsa")
    attr(mem, "mag")[j] <- mag
    return(j)
}
