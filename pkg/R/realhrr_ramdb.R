# Structure of a realhrr ramdb:
# list with components:
#     db: External ptr
#     veclen: integer
#     memsize: integer
#     labels: character
# and attributes:
#     example:
#     mag
# Use attributes for example and mag to be same as other vsamem classes
addmem.vsadb <- function(..., mem=NULL, memsize=NULL, fill=FALSE, labels=NULL, call=match.call(expand.dots=FALSE)) {
    call$...$memsize <- NULL
    call$...$mem <- NULL
    call$...$fill <- NULL
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
        mem <- create.realhrr_ramdb(memsize, example=example, fill=fill, labels=vsamat.labels)
        n <- as.integer(mem$veclen)
        for (i in as.integer(seq(len=vsamat.memsize)))
            attr(mem, "mag")[i] <- .Call("realhrr_ramdb_set", mem$ramdb, n, memsize, i, vsamat[,i], PACKAGE="vsa")
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
        for (i in as.integer(ii))
            attr(mem, "mag")[i] <- .Call("realhrr_ramdb_set", mem$ramdb, n, memsize, i, vsamat[,i], PACKAGE="vsa")
    }
    return(mem)
}

contents.realhrr_ramdb <- function(mem, ii=seq(len=memsize(mem))) {
    veclen <- as.integer(mem$veclen)
    memsize <- as.integer(memsize(mem))
    x <- sapply(as.integer(ii), function(i)
                .Call("realhrr_ramdb_get", mem$ramdb, veclen, memsize, i, PACKAGE="vsa"))
    colnames(x) <- memlabels(mem)[ii]
    x
}

create.realhrr_ramdb <- function(memsize, example=newVec(), fill=FALSE, labels=NULL) {
    if (!inherits(example, "vsa"))
        stop("'example' must be a vsa object")
    if (!is.null(labels) && (length(labels)!=memsize || !is.character(labels)))
        stop("labels must be character vector with length = memsize")
    mem <- list(ramdb=.Call("realhrr_ramdb_create", length(example), as.integer(memsize), NULL, NULL, PACKAGE="vsa"))
    mem$veclen <- length(example)
    mem$memsize <- as.integer(memsize)
    cnorm <- as.logical(getOption("vsacnorm", TRUE))
    if (fill)
        attr(mem, "mag") <- .Call("realhrr_ramdb_set_rand", mem$ramdb, length(example), as.integer(memsize), NULL, cnorm, NULL, PACKAGE="vsa")
    else
        attr(mem, "mag") <- numeric(memsize)
    attr(mem, "example") <- example
    # NA's in labels show empty locations
    mem$labels <- if (!is.null(labels)) labels else rep(as.character(NA), memsize)
    class(mem) <- c("realhrr_ramdb", "vsamem")
    mem
}

memlabels.realhrr_ramdb <- function(mem) return(mem$labels)
memsize.realhrr_ramdb <- function(mem) return(mem$memsize)

print.realhrr_ramdb <- function(x, ...) {
    vsaclass <- "realhrr"
    len <- x$veclen
    m <- x$labels
    n <- x$memsize
    cat("RAM DB (single precision) containing ", n, " ", vsaclass[1], "[", len, "]: ",
        if (is.null(m)) "(no labels)\n"
        else paste(paste(m[seq(1, len=min(length(m), 3))], collapse=", "),
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
    vec <- .Call("realhrr_ramdb_get", mem$ramdb, length(example), as.integer(memsize), as.integer(j), PACKAGE="vsa")
    example[] <- vec
    return(example)
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
    mag <- .Call("realhrr_ramdb_set", mem$ramdb, length(example), as.integer(memsize(mem)), as.integer(j), x, PACKAGE="vsa")
    attr(mem, "mag")[j] <- mag
    return(mem)
}

dotmem.realhrr_ramdb <- function(mem, x, ..., cos=FALSE) {
    conformable(attr(mem, "example"), list(x))
    ii <- which(!is.na(mem$labels))
    memsize <- ii[length(ii)]
    dot <- .Call("realhrr_ramdb_dot", mem$ramdb, length(x), as.integer(memsize), ii, x, PACKAGE="vsa")
    if (cos)
        dot <- dot / (mag(x) * attr(mem, "mag")[ii])
    names(dot) <- mem$labels
    dot
}

