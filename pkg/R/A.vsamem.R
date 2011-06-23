vsamem <- function(..., labels=NULL, type=c("list", "matrix", "db"), call=match.call(expand.dots=FALSE)) {
    type <- match.arg(type)
    if (type=="list") {
        return(addmem.default(..., labels=labels, call=call))
    } else if (type=="matrix") {
        return(addmem.vsamat(..., labels=labels, call=call))
    } else if (type=="db") {
        # addmem.vsadb has additional named arguments, but need compute call at this level to get names
        # addmem.vsadb will have to remove other matched formal arguments from call$...
        if (getOption("vsatype")=="realhrr")
            return(addmem.vsadb(..., labels=labels, call=call))
        else
            stop("db memory not implemented for options('vsatype')=", getOption('vsatype'))
    }
}

addmem <- function(..., labels=NULL, call=match.call(expand.dots=FALSE)) UseMethod("addmem")

addmem.vsamat <- function(..., labels=NULL, call=match.call(expand.dots=FALSE)) {
    items <- list(...)
    # Four modes of operation:
    # (1) convert a matrix or named list of vsa's (first arg) to a vsamat
    # (2) concatenate several vsamat's together,
    # (3) add several vsa's to a vsalist,
    # (4) create a new vsalist from vsa arguments
    # In the cases (3) and (4), take labels from actual argument names.
    if ((is.list(items[[1]]) || is.matrix(items[[1]])) && !inherits(items[[1]], "vsamat")) {
        if (length(items) > 1)
            stop("can only supply one arg when converting an ordinary matrix or list to a vsamat")
        items <- items[[1]]
        if (is.null(labels))
            labels <- names(items)
        if (is.list(items)) {
            if (length(items)>=2)
                conformable(items[[1]], items[-1])
            example <- if (length(items)>0) items[[1]]
            items <- matrix(unlist(items), ncol=length(items), dimnames=list(NULL, labels))
            attr(items, "example") <- example
        }
    } else {
        if (inherits(items[[1]], "vsa")) {
            # args must be all vsa's
            if (length(items)>=2)
                conformable(items[[1]], items[-1])
            if (!is.null(names(items))) {
                # This is the case where args are named
                labels <- names(items)
            } else {
                if (is.null(labels)) {
                    # Take labels from names of actual arguments
                    labels <- sapply(call$..., function(x) if (is.name(x)) as.character(x) else "")
                } else {
                    if (length(labels) != length(items))
                        stop("labels are wrong length for items")
                }
            }
            example <- if (length(items)>0) items[[1]]
            items <- matrix(unlist(items), ncol=length(items), dimnames=list(NULL, labels))
            attr(items, "example") <- example
        } else if (inherits(items[[1]], "vsamat")) {
            if (length(items)==1) {
                # nothing to add ...
                items <- items[[1]]
            } else if (inherits(items[[2]], "vsa")) {
                # remaining args must be all vsa's, this will be checked in the call to conformable() below
                if (is.null(names(items))) {
                    if (is.null(labels))
                        labels <- sapply(call$..., function(x) if (is.name(x)) as.character(x) else "")[-1]
                    else
                        if (length(labels) != (length(items)-1))
                            stop("labels are wrong length for items")
                    names(items)[seq(2, len=length(items)-1)] <- labels
                }
                example <- attr(items[[1]], "example")
                if (is.null(example)) example <- items[[2]]
                conformable(example, items[-1])
                items <- do.call("cbind", items)
                attr(items, "example") <- example
            } else if (inherits(items[[2]], "vsamat")) {
                # remaining args must be all vsamat's
                if (any(!sapply(items, is, "vsamat")))
                    stop("first two args are vsamat objects, but some others are not")
                example <- NULL
                for (i in seq(along=items)) {
                    i.example <- attr(items[[i]], "example")
                    if (!is.null(example) && !is.null(i.example))
                        conformable(example, list(i.example))
                    if (is.null(example))
                        example <- i.example
                }
                items <- do.call("cbind", items)
                attr(items, "example") <- example
            } else {
                stop("when first arg is vsalist, remaining ones must be all vsa or all vsalist")
            }
        } else {
            stop("first arg must be vsa, vsamat or a list or matrix")
        }
    }
    if (ncol(items)>0) {
        if (is.null(colnames(items)) || any(colnames(items)=="") || any(duplicated(colnames(items))))
            stop("all items to store must have unique non-empty names")
    }
    class(items) <- c("vsamat", "vsamem")
    if (ncol(items)>0 && !is.null(example <- attr(items, "example")))
        attr(items, "mag") <- apply(unclass(items), 2, function(x) mag(example, actual=x))
    return(items)
}

addmem.default <- function(..., labels=NULL, call=match.call(expand.dots=FALSE)) {
    items <- list(...)
    # Four modes of operation:
    # (1) convert a named list of vsa's (first arg) to a vsamem (i.e., just add the class)
    # (2) concatenate several vsalist's together,
    # (3) add several vsa's to a vsalist,
    # (4) create a new vsalist from vsas arguments
    # In the cases (3) and (4), take labels from actual argument names.
    if (is.list(items[[1]]) & !inherits(items[[1]], "vsalist")) {
        if (length(items) > 1)
            stop("can only supply one arg when converting an ordinary list to a vsalist")
        items <- items[[1]]
        if (is.null(names(items)))
            names(items) <- labels
    } else {
        if (inherits(items[[1]], "vsa")) {
            # args must be all vsa's, this will be checked in the call to conformable() below
            if (is.null(names(items))) {
                if (is.null(labels))
                    labels <- sapply(call$..., function(x) if (is.name(x)) as.character(x) else "")
                else
                    if (length(labels) != length(items))
                        stop("labels are wrong length for items")
                names(items) <- labels
            }
        } else if (inherits(items[[1]], "vsalist")) {
            if (length(items)==1) {
                # nothing to add ...
                items <- items[[1]]
            } else if (inherits(items[[2]], "vsa")) {
                # remaining args must be all vsa's, this will be checked in the call to conformable() below
                if (is.null(names(items))) {
                    if (is.null(labels))
                        labels <- sapply(call$..., function(x) if (is.name(x)) as.character(x) else "")[-1]
                    else
                        if (length(labels) != (length(items)-1))
                            stop("labels are wrong length for items")
                    names(items)[seq(2, len=length(items)-1)] <- labels
                }
                items <- c(items[[1]], items[-1])
            } else if (inherits(items[[2]], "vsalist")) {
                # remaining args must be all vsalist's
                if (any(!sapply(items, is, "vsalist")))
                    stop("first two args are vsalist objects, but some others are not")
                items <- unlist(items, recursive=FALSE)
            } else {
                stop("when first arg is vsalist, remaining ones must be all vsa or all vsalist")
            }
        } else {
            stop("first arg must be vsa, vsalist or a list")
        }
    }
    if (length(items)>0) {
        if (is.null(names(items)) || any(names(items)=="") || any(duplicated(names(items))))
            stop("all items to store must have unique non-empty names")
        if (length(items) > 1)
            conformable(items[[1]], items[-1])
    }
    class(items) <- c("vsalist", "vsamem")
    if (length(items)>0)
        attr(items, "mag") <- lapply(unclass(items), mag)
    return(items)
}

getmem <- function(mem, i) UseMethod("getmem")

getmem.vsamat <- function(mem, i) {
    example <- attr(mem, "example")
    if (is.character(i) && length(i)==1) {
        j <- match(i, colnames(mem))
        if (is.na(j))
            stop("label '", i, "' not found")
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>ncol(mem))
            stop("index ", i, " out of range (mem has ", ncol(mem), " elements)")
    } else {
        stop("i must be single character or numeric")
    }
    example[] <- mem[,j]
    return(example)
}

getmem.vsalist <- function(mem, i) {
    if (is.character(i) && length(i)==1) {
        j <- match(i, colnames(mem))
        if (is.na(j))
            stop("label '", i, "' not found")
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>length(mem))
            stop("index ", i, " out of range (mem has ", length(mem), " elements)")
    } else {
        stop("i must be single character or numeric")
    }
    return(mem[[j]])
}

setmem <- function(mem, i, x, label=NULL) UseMethod("setmem")

setmem.vsamat <- function(mem, i, x, label=NULL) {
    example <- attr(mem, "example")
    conformable(x, example)
    if (is.character(i) && length(i)==1) {
        j <- match(i, colnames(mem))
        if (is.na(j)) {
            j <- which(is.na(colnames(mem)))[1]
            if (is.na(j))
                stop("label '", i, "' not found, and no empty labels to use")
            colnames(mem)[j] <- i
        }
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>ncol(mem))
            stop("index ", i, " out of range (mem has ", ncol(mem), " elements)")
        if (!is.null(label))
            colnames(mem)[j] <- label
    } else {
        stop("i must be single character or numeric")
    }
    mem[,j] <- x[]
    return(j)
}

setmem.vsalist <- function(mem, i, x, label=NULL) {
    example <- attr(mem, "example")
    conformable(x, example)
    if (is.character(i) && length(i)==1) {
        j <- match(i, names(mem))
        if (is.na(j))
            stop("label '", i, "' not found, and no empty labels to use")
        names(mem)[j] <- i
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>length(mem))
            stop("index ", i, " out of range (mem has ", length(mem), " elements)")
        if (!is.null(label))
            names(mem)[j] <- label
    } else {
        stop("i must be single character or numeric")
    }
    mem[[j]][] <- x[]
    return(j)
}

memlabels <- function(mem) UseMethod("memlabels")

memlabels.vsalist <- function(mem) return(names(mem))
memlabels.vsamat <- function(mem) return(colnames(mem))

memsize <- function(mem) UseMethod("memsize")

memsize.vsalist <- function(mem) return(length(mem))
memsize.vsamat <- function(mem) return(ncol(mem))

memmags <- function(mem) UseMethod("memmags")
memmags.default <- function(mem) attr(mem, "mag")

delmem <- function(mem, items) {
    # 'items' is a vector of the names of items to delete from the memory
    stop("not yet implemented")
}

dotmem <- function(mem, x, ..., cos=FALSE) UseMethod("dotmem")

dotmem.vsalist <- function(mem, x, ..., cos=FALSE) {
    if (!inherits(x, "vsa"))
        stop("x must be a vsa")
    if (!inherits(mem, "vsalist"))
        stop("mem must be a vsalist")
    xmag <- mag(x)
    if (cos)
        res <- sapply(unclass(mem), cosine, e1=x, mag1=xmag)
    else
        res <- sapply(unclass(mem), dot, e1=x)
    names(res) <- memlabels(mem)
    res
}

dotmem.vsamat <- function(mem, x, ..., cos=FALSE) {
    if (!inherits(x, "vsa"))
        stop("x must be a vsa")
    if (!inherits(mem, "vsamat"))
        stop("mem must be a vsamat")
    example <- attr(mem, "example")
    if (!is.null(example))
        conformable(example, list(x))
    if (memsize(mem)==0)
        return(numeric(0))
    if (is.null(example))
        stop("non-empty vsamat memory contains no example vector")
    return(dotmem.vsamat.compute(x, mem, ..., cos=cos))
}

dotmem.vsamat.compute <- function(x, mem, ..., cos) UseMethod("dotmem.vsamat.compute")

# The default will work, but a more efficent version can be supplied that
# dispatches off the vsa subclass (i.e., the type of the vsa vector).
# The method can safely assume the columns of mem conform with x.
dotmem.vsamat.compute.default <- function(x, mem, cos, ...) {
    res <- numeric(memsize(mem))
    names(res) <- colnames(mem)
    xmag <- mag(x)
    y <- x
    for (i in seq(ncol(mem))) {
        y[] <- mem[,i]
        if (cos)
            res[i] <- cosine(x, y, mag1=xmag)
        else
            res[i] <- dot(x, y)
    }
    res
}

cosmem <- function(mem, x, ...) UseMethod("cosmem")

# the default should work in most cases
cosmem.default <- function(mem, x, ...) dotmem(mem, x, cos=TRUE)

cleanup <- function(mem, x, ..., cos=TRUE, threshold=NA) UseMethod("cleanup")

# the default should work in most cases
cleanup.default <- function(mem, x, ..., cos=TRUE, threshold=NA) {
    if (length(mem))
        best <- bestmatch(mem, x, ..., cos=cos, n=1, num=TRUE)
    if (length(mem)==0 || (!is.na(threshold) && attr(best, "scores")<threshold))
        res <- newVec(vsatype=class(x)[1], what="NA", len=length(x))
    else
        res <- mem[[best]]
    res
}

bestmatch <- function(mem, x, ..., cos=TRUE, n=1, num=FALSE) UseMethod("bestmatch")

# the default should work in most cases
bestmatch.default <- function(mem, x, ..., cos=TRUE, n=1, num=FALSE) {
    scores <- dotmem(mem, x, ..., cos=cos)
    best <- order(-scores)[seq(len=min(n, length(mem)))]
    scores <- scores[best]
    if (!num) {
        return(scores)
    } else {
        return(structure(best, scores=scores))
    }
}

print.vsalist <- function(x, ...) {
    vsaclass <- if (length(x)) class(x[[1]]) else "vsa"
    len <- if (length(x)) length(x[[1]]) else "?"
    m <- memlabels(x)
    n <- memsize(x)
    cat("List memory containing ", n, " ", vsaclass[1], "[", len, "]: ",
        paste(m[seq(1, len=min(length(m), 3))], collapse=", "),
        if (length(m)>3) " ...", "\n", sep="")
    invisible(x)
}

print.vsamat <- function(x, ...) {
    vsaclass <- if (!is.null(attr(x, "example"))) class(attr(x, "example")) else "vsa"
    len <- nrow(x)
    m <- memlabels(x)
    n <- memsize(x)
    cat("Matrix memory containing ", n, " ", vsaclass[1], "[", len, "]: ",
        paste(m[seq(1, len=min(length(m), 3))], collapse=", "),
        if (length(m)>3) " ...", "\n", sep="")
    invisible(x)
}

contents <- function(mem, ...) UseMethod("contents")

contents.vsamat <- function(mem, ...) {
    x <- as.matrix(mem)
    class(x) <- NULL
    x
}
