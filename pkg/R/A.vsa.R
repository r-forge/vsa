# setOldClass("vsa")
# setOldClass("new.vsa")
# setOldClass("simval")
# setOldClass("vsamem")

.onLoad <- function(libname, pkgname) {
    options(vsalen=1024)
    options(vsanorm=TRUE)
    options(vsatype="realhrr")
}

# The 'vsa' function uses a dummy object as the first argument --
# just the class of the object is important.  This is not intended
# to be called directly -- it gets called from randVec().

newVec <- function(what=c("rand", "I", "1", "0", "NA"),
                    len=options("vsalen")[[1]],
                    elts=NULL,
                    norm=options("vsanorm")[[1]],
		    vsatype=options("vsatype")[[1]]) UseMethod("newVec")

# Create a dummy object with class == options('vsatype')[[1]]
# in order to dispatch the appropriate 'newVec' method

newVec.default <- function(what=c("rand", "I", "1", "0", "NA"),
                    len=options("vsalen")[[1]],
                    elts=NULL,
                    norm=options("vsanorm")[[1]],
		    vsatype=options("vsatype")[[1]]) {
    if (is.null(vsatype))
	stop("no vsatype specified -- must set options(vsatype=...)")
    if (length(vsatype)!=1)
	stop("must set options('vsatype') to the name of a subclass of 'vsa'")
    if (is.null(len))
	stop("no len specified -- must set options(vsalen=...)")
    what <- match.arg(what)
    newVec(structure(what, class=vsatype), len, elts, norm)
}

# Simval
simval <- function(x) {
    if (!is.numeric(x) || length(x)!=1)
	stop("x must be a scalar")
    structure(as.vector(x), class="simval")
}

# Ops.vsa calls vsa vector generics for all operations:
#   vsa * vsa: vsaprod
#   vsa * scalar: vsascale
#   vsa / scalar: vsascale with 1/scalar
#   vsa + vsa: add ('superimpose' the generic can take any number of arguments)
#   vsa - vsa: add with vsascale with -1
#   vsa ^ scalar: vsapower
#   !vsa: appinv

Ops.vsa <-
function (e1, e2 = NULL)
{
    unary <- nargs() == 1
    lclass <- nchar(.Method[1]) > 0
    rclass <- !unary && (nchar(.Method[2]) > 0)
    if (inherits(e1, "simval") || inherits(e2, "simval"))
        stop("'", .Generic, "' cannot operate on a vsa vector and a simval -- does the expression have ()'s around the arguments of '*'? (if the expression is correct, use scalar(simval) in the expression)")
    if (!is.element(.Generic, c("+", "-", "*", "/", "!", "^", "==", "!=")))
        stop("operator '", .Generic, "' is not implemented for vsa vectors")
    if (rclass && lclass) {
        oclass <- class(e1)
        if (oclass[1] != class(e2)[1])
            stop("operator '", .Generic, "' applied to two arguments with different vsa vector classes")
        olength <- length(e1)
        if (olength != length(e2))
            stop("operator '", .Generic, "' applied to two arguments with different vsa vector lengths")
    } else if (rclass) {
        oclass <- class(e2)
    } else if (lclass) {
        oclass <- class(e1)
    } else {
        stop("odd -- neither arugment has a class?")
    }
    if (.Generic=="*") {
        if (lclass && rclass)
            return(vsaprod(e1, e2))
        else if (lclass)
            return(vsascale(e1, e2))
        else if (rclass)
            return(vsascale(e2, e1))
        else
            stop("odd -- neither arugment has a vsa vector class?")
    } else if (.Generic=="^") {
        if (!lclass || rclass)
            stop("'^' when used with vsa vector objects must have the vsa vector on the left and the scalar on the right")
        return(vsapower(e1, e2))
    } else if (.Generic=="!") {
        return(appinv(e1))
    } else if (.Generic=="/") {
        if (rclass)
	    return(vsaprod(e1, appinv(e2)))
	else
	    return(vsascale(e1, 1/e2))
    } else if (.Generic=="+" || .Generic=="-") {
        if (lclass && rclass) {
            if (.Generic=="-")
                e2 <- vsascale(e2, -1)
            return(add(e1, e2))
        } else {
            stop("'", .Generic, "' requires both arguments to be vsa's")
        }
    } else if (.Generic=="==" || .Generic=="!=") {
        if (lclass && rclass) {
            res <- elts(e1) == elts(e2)
            if (.Generic=="==")
                return(all(res))
            else
                return(any(!res))
        } else {
            stop("'", .Generic, "' requires both arguments to be vsa's")
        }
    }
    stop("should not reach here")
}

vsaprod <- function(e1, e2, method=c("fft", "outer")) UseMethod("vsaprod")

vsapower <- function(e1, e2) UseMethod("vsapower")

appinv <- function(e1) UseMethod("appinv")

dot <- function(e1, e2) UseMethod("dot")

"%.%" <- dot

equiv <- function(e1, e2, tol=1e-6) UseMethod("equiv")

"%==%" <- function(e1, e2) equiv(e1, e2)

# can't use function cos() because it is not generic
cosine <- function(e1, e2, mag1=NULL, mag2=NULL) UseMethod("cosine")

"%cos%" <- function(e1, e2) cosine(e1, e2)

norm <- function(e1) UseMethod("norm")

mag <- function(e1, actual=NULL) UseMethod("mag")

scalar <- function(e1) UseMethod("scalar")

scalar.simval <- function(e1) as.vector(e1)

scalar.default <- function(e1) {
    if (!is.numeric(e1) || length(e1)!=1)
	stop("e1 must be a numeric vector of length 1")
    as.vector(e1)
}

add <- function(e1, ...) UseMethod("add")

vsascale <- function(e1, e2) UseMethod("vsascale")

elts <- function(e1) UseMethod("elts")

elts.vsa <- function(e1) as.vector(e1)

print.simval <- function(x, ...) {
    if (is.null(names(x)) && length(x)==1) {
        cat("simval: ", format(x), "\n", sep="")
    } else {
        cat("simval:\n")
        print(unclass(x))
    }
    invisible(x)
}

print.vsa <- function(x, ..., values=FALSE) {
    cat(class(x)[1], "[", length(x), "] mag=",
        format(mag(x)), " cos(.,I)=", format(cosine(newVec(vsatype=class(x)[1], what="I", len=length(x)), x)),
        if (values) ":", "\n", sep="")
    if (values)
        print(unclass(x))
    return(invisible(x))
}

conformable <- function(x=NULL, list, stop.on.error=TRUE) UseMethod("conformable")

conformable.vsa <- function(x=NULL, list, stop.on.error=TRUE)
{
    if (!is.null(x) && !is(x, "vsa"))
        stop("x is not a vsa object")
    if (length(list)==0)
        return(TRUE)
    if (any(!sapply(list, is, "vsa")))
        stop("some object in 'list' are not vsa objects")
    topClass <- function(x) class(x)[1]
    classes <- unique(sapply(list, topClass))
    lengths <- unique(sapply(list, length))
    if (!is.null(x)) {
        classes <- unique(c(classes, topClass(x)))
        lengths <- unique(c(lengths, length(x)))
    }
    if (length(classes)>1)
        if (stop.on.error)
            stop("not all objects have the same class (", paste(classes, collapse=", "), ")")
        else
            return(FALSE)
    if (length(lengths)>1)
        if (stop.on.error)
            stop("not all objects have the same length (", paste(length, collapse=", "), ")")
        else
            return(FALSE)
    return(TRUE)
}

vsamem <- function(..., labels=NULL, type=c("list", "matrix"), call=match.call(expand.dots=FALSE)) {
    type <- match.arg(type)
    if (type=="list")
        return(addmem.default(..., labels=labels, call=call))
    else
        return(addmem.vsamat(..., labels=labels, call=call))
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
                labels <- names(items)
            } else {
                browser()
                if (is.null(labels))
                    labels <- sapply(call$..., function(x) if (is.name(x)) as.character(x) else "")
                else
                    if (length(labels) != length(items))
                        stop("labels are wrong length for items")
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

memlabels <- function(mem) UseMethod("memlabels")

memlabels.vsalist <- function(mem) return(names(mem))
memlabels.vsamat <- function(mem) return(colnames(mem))

memsize <- function(mem) UseMethod("memsize")

memsize.vsalist <- function(mem) return(length(mem))
memsize.vsamat <- function(mem) return(ncol(mem))

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

cosmem <- function(mem, x) UseMethod("cosmem")

# the default should work in most cases
cosmem.default <- function(mem, x) dotmem(mem, x, cos=TRUE)

cleanup <- function(mem, x, cos=TRUE, threshold=NA) UseMethod("cleanup")

# the default should work in most cases
cleanup.default <- function(mem, x, cos=TRUE, threshold=NA) {
    if (length(mem))
        best <- bestmatch(mem, x, cos=cos, n=1, num=TRUE)
    if (length(mem)==0 || (!is.na(threshold) && attr(best, "scores")<threshold))
        res <- newVec(vsatype=class(x)[1], what="NA", len=length(x))
    else
        res <- mem[[best]]
    res
}

bestmatch <- function(mem, x, cos=TRUE, n=1, num=FALSE) UseMethod("bestmatch")

# the default should work in most cases
bestmatch.default <- function(mem, x, cos=TRUE, n=1, num=FALSE) {
    scores <- dotmem(mem, x, cos=cos)
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

# Ops.simval do standard arithmetic on numbers, but check that we're not trying to do something odd to a simval (which can happen if the user doesn't understand the operator precedence & misspecifies and expression -- this helps to catch such situations)

#   simval * simval: not allowed (?)
#   simval * scalar:
#   simval / scalar:
#   simval + simval:
#   simval - simval:
#   simval ^ scalar: power
#   !simval: not allowed (common operator precedence error: !a %.% b is parsed as !(a %.% b), which returns logical, which is almost certainly not what was intended)

Ops.simval <-
function (e1, e2 = NULL)
{
    unary <- nargs() == 1
    lclass <- nchar(.Method[1]) > 0
    rclass <- !unary && (nchar(.Method[2]) > 0)
    if (lclass || rclass)
            stop("'", .Generic, "' cannot operate on a vsa vector and a simval -- does the expression have ()'s around the arguments of '*'? (if the expression is correct, use scalar(simval) in the expression)")
    if (!is.element(.Generic, c("+", "-", "*", "/", "^", "==", "!=")))
        stop("operator '", .Generic, "' is not implemented for simval's (add parenthesis to force intended operator parsing if your expression is otherwise correct)")
    if (rclass && lclass) {
        oclass <- class(e1)
        if (oclass[1] != class(e2)[1])
            stop("operator '", .Generic, "' applied to two arguments with different simval classes")
    } else if (rclass) {
        oclass <- class(e2)
    } else if (lclass) {
        oclass <- class(e1)
    } else {
        stop("odd -- neither arugment has a class?")
    }
    if (!rclass && !is.numeric(e2))
	stop("can only combine simval's with numbers")
    if (!lclass && !is.numeric(e1))
	stop("can only combine simval's with numbers")
    if (.Generic=="*") {
	if (rclass && lclass)
	    stop("'*' must have a simval as one argument and a number as the other")
	return(simval(unclass(e1) * unclass(e2)))
    } else if (.Generic=="^") {
        if (!inherits(e1, "simval") || inherits(e2, "simval"))
            stop("'^' when used with simval's must have the simval on the left and the scalar on the right")
	return(simval(unclass(e1) ^ e2))
    } else if (.Generic=="/") {
        if (rclass)
            stop("cannot have a simval object on the rhs of '/'")
        return(simval(unclass(e1)/e2))
    } else if (.Generic=="+" || .Generic=="-") {
        if (inherits(e1, "simval") && inherits(e2, "simval")) {
            if (.Generic=="-")
	        return(simval(unclass(e1) - unclass(e2)))
 	    else
	        return(simval(unclass(e1) + unclass(e2)))
        } else {
            stop("'", .Generic, "' requires both arguments to be simval's")
        }
    } else if (.Generic=="==" || .Generic=="!=") {
	if (rclass)
	    e2 <- unclass(e2)
	if (lclass)
	    e1 <- unclass(e1)
        res <- (e1 == e2)
        if (.Generic=="==")
            return(all(res))
        else
            return(any(!res))
    }
    stop("should not reach here")
}

addnorm <- function(...) norm(add(...))
