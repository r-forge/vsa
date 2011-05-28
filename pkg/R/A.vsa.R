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
