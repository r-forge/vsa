create.ramdb.realhrr <- function(memsize, example=newVec(), fill=FALSE, labels=NULL) {
    if (!inherits(example, "vsa"))
        stop("'example' must be a vsa object")
    if (!is.null(labels) && (length(labels)!=memsize || !is.character(labels)))
        stop("labels must be character vector with length = memsize")
    mem <- .Call("realhrr_ramdb_create", length(example), as.integer(memsize), NULL, NULL)
    attr(mem, "veclen") <- length(example)
    attr(mem, "memsize") <- as.integer(memsize)
    cnorm <- as.logical(getOption("vsacnorm", TRUE))
    if (fill)
        attr(mem, "mag") <- .Call("realhrr_ramdb_setrand", mem, length(example), as.integer(memsize), NULL, cnorm, NULL)
    else
        attr(mem, "mag") <- numeric(memsize)
    attr(mem, "example") <- example
    if (!is.null(labels))
        attr(mem, "labels") <- labels
    class(mem) <- c("ramdb.realhrr", "vsamem")
    mem
}

memlabels.ramdb.realhrr <- function(mem) return(attr(mem, "labels"))
memsize.ramdb.realhrr <- function(mem) return(attr(mem, "memsize"))

print.ramdb.realhrr <- function(x, ...) {
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

***** working on getmem and setmem


getmem.ramdb.realhrr <- function(mem, i) {
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
    example[] <- mem[,j]
    return(example)
}

setmem.ramdb.realhrr <- function(mem, i, x, label=NULL) {
    example <- attr(mem, "example")
    conformable(x, example)
    if (is.character(i) && length(i)==1) {
        j <- match(i, memlabels(mem))
        if (is.na(j)) {
            j <- which(is.na(memlabels(mem)))[1]
            if (is.na(j))
                stop("label '", i, "' not found, and no empty labels to use")
            memlabels(mem)[j] <- i
        }
    } else if (is.numeric(i) && length(i)==1) {
        j <- i
        if (i<1 || i>memsize(mem))
            stop("index ", i, " out of range (mem has ", memsize(mem), " elements)")
        if (!is.null(label))
            memlabels(mem)[j] <- label
    } else {
        stop("i must be single character or numeric")
    }
    mem[,j] <- x[]
    return(j)
}
