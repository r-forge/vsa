matprodsna <- function(x, y) {
    if (storage.mode(x)!="double")
        stop("storage.mode(x)!='double'")
    if (storage.mode(y)!="double")
        stop("storage.mode(y)!='double'")
    if (!is.matrix(x))
        stop("x is not a matrix")
    if (!is.matrix(y))
        stop("y is not a matrix")
    if (ncol(x) != nrow(y))
        stop("ncol(x) != nrow(y)")
    z <- matrix(numeric(nrow(x) * ncol(y)), nrow=nrow(x), ncol=ncol(y), dimnames=list(rownames(x), colnames(y)))
    .C("matprod_skipna", x, nrow(x), ncol(x), y, nrow(y), ncol(y), z, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    z
}

# X <- matrix(c(1,2,3,4), ncol=2)
# Y <- matrix(c(3,1,4,2), ncol=2)
# X %*% Y
#     [,1] [,2]
# [1,]    6   10
# [2,]   10   16
# matprodsna(X, Y)
#     [,1] [,2]
# [1,]    6   10
# [2,]   10   16

crossprodsna <- function(x, y) {
    if (storage.mode(x)!="double")
        stop("storage.mode(x)!='double'")
    if (storage.mode(y)!="double")
        stop("storage.mode(y)!='double'")
    if (!is.matrix(x))
        stop("x is not a matrix")
    if (!is.matrix(y))
        stop("y is not a matrix")
    if (nrow(x) != nrow(y))
        stop("nrow(x) != nrow(y)")
    z <- matrix(numeric(ncol(x) * ncol(y)), nrow=ncol(x), ncol=ncol(y), dimnames=list(colnames(x), colnames(y)))
    .C("crossprod_skipna", x, nrow(x), ncol(x), y, nrow(y), ncol(y), z, DUP=FALSE, NAOK=TRUE, PACKAGE="vsa")
    z
}

# crossprod(X, Y)
#      [,1] [,2]
# [1,]    5    8
# [2,]   13   20
# crossprodsna(X, Y)
