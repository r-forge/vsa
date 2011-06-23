/*
 * intvec.c: functions for creating and getting/setting a real vector
 */

#include <Rconfig.h> /**/
#include <R.h>
#include <math.h>
#include <Rmath.h>
/* define USE_RINTERNALS so that vector accessors from Rinternals.h are macros */
#define USE_RINTERNALS 1
#include <Rinternals.h>
/* #include <Rdefines.h> */
/* #include <Defn.h> */
/* #include <Rmath.h> */
/* #include <R_ext/RS.h>     for F77_CALL */
#include <R_ext/Applic.h> /* for dgemm */

/*
 * intvec is an external-pointer vector of ints.
 * Three operations:
 *   (1) create
 *   (2) get (sub vector)
 *   (3) set (sub vector)
 * The length of the vector is stored as an integer
 * at the start. The data starts 2 ints into the memory
 * (to conform with how intvec works).
 */
static void intvec_finalizer(SEXP ptr)
{
    int *mem = R_ExternalPtrAddr(ptr);
    if (!mem)
        return;
    Free(mem);
    R_ClearExternalPtr(ptr);
}

SEXP intvec_create(SEXP n) {
    int *p = 0;
    int *i;
    SEXP ptr;
    p = Calloc(INTEGER(n)[0]+2, int);
    ptr = R_MakeExternalPtr(p, install("vecptr"), R_NilValue);
    PROTECT(ptr);
    i = (int *) p;
    *i = INTEGER(n)[0];
    R_RegisterCFinalizerEx(ptr, intvec_finalizer, TRUE);
    UNPROTECT(1);
    return ptr;
}

/*
 * ptr is a vector of ints
 * ptrlen is the length of ptr
 * idx is a vector of integers (1-based)
 * idxlen is the length of idx
 */
SEXP intvec_get(SEXP ptr, SEXP idx) {
    int *mem = ((int*) R_ExternalPtrAddr(ptr)) + 2;
    int ptrlen = *((int*) R_ExternalPtrAddr(ptr));
    int *idxp = INTEGER(idx); /* vecidx is 1-based, k is 0-based */
    int i, k = length(idx);
    SEXP ans = allocVector(INTSXP, k);
    for (i = 0; i < k; i++) {
	if (idxp[i] > ptrlen)
	    error("idx too large");
	if (idxp[i] < 1)
	    error("vecidx must be +ve (is 1-based)");
	INTEGER(ans)[i] = *(mem+idxp[i]-1);
    }
    return ans;
}

/*
 * ptr is a vector of ints
 * ptrlen is the length of ptr
 * idx is a vector of integers (1-based)
 * idxlen is the length of idx
 */
SEXP intvec_getlen(SEXP ptr) {
    int ptrlen = *((int*) R_ExternalPtrAddr(ptr));
    SEXP ans = allocVector(INTSXP, 1);
    INTEGER(ans)[0] = ptrlen;
    return ans;
}

/*
 * ptr is a vector of ints
 * ptrlen is the length of ptr
 * idx is a vector of integers (1-based)
 * idxlen is the length of idx
 * val is a vector of doubles (values to set into ptr)
 */
SEXP intvec_set(SEXP ptr, SEXP idx, SEXP val) {
    int *mem = (int*) R_ExternalPtrAddr(ptr) + 2;
    int ptrlen = *((int*) R_ExternalPtrAddr(ptr));
    int *idxp = INTEGER(idx); /* vecidx is 1-based, k is 0-based */
    int i, k;
    k = length(idx);
    if (length(val) < length(idx))
	error("val is shorter than idx");
    for (i = 0; i < k; i++) {
	if (idxp[i] > ptrlen)
	    error("idx too large");
	if (idxp[i] < 1)
	    error("vecidx must be +ve (is 1-based)");
	*(mem + idxp[i] - 1) = INTEGER(val)[i];
    }
    return(val);
}
