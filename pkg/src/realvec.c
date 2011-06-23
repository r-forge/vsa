/*
 * realvec.c: functions for creating and getting/setting a real vector
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
 * realvec is an external-pointer vector of floats.
 * Three operations:
 *   (1) create
 *   (2) get (sub vector)
 *   (3) set (sub vector)
 * The length of the vector is stored as an integer
 * at the start. The data starts 2 floats into the memory
 * (to allow for the unlikely possibility that an integer
 * takes twice as many bytes as a float).
 */
static void realvec_finalizer(SEXP ptr)
{
    float *mem = R_ExternalPtrAddr(ptr);
    if (!mem)
        return;
    Free(mem);
    R_ClearExternalPtr(ptr);
}

SEXP realvec_create(SEXP n) {
    float *p = 0;
    int *i;
    SEXP ptr;
    p = Calloc(INTEGER(n)[0]+2, float);
    i = (int *) p;
    *i = INTEGER(n)[0];
    ptr = R_MakeExternalPtr(p, install("vecptr"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, realvec_finalizer, TRUE);
    UNPROTECT(1);
    return ptr;
}

/*
 * ptr is a vector of floats
 * ptrlen is the length of ptr
 * idx is a vector of integers (1-based)
 * idxlen is the length of idx
 */
SEXP realvec_get(SEXP ptr, SEXP idx) {
    float *mem = ((float*) R_ExternalPtrAddr(ptr)) + 2;
    int ptrlen = *((int*) R_ExternalPtrAddr(ptr));
    int *idxp = INTEGER(idx); /* vecidx is 1-based, k is 0-based */
    int i, k = length(idx);
    SEXP ans = allocVector(REALSXP, k);
    for (i = 0; i < k; i++) {
	if (idxp[i] > ptrlen)
	    error("idx too large");
	if (idxp[i] < 1)
	    error("vecidx must be +ve (is 1-based)");
	REAL(ans)[i] = *(mem+idxp[i]-1);
    }
    return ans;
}

/*
 * ptr is a vector of floats
 * ptrlen is the length of ptr
 * idx is a vector of integers (1-based)
 * idxlen is the length of idx
 * val is a vector of doubles (values to set into ptr)
 */
SEXP realvec_set(SEXP ptr, SEXP idx, SEXP val) {
    float *mem = (float*) R_ExternalPtrAddr(ptr) + 2;
    int ptrlen = *((int*) R_ExternalPtrAddr(ptr));
    int *idxp = INTEGER(idx); /* vecidx is 1-based, k is 0-based */
    int i, k = length(idx);
    for (i = 0; i < k; i++) {
	if (idxp[i] > ptrlen)
	    error("idx too large");
	if (idxp[i] < 1)
	    error("vecidx must be +ve (is 1-based)");
	*(mem + idxp[i] - 1) = REAL(val)[i];
    }
    return(val);
}
