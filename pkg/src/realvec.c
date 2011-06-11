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

static void realvec_finalizer(SEXP ptr)
{
    if (!R_ExternalPtrAddr(ptr))
        return;
    /* Free(ptr); */
    R_ClearExternalPtr(ptr);
}

SEXP realvec_create(SEXP n) {
    float *p = 0;
    SEXP ptr;
    p = Calloc(INTEGER(n)[0], float);
    ptr = R_MakeExternalPtr(p, install("vecptr"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, realvec_finalizer, TRUE);
    UNPROTECT(1);
    return ptr;
}

SEXP realvec_get(SEXP ptr, SEXP veclen, SEXP vecidx) {
    SEXP ans;
    float *mem = (float*) R_ExternalPtrAddr(ptr);
    double *x;
    int n = INTEGER(veclen)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    if (k >= n)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    ans = allocVector(REALSXP, 1);
    REAL(ans)[0] = *(mem+k);
    return ans;
}

SEXP realvec_set(SEXP ptr, SEXP veclen, SEXP vecidx, SEXP vec) {
    float *mem = (float*) R_ExternalPtrAddr(ptr);
    double x;
    int n = INTEGER(veclen)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    SEXP ans;
    if (k >= n)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    x = REAL(vec)[0];
    /* vecidx is 0-based */
    REAL(mem)[0] = x;
    ans = allocVector(REALSXP, 1);
    REAL(ans)[0] = x;
    return ans;
}
