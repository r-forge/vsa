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
 * realhrr_ramdb provides a RAM-based vector memory that stores data in floats
 */

static void realhrr_ramdb_finalizer(SEXP ptr)
{
    if (!R_ExternalPtrAddr(ptr))
        return;
    /* close the ramdb */
    /* Free(ptr); */
    R_ClearExternalPtr(ptr);
}

SEXP realhrr_ramdb_create(SEXP n, SEXP m, SEXP tag, SEXP prot) {
    float *p = 0;
    SEXP ptr;
    p = Calloc(INTEGER(n)[0] * INTEGER(m)[0], float);
    ptr = R_MakeExternalPtr(p, install("ramdbptr"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, realhrr_ramdb_finalizer, TRUE);
    UNPROTECT(1);
    return ptr;
}

SEXP realhrr_ramdb_get(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx) {
    SEXP ans;
    float *mem = (float*) R_ExternalPtrAddr(ptr);
    double *x;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i;
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    ans = allocVector(REALSXP, n);
    x = REAL(ans);
    mem += k * n;
    for (i = 0; i<n; i++)
        *(x++) = *(mem++);
    return ans;
}

SEXP realhrr_ramdb_set(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx, SEXP vec) {
    float *mem = (float*) R_ExternalPtrAddr(ptr);
    double *x, s2 = 0;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i;
    SEXP ans;
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    x = REAL(vec);
    /* vecidx is 0-based */
    mem += k * n;
    for (i = 0; i<n; i++) {
        s2 += (*x) * (*x);
        *(mem++) = *(x++);
    }
    ans = allocVector(REALSXP, 1);
    REAL(ans)[0] = sqrt(s2);
    return ans;
}

/*
 * Compute the dot product of a vector with those stored in the memory.
 * 'memsize' can be supplied as less than the actual size to avoid
 * computation for unused locations in the memory.
 * 'active', if non-null, should be a vector of integers specifying
 * which locations in the memory to use.
 */
SEXP realhrr_ramdb_dot(SEXP ptr, SEXP veclen, SEXP memsize, SEXP active, SEXP vec) {
    float *x, *x1, *y, *mem = (float*) R_ExternalPtrAddr(ptr);
    double *x0, *a;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int k, i, j, ma, have_active;
    SEXP ans;
    if (isNull(active)) {
        ans = allocVector(REALSXP, m);
        have_active = 0;
    } else {
        ans = allocVector(REALSXP, length(active));
        have_active = 1;
    }
    PROTECT(ans);
    ma = length(ans);
    x1 = (float *) R_alloc(n, sizeof(float));
    x = x1;
    x0 = REAL(vec);
    for (i = 0; i < n; i++)
        *(x++) = *(x0++);
    for (j = 0; j < ma; j++) {
        if (have_active) {
            k = INTEGER(active)[j] - 1;
            if (k >= m)
                error("element of 'active' too large");
            if (k < 0)
                error("elements of 'active' must be +ve (is 1-based)");
        } else {
            k = j;
        }
        /* k is 0-based */
        y = mem + k * n;
        x = x1;
        a = REAL(ans)+k;
        *a = 0;
        for (i = 0; i < n; i++)
            *a += (*(x++)) * (*(y++));
    }
    UNPROTECT(1);
    return ans;
}

/*
 * Fill memory vectors with appropriate random data.
 * 'memsize' can be supplied as less than the actual size to avoid
 * computation for unused locations in the memory.
 * 'active', if non-null, should be a vector of integers specifying
 * which locations in the memory to use.
 * 'cnormp' is a logical, if TRUE, vectors are normalized to have magnitude exactly 1.
 * 'scalep' is a scale factor, used if cnormp is FALSE and scalep is non-NULL.
 * Return a vector of magnitudes.
 */
SEXP realhrr_ramdb_set_rand(SEXP ptr, SEXP veclen, SEXP memsize, SEXP active, SEXP cnormp, SEXP scalep) {
    float *x, *mem = (float*) R_ExternalPtrAddr(ptr);
    SEXP ans;
    double scale, s2;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int cnorm, k, i, j, active_len, have_active;
    if (isNull(active)) {
        active_len = m;
        have_active = 0;
    } else {
        active_len = length(active);
        have_active = 1;
    }
    ans = allocVector(REALSXP, active_len);
    cnorm = LOGICAL(cnormp)[0];
    if (cnorm || isNull(scalep)) {
        /* keep values at the right scale */
        scale = 1.0 / sqrt((double) n);
    } else {
        scale = REAL(scalep)[0];
    }
    GetRNGstate();
    for (j = 0; j < active_len; j++) {
        if (have_active) {
            k = INTEGER(active)[j] - 1;
            if (k >= m)
                error("element of 'active' too large");
            if (k < 0)
                error("elements of 'active' must be +ve (is 1-based)");
        } else {
            k = j;
        }
        /* k is 0-based */
        x = mem + k * n;
        s2 = 0;
        if (cnorm) {
            for (i = 0; i < n; i++) {
                *x = norm_rand() * scale;
                s2 = (*x)*(*x);
                x++;
            }
            s2 = sqrt(s2);
            if (s2 > 0) {
                x = mem + k * n;
                for (i = 0; i < n; i++) {
                    *x = (*x)/s2;
                    x++;
                }
                REAL(ans)[k] = 1;
            } else {
                REAL(ans)[k] = 0;
            }
        } else {
            for (i = 0; i < n; i++) {
                *(x) = norm_rand() * scale;
                s2 += (*x) * (*x);
                x++;
            }
            REAL(ans)[k] = sqrt(s2);
        }

    }
    PutRNGstate();
    return ans;
}
