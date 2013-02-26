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
 * realhrr_ramdb provides a RAM-based vector memory that stores data in 1-bit integers with no zeros
 * (i.e., no zero value is used)
 * Implement realhrr_ramdb8: 1-bit approximation for HRR ramdb memory.
 * Bits are stored in int's (either 4 or 8 byte, depending on compiler)
 * Order is from LSB first.
 */

static void realhrr_ramdbx1_finalizer(SEXP ptr)
{
    signed char *mem = (signed char*) R_ExternalPtrAddr(ptr);
    if (!mem)
        return;
    /* close the ramdb and free the storage */
    Free(mem);
    R_ClearExternalPtr(ptr);
}

SEXP realhrr_ramdbx1_create(SEXP n, SEXP m, SEXP tag, SEXP prot) {
    unsigned char *p = 0;
    int intbytes = sizeof(int);
    SEXP ptr;
    /* need to round up nup so that nup/8 bytes is a round number of 'int's and 'n' bits fits in nup/8 bytes */
    int nup = ceil(INTEGER(n)[0] / (intbytes * 8.0)) * intbytes * 8;
    p = Calloc(nup * ceil(INTEGER(m)[0] / 8.0), unsigned char);
    ptr = R_MakeExternalPtr(p, install("ramdbptr"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, realhrr_ramdbx1_finalizer, TRUE);
    UNPROTECT(1);
    return ptr;
}

SEXP realhrr_ramdbx1_get(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx) {
    int intbytes = sizeof(int);
    int intbits = sizeof(int) * 8;
    SEXP ans;
    unsigned int *mem = (unsigned int*) R_ExternalPtrAddr(ptr);
    double *x;
    int n = INTEGER(veclen)[0];
    int nup = ceil(INTEGER(veclen)[0] / (intbytes * 8.0)) * intbits;
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i, j;
    unsigned int z;
    /* want scale st n * scale^2 = 1 */
    double scale = 1.0 / sqrt((double) n);
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    ans = allocVector(REALSXP, n);
    x = REAL(ans);
    mem += k * nup / intbits;
    z = *mem;
    j = 0;
    for (i = 0; i<n; i++) {
        *(x++) = (z & 0x1) ? scale : -scale;
        if (++j >= intbits) {
            j = 0;
            z = *(++mem);
        } else {
            z >>= 1;
        }
    }
    return ans;
}

SEXP realhrr_ramdbx1_getraw(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx) {
    int intbytes = sizeof(int);
    int intbits = sizeof(int) * 8;
    SEXP ans;
    unsigned int *mem = (unsigned int*) R_ExternalPtrAddr(ptr);
    unsigned char *memb;
    unsigned char *x;
    int nup = ceil(INTEGER(veclen)[0] / (intbytes * 8.0)) * intbits;
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i, nbytes;
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    nbytes = nup / 8;
    ans = allocVector(RAWSXP, nbytes);
    x = RAW(ans);
    mem += k * nup / intbits;
    memb = (unsigned char *) mem;
    for (i = 0; i<nbytes; i++)
        *(x++) = *(memb++);
    return ans;
}

SEXP realhrr_ramdbx1_setraw(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx, SEXP vec) {
    SEXP ans;
    int intbytes = sizeof(int);
    int intbits = sizeof(int) * 8;
    unsigned int *mem = (unsigned int*) R_ExternalPtrAddr(ptr);
    unsigned char *memb;
    unsigned char *x;
    int nup = ceil(INTEGER(veclen)[0] / (intbytes * 8.0)) * intbits;
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i, nbytes;
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    nbytes = nup / 8;
    ans = allocVector(INTSXP, 1);
    x = RAW(vec);
    mem += k * nup / intbits;
    memb = (unsigned char *) mem;
    if (length(vec) < nbytes)
        nbytes = length(vec);
    for (i = 0; i<nbytes; i++)
        *(memb++) = *(x++);
    INTEGER(ans)[0] = nbytes;
    return ans;
}

SEXP realhrr_ramdbx1_set(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx, SEXP vec) {
    int intbytes = sizeof(int);
    int intbits = sizeof(int) * 8;
    unsigned int *mem = (unsigned int*) R_ExternalPtrAddr(ptr);
    unsigned int hibit = (0x1 << (intbits - 1));
    double *x;
    int n = INTEGER(veclen)[0];
    int nup = ceil(INTEGER(veclen)[0] / (intbytes * 8.0)) * intbits;
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i, j;
    unsigned int z;
    SEXP ans;
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    x = REAL(vec);
    /* vecidx is 0-based */
    mem += k * nup / intbits;
    z = *mem;
    j = 0;
    z = 0;
    for (i = 0; i<n; i++) {
        if (*(x++) > 0)
            z |= hibit;
        j++;
        if (i == n-1) {
            z >>= (intbits - j);
            *(mem++) = z;
        } else if (j >= intbits) {
            j = 0;
            *(mem++) = z;
            z = 0;
        } else {
            z >>= 1;
        }
    }
    ans = allocVector(REALSXP, 1);
    REAL(ans)[0] = 1.0;
    return ans;
}

/*
 * Compute the dot product of a vector with those stored in the memory.
 * 'memsize' can be supplied as less than the actual size to avoid
 * computation for unused locations in the memory.
 * 'active', if non-null, should be a vector of integers specifying
 * which locations in the memory to use.
 */
SEXP realhrr_ramdbx1_dot(SEXP ptr, SEXP veclen, SEXP memsize, SEXP active, SEXP vec, SEXP param, SEXP mlookptr) {
    int intbytes = sizeof(int);
    int intbits = sizeof(int) * 8;
    unsigned int *mem = (unsigned int *) R_ExternalPtrAddr(ptr);
    unsigned int *y;             /* pointer into mem */
    unsigned int *x, *x1;        /* typecast version of vec */
    unsigned int z;
    unsigned int hibit = (0x1 << (intbits - 1));
    int aa;                      /* for accumulating dotprod */
    double *x0, *a;              /* pointers to vec and ans */
    int n = INTEGER(veclen)[0];
    /* want scale st n * scale^2 = 1 */
    double scale = 1.0 / sqrt((double) n);
    int nup = ceil(INTEGER(veclen)[0] / (intbytes * 8.0)) * intbits;
    int nints = nup / intbits;
    int m = INTEGER(memsize)[0];
    int k, i, j, b, ma, have_active;
    int *mlookup = 0;
    int debug = 0;
    int unroll = 1;
    if (length(param)>=2)
        debug = REAL(param)[1];
    if (length(param)>=3)
        unroll = REAL(param)[2];
    SEXP ans;
    if (!isNull(mlookptr)) {
        mlookup = (int *) R_ExternalPtrAddr(mlookptr);
        if (*mlookup == 256 * 256)
            mlookup += 2;
        else
            mlookup = 0;
    }
    if (isNull(active)) {
        ans = allocVector(REALSXP, m);
        have_active = 0;
    } else {
        ans = allocVector(REALSXP, length(active));
        have_active = 1;
    }
    PROTECT(ans);
    /* convert double vec to unsigned ints in x1 */
    x1 = (unsigned int *) R_alloc(nup / intbits, sizeof(int));
    x = x1;
    x0 = REAL(vec);
    j = 0;
    z = 0;
    for (i = 0; i<n; i++) {
        if (*(x0++) > 0)
            z |= hibit;
        j++;
        if (i == n-1) {
            z >>= (intbits - j);
            *(x++) = z;
        } else if (j >= intbits) {
            j = 0;
            *(x++) = z;
            z = 0;
        } else {
            z >>= 1;
        }
    }

    /* iterate over each vector we want to compare to */
    ma = length(ans);
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
        y = mem + k * nup / intbits;
        x = x1;
        a = REAL(ans)+k;
        aa = 0;
        /* count matching bits */
        if (mlookup == 0) {
            for (i = 0; i < nints; i++) {
                z = (*(x++)) ^ (*(y++));
                for (b = 0; b < intbits; b++) {
                    aa += z & 1;
                    z >>= 1;
                }
            }
        } else {
            /* use lookup table for the xor sum */
            /*
            unsigned char *xu = (unsigned char*) x;
            unsigned char *yu = (unsigned char*) y;
            */
        }
        /* aa was the sum of differences, including padding.
         * Convert to -1/1 sum is aa*-1 + (n-aa)*1 = n - 2aa */
        *a = (n - 2 * aa) * scale * scale;
        if (debug)
            *a = aa;
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
SEXP realhrr_ramdbx1_set_rand(SEXP ptr, SEXP veclen, SEXP memsize, SEXP active, SEXP cnormp, SEXP scalep) {
    int intbytes = sizeof(int);
    int intbits = sizeof(int) * 8;
    unsigned char *x, *mem = (unsigned char*) R_ExternalPtrAddr(ptr);
    SEXP ans;
    double rscale, s2;
    double iscale, riscale, scale;
    float *temp = 0, r, rs;
    int n = INTEGER(veclen)[0];
    int nup = ceil(INTEGER(veclen)[0] / (intbytes * 8.0)) * intbits;
    int m = INTEGER(memsize)[0];
    int cnorm, k, i, j, active_len, have_active;
    i = nup;
    iscale = (127.0 * sqrt((double) n)) / 3.0;
    scale = 3.0 / (127.0 * sqrt((double) n));
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
        rscale = 1.0 / sqrt((double) n);
    } else {
        rscale = REAL(scalep)[0];
    }
    riscale = rscale * iscale;
    if (cnorm)
        temp = (float *) R_alloc(n, sizeof(float));
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
        if (cnorm) {
            x = mem + k * n;
            s2 = 0;
            /* create a float version (unrounded) and accumulate sum of squared-rounded elements */
            for (i = 0; i < n; i++) {
                *temp = norm_rand() * riscale;
                r = round(*temp);
                r = r > 127 ? 127 : (r < -127 ? -127 : r);
                s2 += r * r;
                temp++;
            }
            /* apply the scaling to the float version, then round */
            s2 = sqrt(s2) * scale;
            x = mem + k * n;
            rs = (s2 > 0 ? 1/s2 : 1.0);
            s2 = 0;
            for (i = 0; i < n; i++) {
                r = round(*temp * rs);
                r = r > 127 ? 127 : (r < -127 ? -127 : r);
                *x = r;
                s2 += r * r;
                x++;
                temp++;
            }
            /* not guaranteed to get exactly normalized */
            REAL(ans)[k] = sqrt(s2) * scale;
        } else {
            x = mem + k * n;
            s2 = 0;
            for (i = 0; i < n; i++) {
                r = round(norm_rand() * riscale);
                r = r > 127 ? 127 : (r < -127 ? -127 : r);
                *(x) = r;
                s2 += r * r;
                x++;
            }
            REAL(ans)[k] = sqrt(s2) * scale;
        }
    }
    PutRNGstate();
    return ans;
}
