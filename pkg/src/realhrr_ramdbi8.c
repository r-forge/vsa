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
 * realhrr_ramdb provides a RAM-based vector memory that stores data in 8-bit integers with zeros
 * (i.e., a zero value is used, -128 is not used)
 * implement realhrr_ramdb8: 8bit approximation for HRR ramdb memory
 * all ordinary dreps stored as regular doubles
 * 8-bit fractional used for ramdb
 * use 255 signed numbers with range +/-3 sd.
 * sd = sqrt(1/n)
 * scale * 127 = 3*sqrt(1/n)
 * => scale = 3/127 * sqrt(1/n)
 * xm is 8bit signed int, x is double version
 * x = 3/127 * xm * sqrt(1/n)
 * realhrr_ramdbi8_get/_set interface with floats on 3/127 * sqrt(1/n) scale
 */

static void realhrr_ramdbi8_finalizer(SEXP ptr)
{
    signed char *mem = (signed char*) R_ExternalPtrAddr(ptr);
    if (!mem)
        return;
    /* close the ramdb and free the storage */
    Free(mem);
    R_ClearExternalPtr(ptr);
}

SEXP realhrr_ramdbi8_create(SEXP n, SEXP m, SEXP tag, SEXP prot) {
    signed char *p = 0;
    SEXP ptr;
    p = Calloc(INTEGER(n)[0] * INTEGER(m)[0], signed char);
    ptr = R_MakeExternalPtr(p, install("ramdbptr"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, realhrr_ramdbi8_finalizer, TRUE);
    UNPROTECT(1);
    return ptr;
}

SEXP realhrr_ramdbi8_get(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx) {
    SEXP ans;
    signed char *mem = (signed char*) R_ExternalPtrAddr(ptr);
    double *x;
    float scale;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i;
    scale = 3.0 / (127.0 * sqrt((double) n));
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    ans = allocVector(REALSXP, n);
    x = REAL(ans);
    mem += k * n;
    for (i = 0; i<n; i++)
        *(x++) = *(mem++) * scale;
    return ans;
}

SEXP realhrr_ramdbi8_set(SEXP ptr, SEXP veclen, SEXP memsize, SEXP vecidx, SEXP vec) {
    signed char *mem = (signed char*) R_ExternalPtrAddr(ptr);
    double *x, s2 = 0;
    float xi, iscale;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int k = INTEGER(vecidx)[0] - 1; /* vecidx is 1-based, k is 0-based */
    int i;
    SEXP ans;
    iscale = (127.0 * sqrt((double) n)) / 3.0;
    if (k >= m)
        error("vecidx too large");
    if (k < 0)
        error("vecidx must be +ve (is 1-based)");
    x = REAL(vec);
    /* vecidx is 0-based */
    mem += k * n;
    for (i = 0; i<n; i++) {
	xi = (*x) * iscale;
	xi = xi > 127 ? 127 : (xi < -127 ? -127 : xi);
	xi = round(xi);
        s2 += xi * xi;
        *(mem++) = (signed char) xi;
	x++;
    }
    ans = allocVector(REALSXP, 1);
    REAL(ans)[0] = sqrt(s2) / iscale;
    return ans;
}

/*
 * Compute the dot product of a vector with those stored in the memory.
 * 'memsize' can be supplied as less than the actual size to avoid
 * computation for unused locations in the memory.
 * 'active', if non-null, should be a vector of integers specifying
 * which locations in the memory to use.
 */
SEXP realhrr_ramdbi8_dot(SEXP ptr, SEXP veclen, SEXP memsize, SEXP active, SEXP vec, SEXP param, SEXP mlookptr) {
    signed char *mem = (signed char*) R_ExternalPtrAddr(ptr);
    signed char *y;		/* pointer into mem */
    signed char *x, *x1;	/* typecast version of vec */
    int aa;			/* for accumulating dotprod */
    double *x0, *a;		/* pointers to vec and ans */
    float scale, iscale;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int k, i, j, ma, have_active;
    int *mlookup = 0;
    int debug = 0;
    int unroll = 1;
    if (length(param)>=2)
	debug = REAL(param)[1];
    if (length(param)>=3)
	unroll = REAL(param)[2];
    SEXP ans;
    scale = 3.0 / (127.0 * sqrt((double) n));
    iscale = (127.0 * sqrt((double) n)) / 3.0;
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
    /* convert double vec to signed chars in x1 */
    x1 = (signed char *) R_alloc(n, sizeof(signed char));
    x = x1;
    x0 = REAL(vec);
    for (i = 0; i < n; i++) {
        *x = round(iscale * (*(x0++)));
	if (*x > 127)
	    *x = 127;
	else if (*x < -127)
	    *x = -127;
	x++;
    }
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
        y = mem + k * n;
        x = x1;
        a = REAL(ans)+k;
	aa = 0;
	if (mlookup == 0) {
	    if (n % 4 == 0 && unroll) {
		/* unroll loops in steps of 4 */
		for (i = 0; i < n; i += 4) {
		    aa += x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
		    x += 4;
		    y += 4;
		}
	    } else if (1) {
		for (i = 0; i < n; i++)
		    aa += (*(x++)) * (*(y++));
	    } else {
		/* unroll loops in steps of 8 */
		for (i = 0; i < n; i += 8) {
		    aa += (  x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
			     + x[4] * y[4] + x[5] * y[5] + x[6] * y[6] + x[7] * y[7]);
		    x += 8;
		    y += 8;
		}
	    }
	} else {
	    /* use lookup table for signed character multiplication */
	    unsigned char *xu = (unsigned char*) x;
	    unsigned char *yu = (unsigned char*) y;
	    if (debug>0) {
		if (debug==1) {
		    /* lookup table value for first element */
		    aa = mlookup[(*(xu)) << 8 | (*(yu))];
		} else if (debug==2) {
		    /* offset into lookup table value for first element */
		    aa = (*(xu)) << 8 | (*(yu));
		} else {
		    aa = 77;
		}
	    } else if (n % 4 == 0 && unroll) {
		/* unroll loops in steps of 4 */
		/* need to use xu and yu instead of x and y */
		for (i = 0; i < n; i += 4) {
		    aa += (mlookup[xu[0] << 8 | yu[0]] + mlookup[xu[1] << 8 | yu[1]]
			   + mlookup[xu[2] << 8 | yu[2]] + mlookup[xu[3] << 8 | yu[3]]);
		    xu += 4;
		    yu += 4;
		}
	    } else {
		for (i = 0; i < n; i++)
		    aa += mlookup[(*(xu++)) << 8 | (*(yu++))];
	    }
	}
        *a = aa * scale * scale;
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
SEXP realhrr_ramdbi8_set_rand(SEXP ptr, SEXP veclen, SEXP memsize, SEXP active, SEXP cnormp, SEXP scalep) {
    signed char *x, *mem = (signed char*) R_ExternalPtrAddr(ptr);
    SEXP ans;
    double rscale, s2;
    float iscale, riscale, scale;
    float *temp = 0, r, rs;
    int n = INTEGER(veclen)[0];
    int m = INTEGER(memsize)[0];
    int cnorm, k, i, j, active_len, have_active;
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
