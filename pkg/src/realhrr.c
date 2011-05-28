#include <Rconfig.h> /**/
#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
/* #include <Rdefines.h> */
/* #include <Defn.h> */
/* #include <Rmath.h> */
/* #include <R_ext/RS.h>     for F77_CALL */
#include <R_ext/Applic.h> /* for dgemm */

/*
 * x[n] is a vector
 * y[n,m] is a matrix (m columns)
 * z[m] is a vector
 */
void realhrr_dotmem(double *x, int *n, double *y, int *m, double *z)
{
    int i, j;
    for (i = 0; i < *m; i++) {
   	z[i] = 0;
	for (j = 0; j < *n; j++)
	    z[i] += x[j] * y[(i*(*n))+j];
    }
}

/*
 * x[n] is a vector
 * y[n,m] is a matrix (m columns)
 * z[m] is a vector
 */
void realhrr_dotmempp(double *x, int *n, double *y, int *m, double *z)
{
    int i, j;
    double *xx, *yy, *zz;
    for (i = 0; i < *m; i++) {
	xx = x;
	yy = y + (i*(*n));
	zz = z + i;
   	*zz = 0;
	for (j = 0; j < *n; j++)
	    *zz += *(xx++) * (*(yy++));
    }
}

/*
 * x[n] is a vector
 * y[n,m] is a matrix (m columns)
 * z[m] is a vector
 */
void realhrr_dotmempp1(double *x, int *n, double *y, int *m, double *z)
{
    int i;
    double *xx, *yy, *zz, *xx1;
    xx1 = x + (*n);
    for (i = 0; i < *m; i++) {
	xx = x;
	yy = y + (i*(*n));
	zz = z + i;
   	*zz = 0;
	while (xx < xx1)
	    *zz += *(xx++) * (*(yy++));
    }
}

#ifdef SUPPORT_OPENMP
/*
 * x[n] is a vector
 * y[n,m] is a matrix (m columns)
 * z[m] is a vector
 */
void realhrr_dotmemmp1(double *x, int *n, double *y, int *m, double *z, int *coreused)
{
    int i, tid;
    double *xx, *yy, *zz, *xx1;
    xx1 = x + (*n);
#pragma omp parallel for
    for (i = 0; i < *m; i++) {
	xx = x;
	yy = y + (i*(*n));
	zz = z + i;
	*zz = 0;
	while (xx < xx1)
		*zz += *(xx++) * (*(yy++));
    }
}

/*
 * x[n] is a vector
 * y[n,m] is a matrix (m columns)
 * z[m] is a vector
 */
void realhrr_dotmemmp2(double *x, int *n, double *y, int *m, double *z, int *coreused)
{
    int i, tid, nthreads;
    double *xx, *yy, *zz, *xx1;
    xx1 = x + (*n);
#pragma omp parallel private(tid, nthreads)
    {
	tid = omp_get_thread_num();
	nthreads = omp_get_num_threads();
	for (i = tid; i < *m; i += nthreads) {
	    xx = x;
	    yy = y + (i*(*n));
	    zz = z + i;
	    *zz = 0;
	    while (xx < xx1)
		*zz += *(xx++) * (*(yy++));
	}
	if (tid < 4)
	    coreused[tid]++;
	if (tid==0)
	    coreused[4] = nthreads;
    }
}
#endif

SEXP realhrr_create_vecdb(SEXP n, SEXP m, SEXP tag, SEXP prot) {
    int *p = 0;
    SEXP ptr;
    p = Calloc(INTEGER(n)[0] * INTEGER(m)[0], double);
    ptr = R_MakeExternalPtr(p, install("vecdbptr"), R_NilValue);
    /* PROTECT(ptr); */
    return ptr;
}

static void vecdb_finalizer(SEXP ptr)
{
    if (!R_ExternalPtrAddr(ptr))
	return;
    /* close the vec db */
    R_ClearExternalPtr(ptr);
}


