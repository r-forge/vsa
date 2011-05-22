#include <Rconfig.h> /**/
#include <R.h>
#include <math.h>
#include <Rmath.h>
/* #include <Rdefines.h> */
/* #include <Defn.h> */
/* #include <Rmath.h> */
/* #include <R_ext/RS.h>     for F77_CALL */
#include <R_ext/Applic.h> /* for dgemm */

/* from http://radfordneal.wordpress.com/2011/05/21/slowing-down-matrix-multiplication-in-r/ */

/* fast version of matrix multiplication skips NA checks -- lets the BLAS do whatever it wants */
/* core R code for matprod() is in src/main/array.c */
void matprod_skipna(double *x, int nrx, int ncx,
	                    double *y, int nry, int ncy, double *z)
{
    char *transN1 = "N", *transN2 = "N", *transT1 = "T";
    int i, int_1 = 1;
    double one = 1.0, zero = 0.0;

    /* Note: ncx will be equal to nry. */

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        if (nrx==1 && ncy==1)
            z[0] = F77_CALL(ddot)(&ncx, x, &int_1, y, &int_1);
        else if (ncy==1)
            F77_CALL(dgemv)(transN1, &nrx, &ncx, &one, x, &nrx, y, &int_1,
                            &zero, z, &int_1);
        else if (nrx==1)
            F77_CALL(dgemv)(transT1, &nry, &ncy, &one, y, &nry, x, &int_1,
                            &zero, z, &int_1);
        else
            F77_CALL(dgemm)(transN1, transN2, &nrx, &ncy, &ncx, &one,
                            x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
        for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}
