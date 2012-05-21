#include "../gmx_lapack.h"


/* Normally, SSTEVR is the LAPACK wrapper which calls one
 * of the eigenvalue methods. However, our code includes a
 * version of SSTEGR which is never than LAPACK 3.0 and can
 * handle requests for a subset of eigenvalues/vectors too,
 * and it should not need to call SSTEIN.
 * Just in case somebody has a faster version in their lapack
 * library we still call the driver routine, but in our own
 * case this is just a wrapper to sstegr.
 */
void
F77_FUNC(sstevr,SSTEVR)(const char *jobz, 
                        const char *range,
                        int *n,
                        float *d,
                        float *e,
                        float *vl, 
                        float *vu,
                        int *il, 
                        int *iu, 
                        float *abstol,
                        int *m,
                        float *w, 
                        float *z,
                        int *ldz,
                        int *isuppz, 
                        float *work, 
                        int *lwork, 
                        int *iwork,
                        int *liwork, 
                        int *info)
{
  F77_FUNC(sstegr,SSTEGR)(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w,
	  z, ldz, isuppz, work, lwork, iwork, liwork, info);
  

    return;

}


