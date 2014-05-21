/* Internal header file for cminpack, by Frederic Devernay. */
#ifndef __CMINPACKP_H__
#define __CMINPACKP_H__

#ifndef __CMINPACK_H__
#error "cminpackP.h in an internal cminpack header, and must be included after all other headers (including cminpack.h)"
#endif

#if (defined (USE_CBLAS) || defined (USE_LAPACK)) && !defined (__cminpack_double__)
#error "cminpack can use cblas and lapack only in double precision mode"
#endif

#ifdef USE_CBLAS
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
#define __cminpack_enorm__(n,x) cblas_dnrm2(n,x,1)
#else
#define __cminpack_enorm__(n,x) __cminpack_func__(enorm)(n,x)
#endif

#ifdef USE_LAPACK
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int 		__CLPK_integer;
typedef int 		__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef int 		__CLPK_ftnlen;
#else
typedef long int 	__CLPK_integer;
typedef long int 	__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef long int 	__CLPK_ftnlen;
#endif
//extern void dlartg_(double *f, double *g, double *cs, double *sn, double *r__);
int dlartg_(__CLPK_doublereal *f, __CLPK_doublereal *g, __CLPK_doublereal *cs,
            __CLPK_doublereal *sn, __CLPK_doublereal *r__);
//extern void dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
int dgeqp3_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
            lda, __CLPK_integer *jpvt, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork,
            __CLPK_integer *info);
//extern void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
int dgeqrf_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_doublereal *a, __CLPK_integer *
            lda, __CLPK_doublereal *tau, __CLPK_doublereal *work, __CLPK_integer *lwork, __CLPK_integer *info);
#endif
#endif

#define real __cminpack_real__
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE_ (1)
#define FALSE_ (0)

#endif /* !__CMINPACKP_H__ */
