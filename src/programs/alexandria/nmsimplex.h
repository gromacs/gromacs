#ifndef _nmsimplex_h
#define _nmsimplex_h
	
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*nm_target_func)(void *data,double x[]);

/* Nelder & Mead simplex optimizer that does not need derivatives.
 *
 * FILE *fp              Write progress output to if fp != NULL
 * void *data            Additional data to be passed to target function
 * nm_target_func func   Target function to minimize 
 * double start[]        Starting values
 * int n                 Number of variables to optimize
 * double toler          Tolerance for convergence (on the simplex size, not chi2)
 * double scale          Scale factor for vertices
 * int maxiter           Max number of iterations
 * double *chi2          Function value at minimum
 *
 * return value: 1 if converged, 0 otherwise.
 */
extern int nmsimplex(FILE *fp,
		     void *data,
		     nm_target_func func,
		     double start[],
		     int n,
		     double toler,
		     double scale,
		     int    maxiter,
		     double *chi2);

#ifdef __cplusplus
}
#endif
#endif

