#ifndef _poisson_h
#define _poisson_h

#include "typedefs.h"

#define llim2 (-3)
#define ulim2  (3)

/* typedef for poisson solver */
typedef struct {
  int  nx,ny,nz;
  real ***ptr;
} t_PSgrid;

extern void unpack_PSgrid(t_PSgrid *grid,int *nx,int *ny,int *nz,real ****ptr);


extern void calc_nxyz(int nx,int ny,int nz,
		      int **nnx,int **nny,int **nnz);
/* Calculate tables to comput modulo (instead of function call) */
		      
extern real ps_gather_f(FILE *log,bool bVerbose,
			int natoms,rvec x[],rvec f[],real charge[],rvec box,
			real pot[],t_PSgrid *grid,rvec beta,t_nrnb *nrnb);

extern void spread_q_poisson(FILE *log,bool bVerbose,
			     int natoms,rvec x[],real charge[],rvec box,
			     real r1,real rc,t_PSgrid *grid,t_nrnb *nrnb,
			     int ntab,real sftab[],real tabfactor);
			     
extern void solve_poisson(FILE *log,t_PSgrid *pot,t_PSgrid *rho,
			  bool bVerbose,t_nrnb *nrnb,int maxnit,real tol,
			  rvec box);

static void calc_invh_h(rvec box,int nx,int ny,int nz,rvec invh,rvec h)
{
  invh[XX] = nx/box[XX];
  invh[YY] = ny/box[YY];
  invh[ZZ] = nz/box[ZZ];
  h[XX]    = 1.0/invh[XX];
  h[YY]    = 1.0/invh[YY];
  h[ZZ]    = 1.0/invh[ZZ];
}

#endif
