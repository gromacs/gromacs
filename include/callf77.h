#ifndef _callf77_h
#define _callf77_h

#include "typedefs.h"
   
/* Initiate invsqrt calculations in fortran */
extern void fillbuf(void);

/* Fortran versions of shake and settle */
extern void fsettle(int *nshake,int owptr[],
		    real b4[],real after[],
		    real *dOH,real *dHH,real *mO,real *mH,int *error);
		     
extern void fshake(atom_id iatom[],int *ncon,int *nit,int *maxnit,
		   real dist2[],real xp[],real rij[],real m2[],
		   real invmass[],real tt[],int *error);

/* Fortran routines for LINCS algorithm */ 
extern void flincs(real *x,real *xp,int *nc,
		   int *bla1,int *bla2,int *blnr,
		   int *blbnb,real *bllen,real *blc,real *blcc,real *blm,
		   int *nit,int *nrec,real *invmass,real *r,real *temp1,
		   real *temp2,real *temp3,real *wangle,int *warn,
		   real *lambda);

extern void fql77(int *n,real *x, real *d, real *e, int *nmax);

#endif
