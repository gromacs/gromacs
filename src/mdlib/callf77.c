#include "typedefs.h"
#include "callf77.h"
#include "fatal.h"

/* This file provides the interface to fortran routines in a machine
 * independent way.
 */

/* Fortran versions of shake and settle */

#ifdef USE_FORTRAN
extern void F77_FUNC(forsettle,FORSETTLE)(int *nshake,int owptr[],real
					  b4[],real after[],real
					  *dOH,real *dHH,real *mO,real
					  *mH,int *error);
extern void F77_FUNC(forshake,FORSHAKE)(atom_id iatom[],int *ncon,
					int *nit, int *maxnit,
					real dist2[],real xp[],
					real rij[],real m2[],
					real invmass[],real tt[],
					real lambda[],int *error);
extern void F77_FUNC(forlincs,FORLINCS)(real *x,real *xp,int *nc,
					int *bla1,int *bla2,int *blnr,
					int *blbnb,real *bllen,
					real *blc,real *blcc,real *blm,
					int *nit,int *nrec,real *invmass,
					real *r,real *temp1,real *temp2,
					real *temp3,real *wangle,
					int *warn,real *lambda);
#endif 

void fsettle(int *nshake,int owptr[],real b4[],real after[],real *dOH,real *dHH,real *mO,real *mH,int *error)
{
#ifdef USE_FORTRAN
  F77_FUNC(forsettle,FORSETTLE) (nshake,owptr,b4,after,dOH,dHH,mO,mH,error);
#else
  fatal_error(0,"fsettle called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

void fshake(atom_id iatom[],int *ncon,int *nit,int *maxnit,real dist2[],real xp[],real rij[],real m2[],real invmass[],real tt[],real lambda[],int *error)
{
#ifdef USE_FORTRAN
  F77_FUNC(forshake,FORSHAKE)(iatom,ncon,nit,maxnit,dist2,xp,rij,m2,invmass,tt,lambda,error);
#else
  fatal_error(0,"fshake called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

/* LINCS */

void flincs(real *x,real *xp,int *nc,int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,real *blc,real *blcc,real *blm,int *nit,int *nrec,real *invmass,real *r,real *temp1,real *temp2,real *temp3,real *wangle,int *warn,real *lambda)
{
#ifdef USE_FORTRAN
  F77_FUNC(forlincs,FORLINCS)(x,xp,nc,bla1,bla2,blnr,blbnb,bllen,blc,blcc,
  	blm,nit,nrec,invmass,r,temp1,temp2,temp3,wangle,warn,lambda);
#else
  fatal_error(0,"flincs called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}
