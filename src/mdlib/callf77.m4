#include "typedefs.h"
#include "callf77.h"
#include "fatal.h"

/* This file provides the interface to fortran routines in a machine
 * independent way.
 */

/* Fortran versions of shake and settle */
define(`fsettle_args',`(int *nshake,int owptr[],real b4[],real after[],real *dOH,real *dHH,real *mO,real *mH,int *error)')

extern void FUNCTION(forsettle) fsettle_args;
void fsettle fsettle_args
{
#ifdef USEF77
  FUNCTION(forsettle) (nshake,owptr,b4,after,dOH,dHH,mO,mH,error);
#else
  fatal_error(0,"fsettle called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fshake_args',`(atom_id iatom[],int *ncon,int *nit,int *maxnit,real dist2[],real xp[],real rij[],real m2[],real invmass[],real tt[],real lambda[],int *error)')

extern void FUNCTION(forshake) fshake_args;

void fshake fshake_args
{
#ifdef USEF77
  FUNCTION(forshake)(iatom,ncon,nit,maxnit,dist2,xp,rij,m2,invmass,tt,lambda,error);
#else
  fatal_error(0,"fshake called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

/* LINCS */

define(`flincsp_args',`(real *x,real *f,real *fp,int *nc,int *bla1,int *bla2,int *blnr,int *blbnb,real *blc,real *blcc,real *blm,int *nrec,real *invmass,real *r,real *temp1,real *temp2,real *temp3)')

extern void FUNCTION(forlincsp) flincsp_args;

void flincsp flincsp_args
{
#ifdef USEF77
  FUNCTION(forlincsp)(x,f,fp,nc,bla1,bla2,blnr,blbnb,blc,blcc,
  	blm,nrec,invmass,r,temp1,temp2,temp3);
#else
  fatal_error(0,"flincsp called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`flincs_args',`(real *x,real *xp,int *nc,int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,real *blc,real *blcc,real *blm,int *nit,int *nrec,real *invmass,real *r,real *temp1,real *temp2,real *temp3,real *wangle,int *warn,real *lambda)')

extern void FUNCTION(forlincs) flincs_args;

void flincs flincs_args
{
#ifdef USEF77
  FUNCTION(forlincs)(x,xp,nc,bla1,bla2,blnr,blbnb,bllen,blc,blcc,
  	blm,nit,nrec,invmass,r,temp1,temp2,temp3,wangle,warn,lambda);
#else
  fatal_error(0,"flincs called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}
