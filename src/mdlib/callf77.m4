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

define(`fshake_args',`(atom_id iatom[],int *ncon,int *nit,int *maxnit,real dist2[],real xp[],real rij[],real m2[],real invmass[],real tt[],int *error)')

extern void FUNCTION(forshake) fshake_args;

void fshake fshake_args
{
#ifdef USEF77
  FUNCTION(forshake)(iatom,ncon,nit,maxnit,dist2,xp,rij,m2,invmass,tt,error);
#else
  fatal_error(0,"fshake called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

/* LINCS */

define(`flincs_args',`(real *x,real *xp,int *nc,int *ncm,int *cmax,int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,real *blc,real *blcc,real *blm,int *nrec,real *invmass,real *r,real *temp1,real *temp2,real *temp3,real *wangle,int *warn,real *lambda)')

extern void FUNCTION(forlincs) flincs_args;

void flincs flincs_args
{
#ifdef USEF77
  FUNCTION(forlincs)(x,xp,nc,ncm,cmax,bla1,bla2,blnr,blbnb,bllen,blc,blcc,
  	blm,nrec,invmass,r,temp1,temp2,temp3,wangle,warn,lambda);
#else
  fatal_error(0,"flincs called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`flincsld_args',`(real *x,real *xp,int *nc,int *ncm,int *cmax,int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,real *blcc,real *blm,int *nrec,real *r,real *temp1,real *temp2,real *temp3,real *wangle,int *warn)')

extern void FUNCTION(forlincsld) flincsld_args;

void flincsld flincsld_args
{
#ifdef USEF77
  FUNCTION(forlincsld)(x,xp,nc,ncm,cmax,bla1,bla2,blnr,blbnb,bllen,blcc,
  	blm,nrec,r,temp1,temp2,temp3,wangle,warn);
#else
  fatal_error(0,"flincsld called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fconerr_args',`(real *max,real *rms,int *imax,rvec *xprime,int *ncons,int *bla1,int *bla2,real *bllen)')

extern void FUNCTION(forconerr) fconerr_args;

void fconerr fconerr_args
{
#ifdef USEF77
  FUNCTION(forconerr) (max,rms,imax,xprime,ncons,bla1,bla2,bllen);
#else
  fatal_error(0,"fconerr called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

extern void FUNCTION(ffillbuf) (void);

void fillbuf(void)
{
#ifdef USEF77
#ifdef FINVSQRT
  FUNCTION(ffillbuf)();
#endif
#else
  fatal_error(0,"fillbuf called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

/* Matrix diagonalization routine in FORTRAN */ 
define(`fql77_args',`(int *n,real *x, real *d, real *e, int *nmax)')
extern void FUNCTION(forql77) fql77_args;

void fql77 fql77_args
{
#ifdef USEF77
  FUNCTION(forql77)(n,x,d,e,nmax);
#else
  fatal_error(0,"fql77 called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

