#include "typedefs.h"
#include "callf77.h"
#include "fatal.h"

/* This file provides the interface to fortran routines in a machine
 * independent way.
 */

define(`fc_args',`(real *ix,real *iy,real *iz,real *qi,real pos[],int *nj,int jjnr[],real charge[],real faction[],real fip[3],real egcoul[])')
 
extern void FUNCTION(forcoul) fc_args ;

void f77coul fc_args
{
#ifdef USEF77
  FUNCTION(forcoul)(ix,iy,iz,qi,pos,nj,jjnr,charge,faction,fip,egcoul);
#else
  fatal_error(0,"f77coul called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}
  
define(`flj_args',`(real *ix,real *iy,real *iz,real *qi,real pos[],int *nj,int type[],int jjnr[],real charge[],real nbfp[],real faction[],real fip[],real *egcoul,real *egnb)')
extern void FUNCTION(forljc) flj_args;
	    
void f77ljc flj_args
{
#ifdef USEF77
  FUNCTION(forljc)(ix,iy,iz,qi,pos,nj,type,jjnr,
		   charge,nbfp,faction,fip,egcoul,egnb);
#else
  fatal_error(0,"f77ljc called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fbham_args',`(real *ix,real *iy,real *iz,real *qi,real pos[],int *nj,int type[],int jjnr[],real charge[],real nbfp[],real faction[],real fip[],real *egcoul,real *egnb)')
extern void FUNCTION(forbhm) fbham_args;
	    
void f77bham fbham_args
{
#ifdef USEF77
  FUNCTION(forbhm)(ix,iy,iz,qi,pos,nj,type,jjnr,
		   charge,nbfp,faction,fip,egcoul,egnb);
#else
  fatal_error(0,"f77bhm called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fw_args',`(int  *i0,real xw[],real *eps,real pos[],int *nj,int type[],int jjnr[],real charge[],real nbfp[],real faction[],real fw[],real *egcoul,real *egnb)')

extern void FUNCTION(forwater) fw_args;

void f77water fw_args
{
#ifdef USEF77
  FUNCTION(forwater)(i0,xw,eps,pos,nj,type,jjnr,charge,nbfp,faction,fw,
		     egcoul,egnb);
#else
  fatal_error(0,"f77water called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fwtab_args',`(int  *i0,real xw[],real *eps,real pos[],int *nj,int type[],int jjnr[],real charge[],real nbfp[],real faction[],real fw[],real *egcoul,real *egnb,int *ntab,real *tabscale,real VFtab[])')

extern void FUNCTION(forwatertab) fwtab_args;

void f77watertab fwtab_args
{
#ifdef USEF77
  FUNCTION(forwatertab)(i0,xw,eps,pos,nj,type,jjnr,charge,nbfp,faction,fw,
		     egcoul,egnb,ntab,tabscale,VFtab);
#else
  fatal_error(0,"f77watertab called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fwc_args',`(int *i0,real xw[],real *eps,
	      real pos[],int *nj,int jjnr[],
	      real charge[],real faction[],real fw[],
	      real egcoul[])')
			
extern void FUNCTION(forwcoul) fwc_args;

void f77wcoul fwc_args
{
#ifdef USEF77
  FUNCTION(forwcoul)(i0,xw,eps,pos,nj,jjnr,charge,faction,fw,egcoul);
#else
  fatal_error(0,"f77wcoul called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fwcoultab_args',`(int  *i0,real xw[],real *eps,real pos[],int *nj,int jjnr[],real charge[],real faction[],real fw[],real *egcoul,int *ntab,real *tabscale,real VFtab[])')

extern void FUNCTION(forwcoultab) fwcoultab_args;

void f77wcoultab fwcoultab_args
{
#ifdef USEF77
  FUNCTION(forwcoultab)(i0,xw,eps,pos,nj,jjnr,charge,faction,fw,
		     egcoul,ntab,tabscale,VFtab);
#else
  fatal_error(0,"f77wcoultab called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`ff_args',`(real *ix,real *iy,real *iz,int *inr,real pos[],int *nj,int jjnr[],int  typeA[],  int typeB[], real *eps,real chargeA[],real chargeB[],real nbfpA[],  real nbfpB[],real faction[],real fip[],real *Vc,      real *Vnb,real *lambda,  real *dvdlambda,real *krf,     real *crf,real *tfac,    real trunctab[])')  
   
extern void FUNCTION(forfree) ff_args;

void f77free ff_args
{
#ifdef USEF77
  FUNCTION(forfree)(ix,iy,iz,inr,pos,nj,jjnr,typeA,typeB,eps,
		    chargeA,chargeB,nbfpA,nbfpB,faction,fip,
		    Vc,Vnb,lambda,dvdlambda,krf,crf,tfac,trunctab);
#else
  fatal_error(0,"f77free called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`ft_args',`(real *ix,real *iy,real *iz,real *qi,real pos[],int *nj,int type[],t_nl_j jjnr[],real charge[],real nbfp[],real faction[],real fip[],real *Vc,real *Vnb,int  *ntab,real *tabscale,real VFtab[])')

extern void FUNCTION(fortab) ft_args;

void f77tab ft_args
{
#ifdef USEF77
  FUNCTION(fortab)(ix,iy,iz,qi,
		   pos,nj,type,jjnr,charge,nbfp,
		   faction,fip,Vc,Vnb,ntab,tabscale,VFtab);
#else
  fatal_error(0,"f77tab called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

define(`fct_args',`(real *ix,real *iy,real *iz,real *qi,real pos[],int *nj,t_nl_j jjnr[],real charge[],real faction[],real fip[],real *Vc,int  *ntab,real *tabscale,real VFtab[])')

extern void FUNCTION(forcoultab) fct_args;

void f77coultab fct_args
{
#ifdef USEF77
  FUNCTION(forcoultab)(ix,iy,iz,qi,pos,nj,jjnr,charge,
		       faction,fip,Vc,ntab,tabscale,VFtab);
#else
  fatal_error(0,"f77coultab called (Fortran routine from %s %d)",__FILE__,__LINE__);
#endif
}

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

