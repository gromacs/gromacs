#ifndef _callf77_h
#define _callf77_h

#include "typedefs.h"

/*     Where appropriate,
 *     the LJ parameters are stored in a 1D array nbfp as
 *     c6(1),c12(1),c6(2),c12(2),...,c6(n),c12(n)
 */

extern void f77coul(real *ix,real *iy,real *iz,real *qi,
		    real pos[],int *nj,int jjnr[],
		    real charge[],real faction[],real fip[3],real egcoul[]);

extern void f77ljc(real *ix,real *iy,real *iz,real *qi,
		   real pos[],int *nj,int type[],int jjnr[],
		   real charge[],real nbfp[],
		   real faction[],real fip[],
		   real *egcoul,real *egnb);

extern void f77bham(real *ix,real *iy,real *iz,real *qi,
		    real pos[],int *nj,int type[],int jjnr[],
		    real charge[],real nbfp[],
		    real faction[],real fip[],
		    real *egcoul,real *egnb);

extern void f77water(int  *i0,real xw[],real *eps,
		     real pos[],int *nj,int type[],int jjnr[],
		     real charge[],real nbfp[],
		     real faction[],real fw[],
		     real *egcoul,real *egnb);
			
extern void f77wcoul(int *i0,real xw[],real *eps,
		     real pos[],int *nj,int jjnr[],
		     real charge[],real faction[],real fw[],
		     real egcoul[]);
   
extern void f77free(real *ix,real *iy,real *iz,int *inr,
		    real pos[],int *nj,int jjnr[],
		    int  typeA[],  int typeB[], real *eps,
		    real chargeA[],real chargeB[],
		    real nbfpA[],  real nbfpB[],
		    real faction[],real fip[],
		    real *Vc,      real *Vnb,
		    real *lambda,  real *dvdlambda,
		    real *krf,     real *crf,
		    real *tfac,    real trunctab[]);

extern void f77tab(real *ix,real *iy,real *iz,real *qi,
		   real pos[],int *nj,int type[],t_nl_j jjnr[],
		   real charge[],real nbfp[],
		   real faction[],real fip[],
		   real *Vc,real *Vnb,
		   int  *ntab,real *tabscale,
		   real VFtab[]);

extern void f77coultab(real *ix,real *iy,real *iz,real *qi,
		       real pos[],int *nj,int type[],t_nl_j jjnr[],
		       real charge[],real nbfp[],
		       real faction[],real fip[],
		       real *Vc,real *Vnb,
		       int  *ntab,real *tabscale,
		       real VFtab[]);
		       
/* Initiate invsqrt calculations in fortran */
extern void fillbuf(void);

/* Fortran versions of shake and settle */
extern void fsettle(int *nshake,int owptr[],
		    real b4[],real after[],
		    real *dOH,real *dHH,real *mO,real *mH);
		     
extern void fshake(atom_id iatom[],int *ncon,int *nit,int *maxnit,
		   real dist2[],real xp[],real rij[],real m2[],
		   real invmass[],real tt[],int *error);

/* Fortran routines for LINCS algorithm */ 
extern void flincs(real *x,real *xp,int *nc,int *ncm,int *cmax,
		   int *bla1,int *bla2,int *blnr,
		   int *blbnb,real *bllen,real *blc,real *blcc,real *blm,
		   int *nrec,real *invmass,real *r,real *temp1,
		   real *temp2,real *temp3,real *wangle,int *warn,
		   real *lambda);

extern void flincsld(real *x,real *xp,int *nc,int *ncm,int *cmax,
		     int *bla1,int *bla2,int *blnr,
		     int *blbnb,real *bllen,real *blcc,real *blm,int *nrec,
		     real *r,real *temp1,real *temp2,real *temp3,
		     real *wangle,int *warn);

extern void fconerr(real *max,real *rms,int *imax,
		    rvec *xprime,int *ncons,int *bla1,int *bla2,real *bllen);

extern void fql77(int *n,real *x, real *d, real *e, int *nmax);

#endif
