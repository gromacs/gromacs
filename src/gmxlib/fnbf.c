/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_fnbf_c = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "txtdump.h"
#include "smalloc.h"
#include "ns.h"
#include "vec.h"
#include "maths.h"
#include "macros.h"
#include "force.h"
#include "names.h"
#include "main.h"
#include "xvgr.h"
#include "fatal.h"
#include "physics.h"
#include "fnbf.h"
#include "inloop.h"
#include "nrnb.h"
#include "callf77.h"

void do_fnbf(FILE *log,int ftype,t_forcerec *fr,
	     rvec x[],rvec f[],t_mdatoms *mdatoms,
	     real egnb[],real egcoul[],rvec box_size,
	     t_nrnb *nrnb,real lambda,real *dvdlambda,bool bLR)
{
  static bool bFirst = TRUE;
  int      i,itpA,itpB,gid,m,nj,inr,iinr,nri,k;
  rvec     r_i,f_ip,fw[3],xw[3];
  real     qi=0,Vnb,Vc,eps;
  t_nblist *nlist;
  t_nl_i   *nl_i;
  t_nl_j   *nl_j;
  int      *typeA,*typeB;
  real     *chargeA,*chargeB;
  rvec     *svec,*fshift;
  int      nr_ljc,nr_qq,nr_bham,nr_fsum,nr_free;
  bool     bWater,bTab,bWaterTab;
  
  typeA   = mdatoms->typeA;
  chargeA = mdatoms->chargeA;
  typeB   = mdatoms->typeB;
  chargeB = mdatoms->chargeB;

  svec    = fr->shift_vec;
  fshift  = fr->fshift;
  bTab    = ((fr->eeltype != eelCUT) || (fr->vdwtype != evdwCUT));
  eps     = fr->epsfac;
  
  nr_ljc=nr_qq=nr_bham=nr_fsum=nr_free=0;

  switch (ftype) {
  case F_SR:
    nlist = bLR ? fr->coul_lr : fr->coul;
    break;
  case F_DVDL:
    nlist = bLR ? fr->free_lr : fr->free;
    break;
  default:
    nlist = bLR ? fr->vdw_lr  : fr->vdw;
  }

  /* Some macros for "easy" calling of C and Fortran equivalent routines
   * without copying the argument list
   */
#ifdef USEF77
#define SCAL(x) &x
#define FUNC(x) f77##x
#else
#define SCAL(x) x
#define FUNC(x) c_##x
#endif
  
  for(gid=0; (gid<fr->nn); gid++) {
    if (fr->bMask[gid]) {
      nri  = nlist[gid].nri;
      nl_i = nlist[gid].nl_i;
      Vnb  = 0;
      Vc   = 0;
      
      for (i=0; (i<nri); i++) {
	inr      = nl_i[i].i_atom;
	k        = nl_i[i].shift;
	nj       = nl_i[i].nj;
	nl_j     = &(nlist[gid].nl_j[nl_i[i].j_index]);
	itpA     = typeA[inr];
	bWater   = nl_i[i].bWater;
	
	bWaterTab = bWater && bTab;
	if (bFirst && bWaterTab) {
	  fprintf(log,"Using water table routines\n");
	  bFirst = FALSE;
	}
	if (bFirst && bWater && (ftype == F_SR)) {
	  fprintf(log,"Using water coulomb routines\n");
	  bFirst = FALSE;
	}
	
	if (bWater) {
	  for(m=0; (m<3); m++) {
	    rvec_add(x[inr+m],svec[k],xw[m]);
	    clear_rvec(fw[m]);
	  }
	}
	else {
	  qi = chargeA[inr]*eps;
	  rvec_add(x[inr],svec[k],r_i);
	  clear_rvec(f_ip);
	}
	
	switch (ftype) {
	case F_SR:
	  if (bWaterTab) 
	    FUNC(wcoultab)(SCAL(inr),xw[0],SCAL(eps),x[0],SCAL(nj),nl_j,
			   chargeA,f[0],fw[0],&Vc,SCAL(fr->ntab),
			   SCAL(fr->tabscale),fr->VFtab);
	  else if (bTab)
	    FUNC(coultab)(SCAL(r_i[XX]),SCAL(r_i[YY]),SCAL(r_i[ZZ]),SCAL(qi),
			  x[0],SCAL(nj),nl_j,chargeA,f[0],f_ip,&Vc,
			  SCAL(fr->ntab),SCAL(fr->tabscale),fr->VFtab);
	  else if (bWater)
	    FUNC(wcoul)(SCAL(inr),xw[0],SCAL(eps),x[0],SCAL(nj),nl_j,
			chargeA,f[0],fw[0],&Vc);
	  else 
	    FUNC(coul)(SCAL(r_i[XX]),SCAL(r_i[YY]),SCAL(r_i[ZZ]),SCAL(qi),
		       x[0],SCAL(nj),nl_j,chargeA,f[0],f_ip,&Vc);
	  
	  nr_qq+=nj;
	  if (bWater)
	    nr_qq+=2*nj;
	  break;
	  
	case F_LJ:
	  if (bWaterTab) 
	    FUNC(watertab)(SCAL(inr),xw[0],SCAL(eps),x[0],SCAL(nj),typeA,nl_j,
			   chargeA,fr->nbfp[itpA],f[0],fw[0],&Vc,&Vnb,
			   SCAL(fr->ntab),SCAL(fr->tabscale),fr->VFtab);
	  else if (bTab)
	    FUNC(tab)(SCAL(r_i[XX]),SCAL(r_i[YY]),SCAL(r_i[ZZ]),SCAL(qi),
		      x[0],SCAL(nj),typeA,nl_j,
		      chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb,
		      SCAL(fr->ntab),SCAL(fr->tabscale),fr->VFtab);
	  else if (bWater)
	    FUNC(water)(SCAL(inr),xw[0],SCAL(eps),x[0],SCAL(nj),typeA,nl_j,
			chargeA,fr->nbfp[itpA],f[0],fw[0],&Vc,&Vnb);
	  else
	    FUNC(ljc)(SCAL(r_i[XX]),SCAL(r_i[YY]),SCAL(r_i[ZZ]),SCAL(qi),
		      x[0],SCAL(nj),typeA,nl_j,
		      chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb);
	
	  nr_ljc+=nj;
	  if (bWater)
	    nr_qq+=2*nj;
	  break;
	  
	case F_BHAM:
	  if (bTab)
	    c_bhamtab(r_i[XX],r_i[YY],r_i[ZZ],qi,
		      x[0],nj,typeA,nl_j,
		      chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb,
		      fr->ntab,fr->tabscale,fr->tabscale_exp,fr->VFtab);
	  else 
	    FUNC(bham)(SCAL(r_i[XX]),SCAL(r_i[YY]),SCAL(r_i[ZZ]),SCAL(qi),
		       x[0],SCAL(nj),typeA,nl_j,
		       chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb);
	  
	  nr_bham+=nj;
	  break;
	  
	case F_DVDL:
	  if (bWater) 
	    fatal_error(0,"Free energy routines called with water flag.");
	  itpB = typeB[inr];
	  
	  c_free(r_i[XX],r_i[YY],r_i[ZZ],inr,
		 x[0],nj,nl_j,typeA,typeB,eps,
		 chargeA,chargeB,
		 fr->nbfp[itpA],fr->nbfp[itpB],
		 f[0],f_ip,&Vc,&Vnb,lambda,dvdlambda,fr->ntab,
		 fr->tabscale,fr->VFtab);
	  nr_free+=nj;
	  break;
	  
	default:
	  fatal_error(0,"Wrong ftype = %d (%s)",
		      ftype,interaction_function[ftype].longname);
	}
	if (bWater) {
	  /* Copy the forces to force array */
	  for(m=0; (m<3); m++) {
	    rvec_inc(fshift[k],fw[m]);
	    rvec_inc(f[inr+m],fw[m]);
	  }
	  nr_fsum+=18;
	}
	else {
	  rvec_inc(fshift[k],f_ip);
	  rvec_inc(f[inr],f_ip);
	  nr_fsum+=6;
	}
      }
      egcoul[gid] += Vc;
      egnb[gid]   += Vnb;
    }
  }
  if (bTab) {
    inc_nrnb(nrnb,eNR_TAB,nr_ljc);
    inc_nrnb(nrnb,eNR_COULTAB,nr_qq);
  }
  else {
    inc_nrnb(nrnb,eNR_LJC,nr_ljc);
    inc_nrnb(nrnb,eNR_QQ,nr_qq);
    inc_nrnb(nrnb,eNR_BHAM,nr_bham);
  }
  inc_nrnb(nrnb,eNR_FSUM,nr_fsum);
  inc_nrnb(nrnb,eNR_FREE,nr_free);
#undef ARG
#undef FUNC
}

