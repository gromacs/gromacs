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
 
/********************************************************
 *	fnbf.c IS A GENERATED FILE DO NOT EDIT     	*
 *     		edit fnbf.m4 instead			*
 ********************************************************/
 
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
#include "inner.h"
#include "nrnb.h"
#ifdef USEVECTOR
#include "smalloc.h"

static   real *fbuf=NULL;
#endif	

void do_fnbf(FILE *log,t_forcerec *fr,
	     rvec x[],rvec f[],t_mdatoms *mdatoms,
	     real egnb[],real egcoul[],rvec box_size,
	     t_nrnb *nrnb,real lambda,real *dvdlambda,
	     bool bLR,int eNL)
{
  t_nblist *nlist;
  real     *fshift;
  int      i,i0,i1,nrnb_ind;
  bool     bWater;
  
#ifdef USEVECTOR
  snew(fbuf,mdatoms->nr*3);
#endif  
  if (eNL >= 0) {
    i0 = eNL;
    i1 = i0+1;
  }
  else {
    i0 = 0;
    i1 = eNL_NR;
  }
  if (bLR)
    fshift = fr->fshift_lr[0];
  else
    fshift = fr->fshift[0];
  for(i=i0; (i<i1); i++) {
    if (bLR) 
      nlist  = &(fr->nlist_lr[i]);
    else 
      nlist = &(fr->nlist_sr[i]);
      
    if (nlist->nri > 0) {
      nrnb_ind = nlist->il_code;
      bWater   = (nrnb_ind & (1<<4)) != 0;
      
define(`LJC_ARGS',`SCAL(fr->ntype),mdatoms->typeA,fr->nbfp,egnb')
define(`BHAM_ARGS',`LJC_ARGS')
define(`FEP_ARGS',`mdatoms->chargeB,mdatoms->typeB,SCAL(lambda),dvdlambda')
define(`RF_ARGS',`SCAL(fr->k_rf)')
define(`TAB_ARGS',`SCAL(fr->tabscale),fr->VFtab')
define(`BHTAB_ARGS',`SCAL(fr->tabscale_exp)')
ifdef(`USEVECTOR',`define(`ALL_ARGS',`fbuf,SCAL(nlist->nri),nlist->iinr,nlist->`shift',nlist->gid,nlist->jindex,nlist->jjnr,x[0],fshift,SCAL(fr->epsfac),mdatoms->chargeA,f[0],egcoul,fr->shift_vec[0]')',`define(`ALL_ARGS',`SCAL(nlist->nri),nlist->iinr,nlist->`shift',nlist->gid,nlist->jindex,nlist->jjnr,x[0],fshift,SCAL(fr->epsfac),mdatoms->chargeA,f[0],egcoul,fr->shift_vec[0]')')
			
      switch (nrnb_ind) {
      case eNR_LJC:
	FUNC(ljc)(LJC_ARGS,ALL_ARGS);
	break;
      case eNR_QQ:
	FUNC(coul)(ALL_ARGS);
	break;
      case eNR_BHAM:
	FUNC(bham)(BHAM_ARGS,ALL_ARGS);
	break;
      case eNR_LJCRF:
	FUNC(ljcrf)(LJC_ARGS,RF_ARGS,ALL_ARGS);
	break;
      case eNR_QQRF:
	FUNC(coulrf)(RF_ARGS,ALL_ARGS);
	break;
      case eNR_BHAMRF:
	FUNC(bhamrf)(BHAM_ARGS,RF_ARGS,ALL_ARGS);
	break;
      case eNR_TAB:
	FUNC(ljctab)(LJC_ARGS,TAB_ARGS,ALL_ARGS);
	break;
      case eNR_COULTAB:
	FUNC(coultab)(TAB_ARGS,ALL_ARGS);
	break;
      case eNR_LJC_WAT:
	FUNC(ljcwater)(LJC_ARGS,ALL_ARGS);
	break;
      case eNR_BHAMTAB:
	FUNC(bhamtab)(BHAM_ARGS,TAB_ARGS,BHTAB_ARGS,ALL_ARGS);
	break;
      case eNR_QQ_WAT:
	FUNC(coulwater)(ALL_ARGS);
	break;
      case eNR_BHAM_WAT:
	FUNC(bhamwater)(BHAM_ARGS,ALL_ARGS);
	break;
      case eNR_LJCRF_WAT: 
	FUNC(ljcrfwater)(LJC_ARGS,RF_ARGS,ALL_ARGS);
	break;
      case eNR_QQRF_WAT:
	FUNC(coulrfwater)(RF_ARGS,ALL_ARGS);
	break;
      case eNR_BHAMRF_WAT:
	FUNC(bhamrfwater)(BHAM_ARGS,RF_ARGS,ALL_ARGS);
	break;
      case eNR_TAB_WAT:
	FUNC(ljctabwater)(LJC_ARGS,TAB_ARGS,ALL_ARGS);
	break;
      case eNR_COULTAB_WAT:
	FUNC(coultabwater)(TAB_ARGS,ALL_ARGS);
	break;
      case eNR_BHAMTAB_WAT:
	FUNC(bhamtabwater)(BHAM_ARGS,TAB_ARGS,BHTAB_ARGS,ALL_ARGS);
	break;
      case eNR_BHAM_FREE:
	FUNC(bhamtabfree)(BHAM_ARGS,TAB_ARGS,BHTAB_ARGS,FEP_ARGS,ALL_ARGS);
	break;
      case eNR_LJC_FREE:
	FUNC(ljctabfree)(LJC_ARGS,TAB_ARGS,FEP_ARGS,ALL_ARGS);
	break;
      default:
	fatal_error(0,"No function corresponding to %s in %s `line' %d",
		    nrnb_str(nrnb_ind),__FILE__,__LINE__);
      }
    
      /* Mega flops accounting */
      if (bWater)
        inc_nrnb(nrnb,eNR_INL_IATOM,3*nlist->nri);
      else
        inc_nrnb(nrnb,eNR_INL_IATOM,nlist->nri);
      inc_nrnb(nrnb,nrnb_ind,nlist->nrj);
    }
  }
}

static real dist2(rvec x,rvec y)
{
  rvec dx;
  
  rvec_sub(x,y,dx);
  
  return iprod(dx,dx);
}

static real *mk_14parm(int ntype,int nbonds,t_iatom iatoms[],
		       t_iparams *iparams,int type[])
{
  /* This routine fills a matrix with interaction parameters for
   * 1-4 interaction. It is assumed that these are atomtype dependent
   * only... (but this is checked for...)
   */
  real *nbfp,c6sav,c12sav;
  int  i,ip,ti,tj;
  
  snew(nbfp,2*ntype*ntype);
  for(i=0; (i<nbonds); i+= 3) {
    ip = iatoms[i];
    ti = type[iatoms[i+1]];
    tj = type[iatoms[i+2]];
    c6sav  = C6(nbfp,ntype,ti,tj);
    c12sav = C12(nbfp,ntype,ti,tj);
    C6(nbfp,ntype,ti,tj)  = iparams[ip].lj14.c6A;
    C12(nbfp,ntype,ti,tj) = iparams[ip].lj14.c12A;
    if ((c6sav != 0) || (c12sav != 0)) {
      if ((c6sav  !=  C6(nbfp,ntype,ti,tj)) || 
	  (c12sav != C12(nbfp,ntype,ti,tj))) {
	fatal_error(0,"Force field inconsistency: 1-4 interaction parameters "
		    "for atoms %d-%d not the same as for other atoms "
		    "with the same atom type",iatoms[i+1],iatoms[i+2]);
      }
    }
  }
  return nbfp;
}

real do_14(FILE *log,int nbonds,t_iatom iatoms[],t_iparams *iparams,
	   rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	   matrix box,real lambda,real *dvdlambda,
	   t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  static    real *nbfp14=NULL;
  real      eps;
  real      r2,rtab2;
  int       ai,aj,itype;
  t_iatom   *ia0,*iatom;
  int       gid,shift14;
  int       j_index[] = { 0, 1 };
  int       i1=1,i3=3;

#ifdef USEVECTOR
  snew(fbuf,md->nr*3);
#endif  
  if (nbfp14 == NULL) {
    nbfp14 = mk_14parm(fr->ntype,nbonds,iatoms,iparams,md->typeA);
    if (debug)
      pr_rvec(debug,0,"nbfp14",nbfp14,sqr(fr->ntype));
  }
  shift14 = CENTRAL;
  
  /* Reaction field stuff */  
  eps    = fr->epsfac*fr->fudgeQQ;
  
  rtab2 = sqr(fr->rtab);
    
  ia0=iatoms;

  for(iatom=ia0; (iatom<ia0+nbonds); iatom+=3) {
    itype = iatom[0];
    ai    = iatom[1];
    aj    = iatom[2];
    
    r2    = distance2(x[ai],x[aj]);
    
    if (r2 >= rtab2) {
      fprintf(log,"%d %8.3f %8.3f %8.3f\n",(int)ai+1,
	      x[ai][XX],x[ai][YY],x[ai][ZZ]);
      fprintf(log,"%d %8.3f %8.3f %8.3f\n",(int)aj+1,
	      x[aj][XX],x[aj][YY],x[aj][ZZ]);
      fprintf(log,"1-4 (%d,%d) interaction not within cut-off! r=%g\n",
	      (int)ai+1,(int)aj+1,sqrt(r2));
      /* exit(1); */
    }
    
    gid  = GID(md->cENER[ai],md->cENER[aj],ngrp);
#ifdef DEBUG
    fprintf(log,"LJ14: grp-i=%2d, grp-j=%2d, ngrp=%2d, GID=%d\n",
	    md->cENER[ai],md->cENER[aj],ngrp,gid);
#endif

    if (md->bPerturbed[ai] || md->bPerturbed[aj]) {
      int  tiA,tiB,tjA,tjB;
      real nbfp[18];
      
      /* Save old types */
      tiA = md->typeA[ai];
      tiB = md->typeB[ai];
      tjA = md->typeA[aj];
      tjB = md->typeB[aj];
      md->typeA[ai] = 0;
      md->typeB[ai] = 1;
      md->typeA[aj] = 2;
      md->typeB[aj] = 2;
      
      /* Set nonbonded params */
      C6(nbfp,4,0,2)  = iparams[itype].lj14.c6A;
      C6(nbfp,4,1,2)  = iparams[itype].lj14.c6B;
      C12(nbfp,4,0,2) = iparams[itype].lj14.c12A;
      C12(nbfp,4,1,2) = iparams[itype].lj14.c12B;
      
define(`ARG1',`egnb,SCAL(fr->tabscale),fr->VFtab')
define(`ARG2',`SCAL(i1),&ai,&shift14,&gid,j_index,&aj,x[0],fr->fshift[0],SCAL(eps),md->chargeA,f[0],egcoul,fr->shift_vec[0]')
ifdef(`USEVECTOR',`define(`ARG3',`fbuf,'`ARG2')',`define(`ARG3',`ARG2')')

	     
      FUNC(ljctabfree)(SCAL(i3),md->typeA,nbfp,ARG1,
		       md->chargeB,md->typeB,SCAL(lambda),dvdlambda,ARG3);
		 
      /* Restore old types */
      md->typeA[ai] = tiA;
      md->typeB[ai] = tiB;
      md->typeA[aj] = tjA;
      md->typeB[aj] = tjB;
    }
    else 
      FUNC(ljctab)(SCAL(fr->ntype),md->typeA,nbfp14,ARG1,ARG3);
  }
  return 0.0;
}

