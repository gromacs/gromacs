/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Good ROcking Metal Altar for Chronical Sinners
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

void do_fnbf(FILE *log,int ftype,t_forcerec *fr,
	     rvec x[],rvec f[],t_mdatoms *mdatoms,
	     real egnb[],real egcoul[],rvec box_size,
	     t_nrnb *nrnb,real lambda,real *dvdlambda)
{
  int      i,j,itpA,itpB,gid,m,nj,inr,iinr,nri,k;
  rvec     r_i,f_ip,fw[3],xw[3];
  real     qi,Vnb,Vc,eps;
  t_nblist *nlist;
  t_nl_i   *nl_i;
  t_nl_j   *nl_j;
  int      *typeA,*typeB;
  real     *chargeA,*chargeB;
  rvec     *svec,*fshift;
  int      nr_ljc,nr_qq,nr_bham,nr_fsum,nr_free;
  bool     bWater,bTab;
  
  typeA   = mdatoms->typeA;
  chargeA = mdatoms->chargeA;
  typeB   = mdatoms->typeB;
  chargeB = mdatoms->chargeB;

  svec    = fr->shift_vec;
  fshift  = fr->fshift;
  bTab    = (fr->eeltype != eelTWIN);
  eps     = fr->epsfac;
  
  nr_ljc=nr_qq=nr_bham=nr_fsum=nr_free=0;

  switch (ftype) {
  case F_SR:
    nlist=fr->coul;
    break;
  case F_DVDL:
    nlist=fr->free;
    break;
  default:
    nlist=fr->vdw;
  }
  
  for(gid=0; (gid<fr->nn); gid++) {
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
      
      if (bWater && bTab)
	fatal_error(0,"Can't have water loop with tables\n");
	
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
#ifdef USEF77
	if (bTab)
	  CALLF77(forcoultab)(&r_i[XX],&r_i[YY],&r_i[ZZ],&qi,
			      x[0],&nj,typeA,nl_j,
			      chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb,
			      &fr->ntab,&fr->tabscale,fr->VFtab);
	else 
	  CALLF77(forcoul)(&r_i[XX],&r_i[YY],&r_i[ZZ],&qi,
			   x[0],&nj,nl_j,
			   chargeA,f[0],f_ip,&Vc);
#else
	if (bTab) 
	  c_coultab(r_i[XX],r_i[YY],r_i[ZZ],qi,
		    x[0],nj,typeA,nl_j,
		    chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb,
		    fr->ntab,fr->tabscale,fr->VFtab);
	else
	  c_coul(r_i[XX],r_i[YY],r_i[ZZ],qi,
		 x[0],nj,nl_j,chargeA,f[0],f_ip,&Vc);
#endif
	nr_qq+=nj;
	break;
	
      case F_LJ:
#ifdef USEF77
	if (bTab) 
	  CALLF77(fortab)(&r_i[XX],&r_i[YY],&r_i[ZZ],&qi,
			  x[0],&nj,typeA,nl_j,
			  chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb,
			  &fr->ntab,&fr->tabscale,fr->VFtab);
	else if (bWater)
	  CALLF77(forwater)(&inr,xw[0],&eps,x[0],&nj,typeA,nl_j,
			    chargeA,fr->nbfp[itpA],f[0],fw[0],&Vc,&Vnb);
	else
	  CALLF77(forljc)(&r_i[XX],&r_i[YY],&r_i[ZZ],&qi,
			  x[0],&nj,typeA,nl_j,
			  chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb);
#else
	if (bTab)
	  c_tab(r_i[XX],r_i[YY],r_i[ZZ],qi,
		x[0],nj,typeA,nl_j,
		chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb,
		fr->ntab,fr->tabscale,fr->VFtab);
	else if (bWater)
	  c_water(inr,xw[0],eps,x[0],nj,typeA,nl_j,
		  chargeA,fr->nbfp[itpA],f[0],fw[0],&Vc,&Vnb);
	else
	  c_ljc(r_i[XX],r_i[YY],r_i[ZZ],qi,
		x[0],nj,typeA,nl_j,
		chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb);
#endif
	nr_ljc+=nj;
	if (bWater)
	  nr_qq+=2*nj;
	break;
	
      case F_BHAM:
	/* Call C-routine for the time being */
	c_bham(r_i[XX],r_i[YY],r_i[ZZ],qi,
	       x[0],nj,typeA,nl_j,
	       chargeA,fr->nbfp[itpA],f[0],f_ip,&Vc,&Vnb);
	nr_bham+=nj;
	break;
	
      case F_DVDL:
	if (bWater) {
	  /* Call the free enrgy routine three times with different
	   * particles (OW, HW1, HW2)
	   */
	  for(m=0; (m<3); m++) {
	    iinr = inr+m;
	    c_free(xw[m][XX],xw[m][YY],xw[m][ZZ],iinr,
		   x[0],nj,nl_j,typeA,typeB,eps,
		   chargeA,chargeB,
		   fr->nbfp[typeA[iinr]],fr->nbfp[typeB[iinr]],
		   f[0],fw[m],&Vc,&Vnb,lambda,dvdlambda,fr->ntab,
		   fr->tabscale,fr->VFtab);
	    nr_free+=nj;
	  }
	}
	else {
	  itpB     = typeB[inr];
	  c_free(r_i[XX],r_i[YY],r_i[ZZ],inr,
		 x[0],nj,nl_j,typeA,typeB,eps,
		 chargeA,chargeB,
		 fr->nbfp[itpA],fr->nbfp[itpB],
		 f[0],f_ip,&Vc,&Vnb,lambda,dvdlambda,fr->ntab,
		 fr->tabscale,fr->VFtab);
	  nr_free+=nj;
	}
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
}

void fdo_flr(FILE *log,int nri,atom_id i_atoms[],int shift,
	     int njcg,atom_id jcg[],atom_id index[],atom_id acg[],
	     rvec x[],real egcoul[],t_mdatoms *md,int ngener,
	     t_forcerec *fr)
{
#define LR_INC 1024
  static   t_nl_j   **nlj=NULL;
  static   int      *nj,*maxnj;
  atom_id  ip,jp;
  int      i,j,k,k0,k1,jj,m;
  rvec     f_ip,r_i,fw[3],xw[3];
  real     qi,epsje,eps;
  int      igid,jgid,gid,usegid,iaa,nr_fsum;
  rvec     *fshift,*flr,*sv;
  ushort   *cENER;
  int      *type;
  real     Vc,tx,ty,tz;
  real     *charge;
  
  /* Copy some pointers... */
  fshift = fr->fshift_lr;
  flr    = fr->flr;
  sv     = fr->shift_vec;
  cENER  = md->cENER;
  type   = md->typeA;
  charge = md->chargeA;
  
  nr_fsum=0;
  
  if (nlj == NULL) {
    /* Allocate memory for long range */
    snew(nlj,fr->nn);
    snew(nj,fr->nn);
    snew(maxnj,fr->nn);
  }
  else {
    for(j=0; (j<fr->nn); j++)
      nj[j]=0;
  }

  /* Dielectric constant */
  eps  = ONE_4PI_EPS0/sqrt(fr->epsilon_r);
    
  /* Fill the buffer with atoms & group numbers */
  ip   = i_atoms[0];
  igid = cENER[ip];
  
  for(j=0; (j<njcg); j++) {
    k0=index[jcg[j]];
    k1=index[jcg[j]+1];
    for(k=k0; (k<k1); k++) {
      jp           = acg[k];
      jgid         = cENER[jp];
      gid          = GID(igid,jgid,ngener);
      
      if (maxnj[gid] <= nj[gid]) {
	maxnj[gid] += LR_INC;
	fprintf(log,"Increasing buffer size for LR[%d] to %d\n",
		gid,maxnj[gid]);
	srenew(nlj[gid],maxnj[gid]);
      }
      
      nlj[gid][nj[gid]++] = jp;
    }
  }

  for (gid=0; (gid<fr->nn); gid++) {
    Vc     = 0;
    usegid = gid;
    
    if ((type[ip] == fr->nWater) && (nri == 3)) {
      /* Add shift vector to x[ip] to get proper image */
      for(m=0; (m<3); m++)
	rvec_add(x[i_atoms[m]],sv[shift],xw[m]);
      
#ifdef USEF77
      iaa=i_atoms[0];
      CALLF77(forwcoul) (&(iaa),xw[0],&eps,x[0],&nj[gid],nlj[gid],
			charge,flr[0],fw[0],&Vc);
#else
      c_wcoul(i_atoms[0],xw[0],eps,x[0],nj[gid],nlj[gid],
	      charge,flr[0],fw[0],&Vc);
#endif
      
      /* Update force for ip particle and corresponding shift */
      for(j=0; (j<3); j++) {
	ip=i_atoms[j];
	
	tx = fw[j][XX];
	ty = fw[j][YY];
	tz = fw[j][ZZ];
	
	fshift[shift][XX] += tx;
	flr[ip][XX]       += tx;
	fshift[shift][YY] += ty;
	flr[ip][YY]       += ty;
	fshift[shift][ZZ] += tz;
	flr[ip][ZZ]       += tz;
      }
    }
    else { 
      for(i=0; (i<nri); i++) {
	ip = i_atoms[i];
	qi = charge[ip]*eps;
	
	/* Add shift vector to x[ip] to get proper image */
	rvec_add(x[ip],sv[shift],r_i);
	
	if (igid != cENER[ip]) {
	  /* We have multiple energy groups within charge group i.
	   * However, since we have sorted the neighbourlists according
	   * to energy group, we just have to recalculate the entry
	   * in the energy array where this energy contribution is
	   * going to go.
	   */
	  usegid = GID(igid,cENER[nlj[gid][0]],ngener);
	}
	/* Reset force vector for particle ip */
	clear_rvec(f_ip);
#ifdef USEF77
	CALLF77(forcoul) (&r_i[XX],&r_i[YY],&r_i[ZZ],&qi,
			 x[0],&nj[gid],nlj[gid],
			 charge,flr[0],f_ip,&Vc);
#else
	c_coul(r_i[XX],r_i[YY],r_i[ZZ],qi,
	       x[0],nj[gid],nlj[gid],
	       charge,flr[0],f_ip,&Vc);
#endif
	/* Update force for ip particle and corresponding shift */
	for(m=0; (m<DIM); m++) {
	  tx = f_ip[m];
	  fshift[shift][m] += tx;
	  flr[ip][m]       += tx;
	}
      }
    }
    egcoul[usegid] += Vc;
  }
  for(j=0; (j<fr->nn); j++)
    fr->nlr+=nri*nj[j];
}
