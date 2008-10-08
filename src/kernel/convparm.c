/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "physics.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "gmx_fatal.h"
#include "topio.h"
#include "toputil.h"
#include "convparm.h"
#include "names.h"
#include "gpp_atomtype.h"

static int round_check(real r,int limit,int ftype,char *name)
{
  int i;

  if (r >= 0)
    i = (int)(r + 0.5);
  else
    i = (int)(r - 0.5);

  if (r-i > 0.01 || r-i < -0.01)
    gmx_fatal(FARGS,"A non-integer value (%f) was supplied for '%s' in %s",
	      r,name,interaction_function[ftype].longname);

  if (i < limit)
    gmx_fatal(FARGS,"Value of '%s' in %s is %d, which is smaller than the minimum of %d",
	      name,interaction_function[ftype].longname,i,limit);

  return i;
}

static void set_ljparams(int comb,real reppow,real v,real w,real *c6,real *c12)
{
  if (comb == eCOMB_ARITHMETIC || comb == eCOMB_GEOM_SIG_EPS) {
    if (v >= 0) {
      *c6  = 4*w*pow(v,6.0);
      *c12 = 4*w*pow(v,reppow);
    } else {
      /* Interpret negative sigma as c6=0 and c12 with -sigma */
      *c6  = 0;
      *c12 = 4*w*pow(-v,reppow);
    }
  } else {
    *c6  = v;
    *c12 = w;
  }
}

static void assign_param(t_functype ftype,t_iparams *new,
			 real old[MAXFORCEPARAM],int comb,real reppow)
{
  int  i,j;
  real tmp;

  /* Set to zero */
  for(j=0; (j<MAXFORCEPARAM); j++)
    new->generic.buf[j]=0.0;
    
  switch (ftype) {
  case F_G96ANGLES:
    /* Post processing of input data: store cosine iso angle itself */
    new->harmonic.rA =cos(old[0]*DEG2RAD);
    new->harmonic.krA=old[1];
    new->harmonic.rB =cos(old[2]*DEG2RAD);
    new->harmonic.krB=old[3];
    break;
  case F_G96BONDS:
    /* Post processing of input data: store square of length itself */
    new->harmonic.rA =sqr(old[0]);
    new->harmonic.krA=old[1];
    new->harmonic.rB =sqr(old[2]);
    new->harmonic.krB=old[3];
    break;
  case F_FENEBONDS:
    new->fene.bm=old[0];
    new->fene.kb=old[1];
    break;
  case F_TABBONDS:
  case F_TABBONDSNC:
  case F_TABANGLES:
  case F_TABDIHS:
    new->tab.table = round_check(old[0],0,ftype,"table index");
    new->tab.kA    = old[1];
    new->tab.kB    = old[3];
    break;
  case F_CROSS_BOND_BONDS:
    new->cross_bb.r1e=old[0];
    new->cross_bb.r2e=old[1];
    new->cross_bb.krr=old[2];
    break;
  case F_CROSS_BOND_ANGLES:
    new->cross_ba.r1e=old[0];
    new->cross_ba.r2e=old[1];
    new->cross_ba.r3e=old[2];
    new->cross_ba.krt=old[3];
    break;
  case F_UREY_BRADLEY:
    new->u_b.theta=old[0];
    new->u_b.ktheta=old[1];
    new->u_b.r13=old[2];
    new->u_b.kUB=old[3];
    break;
  case F_QUARTIC_ANGLES:
    new->qangle.theta=old[0];
    for(i=0; i<5; i++)
      new->qangle.c[i]=old[i+1];
    break;
  case F_ANGLES:
  case F_BONDS:
  case F_HARMONIC:
  case F_IDIHS:
    new->harmonic.rA =old[0];
    new->harmonic.krA=old[1];
    new->harmonic.rB =old[2];
    new->harmonic.krB=old[3];
    break;
  case F_MORSE:
    new->morse.b0    =old[0];
    new->morse.cb    =old[1];
    new->morse.beta  =old[2];
    break;
  case F_CUBICBONDS:
    new->cubic.b0    =old[0];
    new->cubic.kb    =old[1];
    new->cubic.kcub  =old[2];
    break;
  case F_CONNBONDS:
    break;
  case F_POLARIZATION:
    new->polarize.alpha = old[0];
    break;
  case F_WATER_POL:
    new->wpol.al_x   =old[0];
    new->wpol.al_y   =old[1];
    new->wpol.al_z   =old[2];
    new->wpol.rOH    =old[3];
    new->wpol.rHH    =old[4];
    new->wpol.rOD    =old[5];
    break;
  case F_THOLE_POL:
    new->thole.a      = old[0];
    new->thole.alpha1 = old[1];
    new->thole.alpha2 = old[2];
    if ((old[1] > 0) && (old[2] > 0))
      new->thole.rfac = old[0]*pow(old[1]*old[2],-1.0/6.0);
    else
      new->thole.rfac = 1;
    break;
  case F_BHAM:
    new->bham.a = old[0];
    new->bham.b = old[1];
    new->bham.c = old[2];
    break;
  case F_LJ14:
    set_ljparams(comb,reppow,old[0],old[1],&new->lj14.c6A,&new->lj14.c12A);
    set_ljparams(comb,reppow,old[2],old[3],&new->lj14.c6B,&new->lj14.c12B);
    break;
  case F_LJC14_Q:
    new->ljc14.fqq = old[0];
    new->ljc14.qi  = old[1];
    new->ljc14.qj  = old[2];
    set_ljparams(comb,reppow,old[3],old[4],&new->ljc14.c6,&new->ljc14.c12);
    break;
  case F_LJC_PAIRS_NB:
    new->ljcnb.qi = old[0];
    new->ljcnb.qj = old[1];
    set_ljparams(comb,reppow,old[2],old[3],&new->ljcnb.c6,&new->ljcnb.c12);
    break;
  case F_LJ:
    set_ljparams(comb,reppow,old[0],old[1],&new->lj.c6,&new->lj.c12);
    break;
  case F_PDIHS:
  case F_ANGRES:
  case F_ANGRESZ:
    new->pdihs.phiA = old[0];
    new->pdihs.cpA  = old[1];
		  
    /* Dont do any checks if all parameters are zero (such interactions will be removed) */
    tmp=fabs(old[0])+fabs(old[1])+fabs(old[2])+fabs(old[3])+fabs(old[4]);
    new->pdihs.mult = (tmp < GMX_REAL_MIN) ? 0 : round_check(old[2],1,ftype,"multiplicity");
		  
    new->pdihs.phiB = old[3];
    new->pdihs.cpB  = old[4];
    break;
  case F_POSRES:
    new->posres.fcA[XX]   = old[0];
    new->posres.fcA[YY]   = old[1];
    new->posres.fcA[ZZ]   = old[2];
    new->posres.fcB[XX]   = old[3];
    new->posres.fcB[YY]   = old[4];
    new->posres.fcB[ZZ]   = old[5];
    new->posres.pos0A[XX] = old[6];
    new->posres.pos0A[YY] = old[7];
    new->posres.pos0A[ZZ] = old[8];
    new->posres.pos0B[XX] = old[9];
    new->posres.pos0B[YY] = old[10];
    new->posres.pos0B[ZZ] = old[11];
    break;
  case F_DISRES:
    new->disres.label = round_check(old[0],0,ftype,"label");
    new->disres.type  = round_check(old[1],1,ftype,"type'");
    new->disres.low   = old[2];
    new->disres.up1   = old[3];
    new->disres.up2   = old[4];
    new->disres.kfac  = old[5];
    break;
  case F_ORIRES:
    new->orires.ex    = round_check(old[0],1,ftype,"experiment") - 1;
    new->orires.label = round_check(old[1],1,ftype,"label");
    new->orires.power = round_check(old[2],0,ftype,"power");
    new->orires.c     = old[3];
    new->orires.obs   = old[4];
    new->orires.kfac  = old[5];
    break;
  case F_DIHRES:
    new->dihres.label = round_check(old[0],0,ftype,"label");
    new->dihres.phi   = old[1];
    new->dihres.dphi  = old[2];
    new->dihres.kfac  = old[3];
    new->dihres.power = round_check(old[4],0,ftype,"power");
    break;
  case F_RBDIHS:
    for (i=0; (i<NR_RBDIHS); i++) {
      new->rbdihs.rbcA[i]=old[i]; 
      new->rbdihs.rbcB[i]=old[NR_RBDIHS+i]; 
    }
    break;
  case F_FOURDIHS:
    /* Read the dihedral parameters to temporary arrays,
     * and convert them to the computationally faster
     * Ryckaert-Bellemans form.
     */   
    /* Use conversion formula for OPLS to Ryckaert-Bellemans: */
    new->rbdihs.rbcA[0]=old[1]+0.5*(old[0]+old[2]);
    new->rbdihs.rbcA[1]=0.5*(3.0*old[2]-old[0]);
    new->rbdihs.rbcA[2]=4.0*old[3]-old[1];
    new->rbdihs.rbcA[3]=-2.0*old[2];
    new->rbdihs.rbcA[4]=-4.0*old[3];
    new->rbdihs.rbcA[5]=0.0;

    new->rbdihs.rbcB[0]=old[NR_FOURDIHS+1]+0.5*(old[NR_FOURDIHS+0]+old[NR_FOURDIHS+2]);
    new->rbdihs.rbcB[1]=0.5*(3.0*old[NR_FOURDIHS+2]-old[NR_FOURDIHS+0]);
    new->rbdihs.rbcB[2]=4.0*old[NR_FOURDIHS+3]-old[NR_FOURDIHS+1];
    new->rbdihs.rbcB[3]=-2.0*old[NR_FOURDIHS+2];
    new->rbdihs.rbcB[4]=-4.0*old[NR_FOURDIHS+3];
    new->rbdihs.rbcB[5]=0.0;
    break;    
  case F_CONSTR:
  case F_CONSTRNC:
    new->constr.dA = old[0];
    new->constr.dB = old[1];
    break;
  case F_SETTLE:
    new->settle.doh=old[0];
    new->settle.dhh=old[1];
    break;
  case F_VSITE2:
  case F_VSITE3:
  case F_VSITE3FD:
  case F_VSITE3OUT:
  case F_VSITE4FD:
  case F_VSITE4FDN:
    new->vsite.a=old[0];
    new->vsite.b=old[1];
    new->vsite.c=old[2];
    new->vsite.d=old[3];
    new->vsite.e=old[4];
    new->vsite.f=old[5];
    break;
  case F_VSITE3FAD:
    new->vsite.a=old[1] * cos(DEG2RAD * old[0]);
    new->vsite.b=old[1] * sin(DEG2RAD * old[0]);
    new->vsite.c=old[2];
    new->vsite.d=old[3];
    new->vsite.e=old[4];
    new->vsite.f=old[5];
    break;
  case F_VSITEN:
    new->vsiten.n = round_check(old[0],1,ftype,"number of atoms");
    new->vsiten.a = old[1];
    break;
  default:
    gmx_fatal(FARGS,"unknown function type %d in %s line %d",
		ftype,__FILE__,__LINE__);
  }
}

static int enter_params(gmx_ffparams_t *ffparams, t_functype ftype,
			real forceparams[MAXFORCEPARAM],int comb,real reppow,
			int start,bool bAppend)
{
  t_iparams new;
  int       type;
  
  assign_param(ftype,&new,forceparams,comb,reppow);
  if (!bAppend) {
    for (type=start; (type<ffparams->ntypes); type++) {
      if (ffparams->functype[type]==ftype) {
	if (memcmp(&new,&ffparams->iparams[type],(size_t)sizeof(new)) == 0)
	  return type;
      }
    }
  }
  else {
    type = ffparams->ntypes;
  }
  if (debug)
    fprintf(debug,"copying new to ffparams->iparams[%d] (ntypes=%d)\n",
	    type,ffparams->ntypes);
  memcpy(&ffparams->iparams[type],&new,(size_t)sizeof(new));
  
  ffparams->ntypes++;
  ffparams->functype[type]=ftype;

  return type;
}

static void append_interaction(t_ilist *ilist,
                               int type,int nral,atom_id a[MAXATOMLIST])
{
  int i,where1;
  
  where1     = ilist->nr;
  ilist->nr += nral+1;

  ilist->iatoms[where1++]=type;
  for (i=0; (i<nral); i++) 
    ilist->iatoms[where1++]=a[i];
}

static void enter_function(t_params *p,t_functype ftype,int comb,real reppow,
                           gmx_ffparams_t *ffparams,t_ilist *il,
			   int *maxtypes,
			   bool bNB,bool bAppend)
{
  int     k,type,nr,nral,delta,start;
  
  start = ffparams->ntypes;
  nr    = p->nr;
  
  for (k=0; k<nr; k++) {
    if (*maxtypes <= ffparams->ntypes) {
      *maxtypes += 1000;
      srenew(ffparams->functype,*maxtypes);
      srenew(ffparams->iparams, *maxtypes);
      if (debug) 
	fprintf(debug,"%s, line %d: srenewed idef->functype and idef->iparams to %d\n",
		__FILE__,__LINE__,*maxtypes);
    }
    type = enter_params(ffparams,ftype,p->param[k].c,comb,reppow,start,bAppend);
    if (!bNB) {
      nral  = NRAL(ftype);
      delta = nr*(nral+1);
      srenew(il->iatoms,il->nr+delta);
      append_interaction(il,type,nral,p->param[k].a);
    }
  }
}

static void new_interaction_list(t_ilist *ilist)
{
  int i;
  
  ilist->nr=0;
  ilist->iatoms=NULL;
}

void convert_params(int atnr,t_params nbtypes[],
		    t_molinfo *mi,int comb,real reppow,real fudgeQQ,
		    gmx_mtop_t *mtop)
{
  int    i,j,maxtypes,mt;
  unsigned long  flags;
  gmx_ffparams_t *ffp;
  gmx_moltype_t *molt;
  t_params *plist;

  maxtypes=0;
  
  ffp = &mtop->ffparams;
  ffp->ntypes   = 0;
  ffp->atnr     = atnr;
  ffp->functype = NULL;
  ffp->iparams  = NULL;

  enter_function(&(nbtypes[F_LJ]),  (t_functype)F_LJ,    comb,reppow,ffp,NULL,
		 &maxtypes,TRUE,TRUE);
  enter_function(&(nbtypes[F_BHAM]),(t_functype)F_BHAM,  comb,reppow,ffp,NULL,
		 &maxtypes,TRUE,TRUE);

  for(mt=0; mt<mtop->nmoltype; mt++) {
    molt = &mtop->moltype[mt];
    for(i=0; (i<F_NRE); i++) {
      molt->ilist[i].nr     = 0;
      molt->ilist[i].iatoms = NULL;
      
      plist = mi[mt].plist;

      flags = interaction_function[i].flags;
      if ((i != F_LJ) && (i != F_BHAM) && ((flags & IF_BOND) ||
					   (flags & IF_VSITE) ||
					   (flags & IF_CONSTRAINT))) {
	enter_function(&(plist[i]),(t_functype)i,comb,reppow,
		       ffp,&molt->ilist[i],
		       &maxtypes,FALSE,(i == F_POSRES));
      }
    }
  }
  if (debug) {
    fprintf(debug,"%s, line %d: There are %d functypes in idef\n",
	    __FILE__,__LINE__,ffp->ntypes);
  }

  ffp->fudgeQQ = fudgeQQ;
}
