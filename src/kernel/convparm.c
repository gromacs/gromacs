/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_convparm_c = "$Id$";
#include <math.h>
#include "sysstuff.h"
#include "physics.h"
#include "vec.h"
#include "assert.h"
#include "smalloc.h"
#include "typedefs.h"
#include "fatal.h"
#include "topio.h"
#include "toputil.h"
#include "convparm.h"
#include "names.h"

static void assign_param(t_functype ftype,t_iparams *new,
			 real old[MAXFORCEPARAM])
{
  int i,j;
  
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
  case F_WPOL:
    new->wpol.kx     =old[0];
    new->wpol.ky     =old[1];
    new->wpol.kz     =old[2];
    new->wpol.rOH    =old[3];
    new->wpol.rHH    =old[4];
    new->wpol.rOD    =old[5];
    break;
  case F_BHAM:
    new->bham.a = old[0];
    new->bham.b = old[1];
    new->bham.c = old[2];
    break;
  case F_LJ14:
    new->lj14.c6A  = old[0]; 
    new->lj14.c12A = old[1];
    new->lj14.c6B  = old[2]; 
    new->lj14.c12B = old[3];
    break;
  case F_LJ:
    new->lj.c6  = old[0]; 
    new->lj.c12 = old[1];
    break;
  case F_PDIHS:
  case F_ANGRES:
  case F_ANGRESZ:
    new->pdihs.phiA=old[0];
    new->pdihs.cpA =old[1];
    new->pdihs.mult=old[2];
    new->pdihs.phiB=old[3];
    new->pdihs.cpB =old[4];
    break;
  case F_POSRES:
    new->posres.fc[XX]   = old[0];
    new->posres.fc[YY]   = old[1];
    new->posres.fc[ZZ]   = old[2];
    new->posres.pos0[XX] = old[3];
    new->posres.pos0[YY] = old[4];
    new->posres.pos0[ZZ] = old[5];
    break;
  case F_DISRES:
    new->disres.label = old[0];
    new->disres.type  = old[1];
    new->disres.low   = old[2];
    new->disres.up1   = old[3];
    new->disres.up2   = old[4];
    new->disres.kfac  = old[5];
    break;
  case F_ORIRES:
    if (old[0] < 0)
      fatal_error(0,"Found experiment number for orientation restraints which is smaller than 1 (%d)",old[0]);
    new->orires.ex    = old[0] - 1;
    new->orires.label = old[1];
    new->orires.power = old[2];
    new->orires.c     = old[3];
    new->orires.obs   = old[4];
    new->orires.kfac  = old[5];
    break;
  case F_RBDIHS:
    for (i=0; (i<NR_RBDIHS); i++) 
      new->rbdihs.rbc[i]=old[i]; 
    break;
  case F_SHAKE:
  case F_SHAKENC:
    new->shake.dA = old[0];
    new->shake.dB = old[1];
    break;
  case F_SETTLE:
    new->settle.doh=old[0];
    new->settle.dhh=old[1];
    break;
  case F_DUMMY2:
  case F_DUMMY3:
  case F_DUMMY3FD:
  case F_DUMMY3OUT:
  case F_DUMMY4FD:
    new->dummy.a=old[0];
    new->dummy.b=old[1];
    new->dummy.c=old[2];
    new->dummy.d=old[3];
    new->dummy.e=old[4];
    new->dummy.f=old[5];
    break;
  case F_DUMMY3FAD:
    new->dummy.a=old[1] * cos(DEG2RAD * old[0]);
    new->dummy.b=old[1] * sin(DEG2RAD * old[0]);
    new->dummy.c=old[2];
    new->dummy.d=old[3];
    new->dummy.e=old[4];
    new->dummy.f=old[5];
    break;
  default:
    fatal_error(0,"unknown function type %d in %s line %d",
		ftype,__FILE__,__LINE__);
  }
}

static int enter_params(t_idef *idef, t_functype ftype,
			real forceparams[MAXFORCEPARAM],int start,bool bAppend)
{
  t_iparams new;
  int       type;
  
  assign_param(ftype,&new,forceparams);
  if (!bAppend) {
    for (type=start; (type<idef->ntypes); type++) {
      if (idef->functype[type]==ftype) {
	if (memcmp(&new,&idef->iparams[type],(size_t)sizeof(new)) == 0)
	  return type;
      }
    }
  }
  else
    type=idef->ntypes;
  if (debug)
    fprintf(debug,"copying new to idef->iparams[%d] (ntypes=%d)\n",
	    type,idef->ntypes);
  memcpy(&idef->iparams[type],&new,(size_t)sizeof(new));
  
  idef->ntypes++;
  idef->functype[type]=ftype;

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

static void enter_function(t_params *p,t_functype ftype,
                           t_idef *idef,int *maxtypes,bool bNB,bool bAppend)
{
  int     k,type,nr,nral,delta,start;
  t_ilist *il;
  
  il    = &(idef->il[ftype]);
  start = idef->ntypes;
  nr    = p->nr;
  nral  = NRAL(ftype);
  delta = nr*(nral+1);
  srenew(il->iatoms,il->nr+delta);
  
  for (k=0; k<nr; k++) {
    if (*maxtypes <= idef->ntypes) {
      *maxtypes += 1000;
      srenew(idef->functype,*maxtypes);
      srenew(idef->iparams, *maxtypes);
      if (debug) 
	fprintf(debug,"%s, line %d: srenewed idef->functype and idef->iparams to %d\n",
		__FILE__,__LINE__,*maxtypes);
    }
    type = enter_params(idef,ftype,p->param[k].c,start,bAppend);
    if (!bNB)
      append_interaction(il,type,nral,p->param[k].a);
  }
}

static void new_interaction_list(t_ilist *ilist)
{
  int i;
  
  ilist->nr=0;
  for(i=0; (i<MAXNODES); i++) 
    ilist->multinr[i]=0;
  ilist->iatoms=NULL;
}

void convert_params(int atnr,t_params nbtypes[],
		    t_params plist[],t_idef *idef)
{
  int    i,j,maxtypes;
  unsigned long  flags;

  maxtypes=0;
  
  idef->ntypes   = 0;
  idef->atnr     = atnr;
  idef->nodeid   = 0;  
  idef->functype = NULL;
  idef->iparams  = NULL;
  for(i=0; (i<F_NRE); i++) {
    idef->il[i].nr=0;
    idef->il[i].iatoms=NULL;
  }
  enter_function(&(nbtypes[F_LJ]),  (t_functype)F_LJ,  idef,
		 &maxtypes,TRUE,TRUE);
  enter_function(&(nbtypes[F_BHAM]),(t_functype)F_BHAM,idef,
		 &maxtypes,TRUE,TRUE);
  enter_function(&(plist[F_POSRES]),(t_functype)F_POSRES,idef,
		 &maxtypes,FALSE,TRUE);
		 
  for(i=0; (i<F_NRE); i++) {
    flags = interaction_function[i].flags;
    if ((i != F_LJ) && (i != F_BHAM) && (i != F_POSRES) &&
	((flags & IF_BOND) || (flags & IF_DUMMY) || (flags & IF_CONSTRAINT)))
      enter_function(&(plist[i]),(t_functype)i,idef,&maxtypes,FALSE,FALSE);
  }
  if (debug)
    fprintf(debug,"%s, line %d: There are %d functypes in idef\n",
	    __FILE__,__LINE__,idef->ntypes);
  for(j=0; (j<F_NRE); j++) {
    for (i=0; (i<MAXNODES); i++) 
      idef->il[j].multinr[i]=idef->il[j].nr;
    
    if (idef->il[j].nr > 0)
      printf("# %10s:   %d\n",interaction_function[j].name,idef->il[j].nr);
  }
}
