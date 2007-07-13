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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "gmx_fatal.h"
#include "string2.h"
#include "ebin.h"
#include "main.h"
#include "maths.h"
#include "vec.h"

static real rms_ener(t_energy *e,int nsteps)
{
  return sqrt(e->eav/nsteps);
}

t_ebin *mk_ebin(void)
{
  t_ebin *eb;
  
  snew(eb,1);
  
  return eb;
}

int get_ebin_space(t_ebin *eb,int nener,char *enm[])
{
  int index;
  int i;
  
  index=eb->nener;
  eb->nener+=nener;
  srenew(eb->e,eb->nener);
  srenew(eb->enm,eb->nener);
  for(i=index; (i<eb->nener); i++) {
    eb->e[i].e=0;
    eb->e[i].eav=0;
    eb->e[i].esum=0;
    eb->e[i].e2sum=0;
    eb->enm[i]=strdup(enm[i-index]);
  }
  return index;
}

void add_ebin(t_ebin *eb,int index,int nener,real ener[],int step)
{
  int      i,m;
  double   e,sum,sigma,invmm,diff;
  t_energy *eg;
  
  if ((index+nener > eb->nener) || (index < 0))
    gmx_fatal(FARGS,"%s-%d: Energies out of range: index=%d nener=%d maxener=%d",
		__FILE__,__LINE__,index,nener,eb->nener);
    
  m      = step;
  if (m > 0) 
    invmm = (1.0/(double)m)/((double)m+1.0);
  else
    invmm = 0.0;
    
  eg=&(eb->e[index]);
  
  for(i=0; (i<nener); i++) {
    /* Value for this component */
    e      = ener[i];
    
    /* first update sigma, then sum */
    eg[i].e    = e;
    diff       = eg[i].esum - m*e;
    eg[i].eav  += diff*diff*invmm;
    eg[i].esum += e;
  }
}

void add_ebin_nosum(t_ebin *eb,int index,int nener,real ener[])
{
  int      i;
  t_energy *eg;
  
  if ((index+nener > eb->nener) || (index < 0))
    gmx_fatal(FARGS,"%s-%d: Energies out of range: index=%d nener=%d maxener=%d",
		__FILE__,__LINE__,index,nener,eb->nener);
    
  eg=&(eb->e[index]);
  
  for(i=0; (i<nener); i++) {
    /* Value for this component */
    eg[i].e = ener[i];
  }
}

void pr_ebin(FILE *fp,t_ebin *eb,int index,int nener,int nperline,
	     int prmode,int tsteps,bool bPrHead)
{
  int  i,j,i0;
  real ee=0;
    
  if (index < 0)
    gmx_fatal(FARGS,"Invalid index in pr_ebin: %d",index);
  if (nener == -1)
    nener=eb->nener;
  else
    nener=index+nener;
  for(i=index; (i<nener); ) {
    if (bPrHead) {
      i0=i;
      for(j=0; (j<nperline) && (i<nener); j++,i++)
	fprintf(fp,"%15s",eb->enm[i]);
      fprintf(fp,"\n");
      i=i0;
    }
    for(j=0; (j<nperline) && (i<nener); j++,i++) {
      if (prmode == eprNORMAL)
	ee=eb->e[i].e;
      else if (prmode == eprRMS)
	ee=rms_ener(&(eb->e[i]),tsteps);
      else if (prmode == eprAVER)
	ee=eb->e[i].esum/tsteps;
      else
	gmx_fatal(FARGS,"Invalid print mode %d in pr_ebin",prmode);
      
      fprintf(fp,"   %12.5e",ee);
    }
    fprintf(fp,"\n");
  }
}

#ifdef DEBUGEBIN
int main(int argc,char *argv[])
{
#define NE 12
#define NT 7
#define NS 5

  t_ebin *eb;
  int    i;
  char   buf[25];
  char   *ce[NE],*ct[NT],*cs[NS];
  real   e[NE],t[NT],s[NS];
  int    ie,it,is;
  
  eb=mk_ebin();
  for(i=0; (i<NE); i++) {
    e[i]=i;
    sprintf(buf,"e%d",i);
    ce[i]=strdup(buf);
  }
  ie=get_ebin_space(eb,NE,ce);
  add_ebin(eb,ie,NE,e,0);
  for(i=0; (i<NS); i++) {
    s[i]=i;
    sprintf(buf,"s%d",i);
    cs[i]=strdup(buf);
  }
  is=get_ebin_space(eb,NS,cs);
  add_ebin(eb,is,NS,s,0);
  for(i=0; (i<NT); i++) {
    t[i]=i;
    sprintf(buf,"t%d",i);
    ct[i]=strdup(buf);
  }
  it=get_ebin_space(eb,NT,ct);
  add_ebin(eb,it,NT,t,0);
  
  printf("Normal:\n");
  pr_ebin(stdout,eb,0,-1,5,eprNORMAL,1);

  printf("Average:\n");
  pr_ebin(stdout,eb,ie,NE,5,eprAVER,1);
  pr_ebin(stdout,eb,is,NS,3,eprAVER,1);
  pr_ebin(stdout,eb,it,NT,4,eprAVER,1);

  printf("RMS:\n");
  pr_ebin(stdout,eb,0,-1,5,eprRMS,1);
}
#endif
