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
static char *SRCID_ebin_c = "$Id$";

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "fatal.h"
#include "string2.h"
#include "ebin.h"
#include "main.h"
#include "maths.h"

static real rms_ener(t_energy *e,int nsteps)
{
  double eav,rms;
  
  eav=e->esum/nsteps;
  rms=sqrt(e->eav/nsteps);

  if (eav != 0.0)  
    if (fabs(rms/eav) < 1e-6)
      rms=0.0;
    
  return rms;
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
  double   e,sum,sigma,invmm;
  t_energy *eg;
  
  if ((index+nener > eb->nener) || (index < 0))
    fatal_error(0,"Energies out of range: index=%d nener=%d",index,nener);
    
  m=step+1;
  invmm=1.0/m;
  invmm/=(m+1);
  eg=&(eb->e[index]);
  
  for(i=0; (i<nener); i++) {
    e      = ener[i];
    sum    = eg[i].esum;
    sigma  = eg[i].eav;
    sum   += e;
    sigma += sqr(sum - m*e)*invmm;
    
    eg[i].e=e;
    eg[i].esum=sum;
    eg[i].eav=sigma;
  }
}

void pr_ebin(FILE *fp,t_ebin *eb,int index,int nener,int nperline,
	     int prmode,int tsteps,bool bPrHead)
{
  int  i,j,i0;
  real ee=0;
    
  if (index < 0)
    fatal_error(0,"Invalid index in pr_ebin: %d",index);
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
	fatal_error(0,"Invalid print mode %d in pr_ebin",prmode);
      
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
