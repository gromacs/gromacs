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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_add_par_c = "$Id$";

#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "grompp.h"
#include "pdb2gmx.h"

void add_param(t_params *ps,int ai,int aj, real *c)
{
  int i;
  
  if ((ai < 0) || (aj < 0)) 
    fatal_error(0,"Trying to add impossible atoms: ai=%d, aj=%d",ai,aj);
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=-1;
  ps->param[ps->nr-1].AL=-1;
  if (c==NULL) 
    for(i=0; (i<MAXFORCEPARAM); i++)
      ps->param[ps->nr-1].c[i]=NOTSET;
  else
    for(i=0; (i<MAXFORCEPARAM); i++)
      ps->param[ps->nr-1].c[i]=c[i];
}

void add_imp_param(t_params *ps,int ai,int aj,int ak,int al,real c0, real c1)
{
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=ak;
  ps->param[ps->nr-1].AL=al;
  ps->param[ps->nr-1].C0=c0;
  ps->param[ps->nr-1].C1=c1;
}

void add_dum1_param(t_params *ps,int ai,int aj,int ak,real a)
{
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=ak;
  ps->param[ps->nr-1].C0=a;
}

void add_dum2_param(t_params *ps,int ai,int aj,int ak,int al,real a, real b)
{
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=ak;
  ps->param[ps->nr-1].AL=al;
  ps->param[ps->nr-1].C0=a;
  ps->param[ps->nr-1].C1=b;
}

void add_dum3_param(t_params *ps,int ai,int aj,int ak,int al,
		    real a, real b, real c)
{
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=ak;
  ps->param[ps->nr-1].AL=al;
  ps->param[ps->nr-1].C0=a;
  ps->param[ps->nr-1].C1=b;
  ps->param[ps->nr-1].C2=c;
}

int search_jtype(t_restp *rtp,char *name,bool bNterm)
{
  int  j,k,mn,kmax,jmax;
  char *rtpname,searchname[12];
  
  strcpy(searchname,name);
  
  /* Do a best match comparison */
  /* for protein N-terminus, rename H1, H2 and H3 to H */
  if ( bNterm && (strlen(searchname)==2) && (searchname[0] == 'H') && 
       ( (searchname[1] == '1') || (searchname[1] == '2') || 
	 (searchname[1] == '3') ) )
    searchname[1]='\0';
  kmax=0;
  jmax=-1;
  for(j=0; (j<rtp->natom); j++) {
    rtpname=*(rtp->atomname[j]);
    if (strcasecmp(searchname,rtpname) == 0) {
      jmax=j;
      kmax=strlen(searchname);
      break;
    }
    mn=min((int)strlen(searchname),(int)strlen(rtpname));
    for(k=0; (k<mn); k++) 
      if (searchname[k] != rtpname[k])
	break;
    if (k > kmax) {
      kmax=k;
      jmax=j;
    }
  }
  if (jmax == -1)
    fatal_error(0,"Atom %s not found in database in residue %s",
		searchname,rtp->resname);
  if (kmax != strlen(searchname))
    fatal_error(0,"Atom %s not found in database in residue %s, "
		"it looks a bit like %s",
		searchname,rtp->resname,*(rtp->atomname[jmax]));
  return jmax;
}

