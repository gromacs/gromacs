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
static char *SRCID_pgutil_c = "$Id$";

#include "pgutil.h"
#include "string.h"
	
int search_atom(char *type,int start,int natoms,t_atom at[],char **anm[])
{
  int  i,resnr;
  bool bPrevious,bNext;

  bPrevious = (strchr(type,'-') != NULL);

  if (!bPrevious) {
    resnr = at[start].resnr;
    if (strchr(type,'+') != NULL) {
      /* The next residue */
      type++;
      while ((start<natoms) && (at[start].resnr == resnr))
	start++;
      if (start < natoms)
	resnr = at[start].resnr;
    }
    for(i=start; (i<natoms) && (at[i].resnr == resnr); i++)
      if (strcasecmp(type,*(anm[i]))==0)
	return i;
  }
  else {
    /* The previous residue */
    type++;
    if (start > 0)
      resnr = at[start-1].resnr;
    for(i=start-1; (i>=0) && (at[i].resnr == resnr); i--)
      if (strcasecmp(type,*(anm[i]))==0)
	return i;
  }
  return -1;
}

void set_at(t_atom *at,real m,real q,int type,int resnr)
{
  at->m=m;
  at->q=q;
  at->type=type;
  at->resnr=resnr;
}

