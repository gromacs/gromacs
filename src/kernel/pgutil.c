/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_pgutil_c = "$Id$";

#include "pgutil.h"
#include "string.h"
	
atom_id search_atom(char *type,int start,int natoms,t_atom at[],char **anm[])
{
  int     i,resnr=-1;
  bool    bPrevious,bNext;

  bPrevious = (strchr(type,'-') != NULL);
  bNext     = (strchr(type,'+') != NULL);

  if (!bPrevious) {
    resnr = at[start].resnr;
    if (bNext) {
      /* The next residue */
      type++;
      while ((start<natoms) && (at[start].resnr == resnr))
	start++;
      if (start < natoms)
	resnr = at[start].resnr;
    }
    for(i=start; (i<natoms) && (bNext || (at[i].resnr == resnr)); i++)
      if (strcasecmp(type,*(anm[i]))==0)
	return (atom_id) i;
  }
  else {
    /* The previous residue */
    type++;
    if (start > 0)
      resnr = at[start-1].resnr;
    for(i=start-1; (i>=0) /*&& (at[i].resnr == resnr)*/; i--)
      if (strcasecmp(type,*(anm[i]))==0)
	return (atom_id) i;
  }
  return NO_ATID;
}

void set_at(t_atom *at,real m,real q,int type,int resnr)
{
  at->m=m;
  at->q=q;
  at->type=type;
  at->resnr=resnr;
}

