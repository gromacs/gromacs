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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_pgutil_c = "$Id$";

#include "pgutil.h"
#include "string.h"
	
int search_atom(char *type,int start,int natoms,char **atom[])
{
  int  i;
  bool bForwards;

  bForwards=(strchr(type,'-') == NULL);

  if (bForwards) {
    for(i=start; (i<natoms); i++)
      if (strcasecmp(type,*(atom[i]))==0)
	return i;
  }
  else {
    type++;
    for(i=start-1; (i>=0); i--)
      if (strcasecmp(type,*(atom[i]))==0)
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

