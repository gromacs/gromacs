/*
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "sysstuff.h"
#include "macros.h"
#include "symtab.h"
#include "string2.h"
#include "gpp_bond_atomtype.h"

typedef struct {
  int           nr;		/* The number of atomtypes		*/
  char          ***atomname;	/* Names of the atomtypes		*/
} gpp_bond_atomtype;

int get_bond_atomtype_type(char *str,t_bond_atomtype at)
{
  gpp_bond_atomtype *ga = (gpp_bond_atomtype *) at;
  
  int i;

  for (i=0; (i<ga->nr); i++)
    if (gmx_strcasecmp(str,*(ga->atomname[i])) == 0)
      return i;
  
  return NOTSET;
}

char *get_bond_atomtype_name(int nt, t_bond_atomtype at)
{
  gpp_bond_atomtype *ga = (gpp_bond_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NULL;
  
  return *(ga->atomname[nt]);
}

t_bond_atomtype init_bond_atomtype(void)
{
  gpp_bond_atomtype *ga;
  
  snew(ga,1);
  
  return (t_bond_atomtype ) ga;
}

void add_bond_atomtype(t_bond_atomtype at,t_symtab *tab,
		       char *name)
{
  gpp_bond_atomtype *ga = (gpp_bond_atomtype *) at;
  
  ga->nr++;
  srenew(ga->atomname,ga->nr); 
  ga->atomname[ga->nr-1] = put_symtab(tab,name);
}

void done_bond_atomtype(t_bond_atomtype *at)
{
  gpp_bond_atomtype *ga = (gpp_bond_atomtype *) *at;

  sfree(ga->atomname);
  ga->nr = 0;
  sfree(ga);
  
  *at = NULL;
}
