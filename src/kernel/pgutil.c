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
/* This file is completely threadsafe - keep it that way! */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "string2.h"
#include "pgutil.h"
#include "string.h"
#include "gmx_fatal.h"

static void atom_not_found(int fatal_errno,const char *file,int line,
			   const char *atomname,int resind,
			   const char *bondtype,gmx_bool bDontQuit)
{
  if (strcmp(bondtype,"check") != 0) {
    if (bDontQuit) {
      gmx_warning("WARNING: Atom %s is used in an interaction of type %s in the\n"
		  "topology database, but an atom of that name was not found in\n"
		  "residue number %d.\n",
		  atomname,bondtype,resind+1);
    } else {
      gmx_fatal(fatal_errno,file,line,
		"Atom %s is used in an interaction of type %s in the topology\n"
		"database, but an atom of that name was not found in residue\n"
		"number %d.\n",
		atomname,bondtype,resind+1);
    }
  }
}
	
atom_id search_atom(const char *type,int start,int natoms,t_atom at[],
		    char ** const * anm,
		    const char *bondtype,gmx_bool bDontQuit)
{
  int     i,resind=-1;
  gmx_bool    bPrevious,bNext;

  bPrevious = (strchr(type,'-') != NULL);
  bNext     = (strchr(type,'+') != NULL);

  if (!bPrevious) {
    resind = at[start].resind;
    if (bNext) {
      /* The next residue */
      type++;
      while ((start<natoms) && (at[start].resind == resind))
	start++;
      if (start < natoms)
	resind = at[start].resind;
    }
    
    for(i=start; (i<natoms) && (bNext || (at[i].resind == resind)); i++) {
      if (anm[i] && gmx_strcasecmp(type,*(anm[i]))==0)
	return (atom_id) i;
    }
    if (!(bNext && at[start].resind==at[natoms-1].resind))
      atom_not_found(FARGS,type,at[start].resind,bondtype,bDontQuit);
  }
  else {
    /* The previous residue */
    type++;
    if (start > 0)
      resind = at[start-1].resind;
    for(i=start-1; (i>=0) /*&& (at[i].resind == resind)*/; i--)
      if (gmx_strcasecmp(type,*(anm[i]))==0)
	return (atom_id) i;
    if (start > 0)
      atom_not_found(FARGS,type,at[start].resind,bondtype,bDontQuit);
  }
  return NO_ATID;
}

void set_at(t_atom *at,real m,real q,int type,int resind)
{
  at->m=m;
  at->q=q;
  at->type=type;
  at->resind=resind;
}
