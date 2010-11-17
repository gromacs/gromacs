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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "invblock.h"
#include "gmx_fatal.h"

atom_id *make_invblock(const t_block *block,int nr)
{
  int i,j;
  atom_id *invblock;
  
  snew(invblock,nr+1);
  /* Mark unused numbers */
  for (i=0; i<=nr; i++) 
    invblock[i]=NO_ATID; 
  for (i=0; (i<block->nr); i++)
    for (j=block->index[i]; (j<block->index[i+1]); j++) 
      if (invblock[j] == NO_ATID)
	invblock[j]=i;
      else
	gmx_fatal(FARGS,"Double entries in block structure. Item %d is in blocks %d and %d\n"
		  " Cannot make an unambiguous inverse block.",
		  j,i,invblock[j]);
  return invblock;
}

atom_id *make_invblocka(const t_blocka *block,int nr)
{
  int i,j;
  atom_id *invblock;
  
  snew(invblock,nr+1);
  /* Mark unused numbers */
  for (i=0; i<=nr; i++) 
    invblock[i]=NO_ATID; 
  for (i=0; (i<block->nr); i++)
    for (j=block->index[i]; (j<block->index[i+1]); j++) 
      if (invblock[block->a[j]] == NO_ATID)
	invblock[block->a[j]]=i;
      else
	gmx_fatal(FARGS,"Double entries in block structure. Item %d is in blocks %d and %d\n"
		  " Cannot make an unambiguous inverse block.",
		  j,i,invblock[block->a[j]]);
  return invblock;
}

