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

#include "sysstuff.h"
#include "typedefs.h"
#include "splittop.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "main.h"
#include "vsite.h"
#include "mdrun.h"

static void split_ilist(FILE *log,t_ilist *il,t_commrec *cr)
{
  t_iatom *ia;
  int     i,start,end,nr;
  
  if (cr->nodeid == 0)
    start=0;
  else
    start=il->multinr[cr->nodeid-1];
  end=il->multinr[cr->nodeid];
  
  nr=end-start;
  if (nr < 0)
    gmx_fatal(FARGS,"Negative number of atoms (%d) on node %d\n"
		"You have probably not used the same value for -np with grompp"
		" and mdrun",
		nr,cr->nodeid);
  snew(ia,nr);

  for(i=0; (i<nr); i++)
    ia[i]=il->iatoms[start+i];

  sfree(il->iatoms);
  il->iatoms=ia;
  
  for(i=0; (i<MAXNODES); i++)
    il->multinr[i]=nr;
  il->nr=nr;
}

static void split_idef(FILE *log,t_idef *idef,t_commrec *cr)
{
  int i;
  
  for(i=0; (i<F_NRE); i++)
    split_ilist(log,&idef->il[i],cr);
}
	
void mdsplit_top(FILE *log,t_topology *top,t_commrec *cr,
		 t_nsborder *nsb, bool *bParallelVsites,
		 t_comm_vsites *vsitecomm)
{
  if (cr->nnodes - cr->npmenodes < 2)
    return;

  *bParallelVsites=setup_parallel_vsites(&(top->idef),cr,nsb,vsitecomm);
  
  split_idef(log,&top->idef,cr);
#ifdef DEBUG
  pr_idef(log,0,"After Split",&(top->idef));
#endif
}
