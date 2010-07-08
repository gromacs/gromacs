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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "mshift.h"
#include "pbc.h"
#include "gstat.h"
#include "futil.h"
#include "vec.h"	

typedef struct gmx_rmpbc {
  int     natoms;
  int     ePBC;
  t_graph *gr;
} koeiepoep;

gmx_rmpbc_t gmx_rmpbc_init(t_idef *idef,int ePBC,int natoms,
			   matrix box)
{
  gmx_rmpbc_t gpbc;
  
  snew(gpbc,1);
  
  gpbc->natoms=natoms;

  if (ePBC == -1)
    gpbc->ePBC = guess_ePBC(box);
  else
    gpbc->ePBC = ePBC;
    
  if (idef->ntypes <= 0)
    fprintf(stderr,
	    "\nWarning: if there are broken molecules in the trajectory file,\n"
	    "         they can not be made whole without a run input file\n\n");
  else if (ePBC == epbcNONE) 
    fprintf(stderr,"\nNot treating periodicity since it is turned off in the input file\n");
  else
    gpbc->gr = mk_graph(NULL,idef,0,natoms,FALSE,FALSE);

  return gpbc;
}

void gmx_rmpbc_done(gmx_rmpbc_t gpbc)
{
  if (NULL != gpbc->gr)
    done_graph(gpbc->gr);
}

void gmx_rmpbc(gmx_rmpbc_t gpbc,matrix box,rvec x[],rvec x_s[])
{
  int    i;

  if (NULL != gpbc->gr) {
    mk_mshift(stdout,gpbc->gr,gpbc->ePBC,box,x);
    shift_x(gpbc->gr,box,x,x_s);
  }
  if (x != x_s)
    for (i=0; i<gpbc->natoms; i++)
      copy_rvec(x[i],x_s[i]);
}

void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
  real dist;
  int  n,m,d;
  
  /* check periodic boundary */
  for(n=1;(n<atoms->nr);n++) {
    for(m=DIM-1; m>=0; m--) {
      dist = x[n][m]-x[n-1][m];
      if (fabs(dist) > 0.9*box[m][m]) { 
	if ( dist >  0 )
	  for(d=0; d<=m; d++)
	    x[n][d] -= box[m][d];
	else
	  for(d=0; d<=m; d++)
	    x[n][d] += box[m][d];
      } 	
    }
  }
}

