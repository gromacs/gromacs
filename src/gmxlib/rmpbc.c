/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_rmpbc_c = "$Id$";
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "mshift.h"
#include "pbc.h"
#include "gstat.h"
#include "futil.h"
#include "vec.h"	

void rm_pbc(t_idef *idef,int natoms,matrix box,rvec x[],rvec x_s[])
{
  typedef struct {
    int     natoms;
    t_graph *gr;
  } multi_graph;
  
  static int ngraph=0;
  static multi_graph *mgraph=NULL;
  static bool bFirst=TRUE;
  rvec   sv[SHIFTS],box_size;
  int    n,i;
  bool   bNeedToCopy;

  bNeedToCopy = (x != x_s);

  if (box[0][0]) {
    if (idef->ntypes!=-1) {
      n=-1;
      for(i=0; i<ngraph; i++)
	if (mgraph[i].natoms==natoms)
	  n=i;
      if (n==-1) {
	/* make a new graph if there isn't one with this number of atoms */
	n=ngraph;
	ngraph++;
	srenew(mgraph,ngraph);
	mgraph[n].natoms=natoms;
	mgraph[n].gr=mk_graph(idef,natoms,FALSE,FALSE);
      }
      mk_mshift(stdout,mgraph[n].gr,box,x);
      calc_shifts(box,box_size,sv);
      shift_x(mgraph[n].gr,box,x,x_s);
      bNeedToCopy=FALSE;
    } else if (bFirst) {
      fprintf(stderr,
	      "\nWarning: can not make broken molecules whole without a run input file,\n         don't worry, mdrun doesn't write broken molecules\n\n");
      bFirst=FALSE;
    }
  }
  if (bNeedToCopy)
    for (i=0; i<natoms; i++)
      copy_rvec(x[i],x_s[i]);
}

