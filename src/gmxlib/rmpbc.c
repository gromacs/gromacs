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
 * Grunge ROck MAChoS
 */
static char *SRCID_rmpbc_c = "$Id$";

#include "sysstuff.h"
#include "typedefs.h"
#include "mshift.h"
#include "pbc.h"
#include "gstat.h"
#include "futil.h"
	
void rm_pbc(t_idef *idef,int natoms,matrix box,rvec x[],rvec x_s[])
{
  static bool    bFirst=TRUE;
  static t_graph *graph;
  rvec   sv[SHIFTS],box_size;
  
  if (bFirst) {
    graph=mk_graph(idef,natoms,FALSE);
    bFirst=FALSE;
  }
  mk_mshift(stdout,graph,box,x);
  calc_shifts(box,box_size,sv,FALSE);
  shift_x(graph,sv,x,x_s);
}

