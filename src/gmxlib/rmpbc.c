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
  
  static int ngraph;
  static multi_graph *mgraph=NULL;
  static t_graph *graph;
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
	mgraph[n].gr=mk_graph(idef,natoms,FALSE);
      }
      mk_mshift(stdout,mgraph[n].gr,box,x);
      calc_shifts(box,box_size,sv,FALSE);
      shift_x(mgraph[n].gr,sv,x,x_s);
      bNeedToCopy=FALSE;
    } else
      if (bFirst) {
	fprintf(stderr,
		"\nWarning: can not make broken molecules whole without a run input file\n");
	bFirst=FALSE;
      }
  }
  if (bNeedToCopy)
    for (i=0; i<natoms; i++)
      copy_rvec(x[i],x_s[i]);
}

