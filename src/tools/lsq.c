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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_lsq_c = "$Id$";

#include "typedefs.h"
#include "gstat.h"
#include "vec.h"

void init_lsq(t_lsq *lsq)
{
  lsq->yy=lsq->yx=lsq->xx=lsq->sx=lsq->sy=0.0;
  lsq->np=0;
}

int npoints_lsq(t_lsq *lsq)
{
  return lsq->np;
}

void done_lsq(t_lsq *lsq)
{
  init_lsq(lsq);
}

void add_lsq(t_lsq *lsq,real x,real y)
{
  lsq->yy+=y*y;
  lsq->yx+=y*x;
  lsq->xx+=x*x;
  lsq->sx+=x;
  lsq->sy+=y;
  lsq->np++;
}

void get_lsq_ab(t_lsq *lsq,real *a,real *b)
{
  double yx,xx,sx,sy;
  
  yx=lsq->yx/lsq->np;
  xx=lsq->xx/lsq->np;
  sx=lsq->sx/lsq->np;
  sy=lsq->sy/lsq->np;
  
  (*a)=(yx-sx*sy)/(xx-sx*sx);
  (*b)=(sy)-(*a)*(sx);
}

real aver_lsq(t_lsq *lsq)
{
  if (lsq->np == 0)
    fatal_error(0,"No points in distribution\n");
  
  return (lsq->sy/lsq->np);
}

real sigma_lsq(t_lsq *lsq)
{
  if (lsq->np == 0)
    fatal_error(0,"No points in distribution\n");
    
  return sqrt(lsq->yy/lsq->np - sqr(lsq->sy/lsq->np));
}

real error_lsq(t_lsq *lsq)
{
  return sigma_lsq(lsq)/sqrt(lsq->np);
}

