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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _nrama_h
#define _nrama_h

static char *SRCID_nrama_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) nrama.h 1.9 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "statutil.h"
#include "mshift.h"

typedef struct {
  bool bShow;
  char *label;
  int  iphi,ipsi; /* point in the dih array of xr... */
} t_phipsi;

typedef struct {
  atom_id ai[4];
  int     mult;
  real    phi0;
  real    ang;
} t_dih;

typedef struct {
  int       ndih;
  t_dih     *dih;
  int       npp;
  t_phipsi  *pp;
  int       traj;
  int       natoms;
  int       amin,amax;
  real      t;
  rvec      *x;
  matrix    box;
  t_idef    *idef;
} t_xrama;

extern void init_rama(char *infile,char *topfile,t_xrama *xr);

extern bool new_data(t_xrama *xr);

#endif	/* _nrama_h */
