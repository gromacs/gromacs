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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _struc2_h
#define _struc2_h

static char *SRCID_struc2_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) struc2.h 1.2 15 Sep 1993"
#endif /* HAVE_IDENT */
#include "typedefs.h"

typedef struct {
  atom_id d;
  atom_id a;
  atom_id h;
} t_dah;

typedef struct {
  real  a_dist; 
  real  h_dist;
  real  d_dist;
  int   nr;	
} t_water;		
	
	
typedef struct {
  real    distance_ad;
  real    distance_ah;
  real    angle;
  bool    h_bond;
  t_water water;
} t_info;

typedef struct {
  t_info *info;
  atom_id  a;
  atom_id  d;
  atom_id  h;
  int    monitor;
} t_hbond;

typedef struct {
  int     nrt;
  int     nr;
  t_hbond *hbond;
  real    *time;
} t_list;

extern void determine_2struc(FILE *out,t_topology *top,t_list *list);

#endif	/* _struc2_h */
