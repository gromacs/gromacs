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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* The bools indicate whether a field was read from the trajectory.
 * Do not try to use a pointer when its bool is FALSE, as memory might
 * not be allocated.
 */ 
typedef struct
{
  int  flags;     /* flags for read_first/next_frame  */
  int  not_ok;    /* integrity flags (see statutil.h  */
  int  natoms;    /* number of atoms (atoms, x, v, f) */
  real t0;        /* time of the first frame, needed  *
		   * for skipping frames with -dt     */
  bool bTitle;
  char *title;    /* title of the frame               */
  bool bStep;
  int  step;      /* MD step number                   */
  bool bTime;
  real time;      /* time of the frame                */
  bool bLambda;
  real lambda;    /* free energy perturbation lambda  */
  bool bAtoms;
  t_atoms *atoms; /* atoms struct (natoms)            */
  bool bPrec;
  real prec;      /* precision of x, fraction of 1 nm */
  bool bX;
  rvec *x;        /* coordinates (natoms)             */
  bool bV;
  rvec *v;        /* velocities (natoms)              */
  bool bF;
  rvec *f;        /* forces (natoms)                  */
  bool bBox;
  matrix box;     /* the 3 box vectors                */
} t_trxframe;

