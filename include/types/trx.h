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
 * Good ROcking Metal Altar for Chronical Sinners
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

