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
 * GRoups of Organic Molecules in ACtion for Science
 */

/* The gmx_bools indicate whether a field was read from the trajectory.
 * Do not try to use a pointer when its gmx_bool is FALSE, as memory might
 * not be allocated.
 */ 

#include "atoms.h"
#include "../molfile_plugin.h"
#include "../vmdio.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct trxframe
{
  int  flags;     /* flags for read_first/next_frame  */
  int  not_ok;    /* integrity flags (see statutil.h  */
  gmx_bool bDouble;   /* Double precision?                */
  int  natoms;    /* number of atoms (atoms, x, v, f) */
  real t0;        /* time of the first frame, needed  *
		   * for skipping frames with -dt     */
  real tpf;       /* time of the previous frame, not  */
                  /* the read, but real file frames   */
  real tppf;      /* time of two frames ago           */
                  /* tpf and tppf are needed to       */
                  /* correct rounding errors for -e   */
  gmx_bool bTitle;
  const char *title; /* title of the frame            */
  gmx_bool bStep;
  int  step;      /* MD step number                   */
  gmx_bool bTime;
  real time;      /* time of the frame                */
  gmx_bool bLambda;
  real lambda;    /* free energy perturbation lambda  */
  gmx_bool bAtoms;
  t_atoms *atoms; /* atoms struct (natoms)            */
  gmx_bool bPrec;
  real prec;      /* precision of x, fraction of 1 nm */
  gmx_bool bX;
  rvec *x;        /* coordinates (natoms)             */
  gmx_bool bV;
  rvec *v;        /* velocities (natoms)              */
  gmx_bool bF;
  rvec *f;        /* forces (natoms)                  */
  gmx_bool bBox;
  matrix box;     /* the 3 box vectors                */
  gmx_bool bPBC;
  int  ePBC;      /* the type of pbc                  */
  t_gmxvmdplugin vmdplugin;
} t_trxframe;

#ifdef __cplusplus
}
#endif

