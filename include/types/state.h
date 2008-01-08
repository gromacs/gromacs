/*
 * $Id$
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 * The t_state struct should contain all the (possibly) non-static
 * information required to define the state of the system.
 * Currently the random seeds for SD and BD are missing.
 */

#define STATE_HAS_X   (1<<0)
#define STATE_HAS_V   (1<<1)
#define STATE_HAS_SDX (1<<2)
#define STATE_HAS_CGP (1<<3)

typedef struct
{
  int           natoms;
  int           ngtc;
  int           flags;  /* Flags telling which entries are present      */
  real          lambda; /* the free energy switching parameter          */
  matrix 	box;    /* box vector coordinates                      	*/
  matrix 	boxv;   /* box velocitites for Parrinello-Rahman pcoupl */
  matrix        pcoupl_mu; /* for Berendsen pcoupl                      */
  real          *nosehoover_xi; /* for Nose-Hoover tcoupl (ngtc)        */
  int           nalloc; /* Allocation size for x, v and sd_x when !=NULL*/
  rvec          *x;     /* the coordinates (natoms)                     */
  rvec          *v;     /* the velocities (natoms)                      */
  rvec          *sd_X;  /* random part of the x update for stoch. dyn.  */
  rvec          *cg_p;  /* p vector for conjugate gradient minimization */
  int           ddp_count; /* The DD partitioning count for this state  */
  int           ddp_count_cg_gl; /* The DD part. count for index_gl     */
  int           ncg_gl; /* The number of local charge groups            */
  int           *cg_gl; /* The global cg number of the local cgs        */
  int           cg_gl_nalloc; /* Allocation size of cg_gl;              */
} t_state;
