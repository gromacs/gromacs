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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"

double compute_io(t_inputrec *ir,int natoms,int nrener,int nrepl)
{
  int nsteps = ir->nsteps-ir->init_step;
  int nstx=0,nstv=0,nstf=0,nste=0,nstlog=0,nstxtc=0,nfep=0;
  double cio;
  
  if (ir->nstxout > 0)
    nstx = 1 + nsteps / ir->nstxout;
  if (ir->nstvout > 0)
    nstv = 1 + nsteps / ir->nstvout;
  if (ir->nstfout > 0)
    nstf = 1 + nsteps / ir->nstfout;
  if (ir->nstxtcout > 0)
    nstxtc = 1 + nsteps / ir->nstxtcout;
  if (ir->nstlog > 0)
    nstlog = 1 + nsteps / ir->nstlog;
  if (ir->nstenergy > 0)
    nste = 3 + nsteps % ir->nstenergy;
  cio  = 80*natoms;
  cio += (nstx+nstf+nstv)*sizeof(real)*3*natoms;
  cio += nstxtc*natoms*5; /* roughly 5 bytes per atom */
  cio += nstlog*nrener*16*2; /* 16 bytes per energy term plus header */
  cio += nste*nrener*sizeof(t_energy);
  
  return 1e-6*cio*nrepl;
}
