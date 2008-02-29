/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"

double compute_io(t_inputrec *ir,t_atoms *atoms,int nrener,int nrepl)
{
  int nsteps = ir->nsteps;
  int natoms = atoms->nr;
  int i,nxtcatoms=0;
  int nstx=0,nstv=0,nstf=0,nste=0,nstlog=0,nstxtc=0,nfep=0;
  double cio;
  
  if (ir->nstxout > 0)
    nstx = 1 + nsteps / ir->nstxout;
  if (ir->nstvout > 0)
    nstv = 1 + nsteps / ir->nstvout;
  if (ir->nstfout > 0)
    nstf = (1 + nsteps) / ir->nstfout;
  if (ir->nstxtcout > 0) {
    for(i=0; i<natoms; i++)
      if (atoms->atom[i].grpnr[egcXTC] == 0)
	nxtcatoms++;
    nstxtc = (1 + nsteps) / ir->nstxtcout;
  }
  if (ir->nstlog > 0)
    nstlog = 1 + nsteps / ir->nstlog;
  if (ir->nstenergy > 0)
    nste = 3 + nsteps % ir->nstenergy;
  cio  = 80*natoms;
  cio += (nstx+nstf+nstv)*sizeof(real)*(3.0*natoms);
  cio += nstxtc*(14*4 + nxtcatoms*5.0); /* roughly 5 bytes per atom */
  cio += nstlog*(nrener*16*2.0); /* 16 bytes per energy term plus header */
  cio += (1.0*nste)*nrener*sizeof(t_energy);
  
  return cio*nrepl/(1024*1024);
}
