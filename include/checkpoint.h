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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _checkpoint_h
#define _checkpoint_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "gmxfio.h"

/* Write a checkpoint to fn */
extern void write_checkpoint(char *fn,FILE *fplog,t_commrec *cr,
			     int eIntegrator,int simulation_part,int step,double t,
			     t_state *state);

/* Loads a checkpoint from fn for run continuation.
 * Generates a fatal error on system size mismatch.
 * The master node reads the file
 * and communicates all the modified number of steps and the parallel setup,
 * but not the state itself.
 */
extern void load_checkpoint(char *fn,FILE *fplog,
			    t_commrec *cr,bool bPartDecomp,ivec dd_nc,
			    t_inputrec *ir,t_state *state,bool *bReadRNG, 
			    bool *bReadEkin,
			    bool bTruncateOutputFiles);

/* Read the state from checkpoint file.
 * Arrays in state that are NULL are allocated.
 * If bReadRNG=TRUE a RNG state compatible with the current
 * number of nodes was read.
 */
extern void read_checkpoint_state(char *fn,int *simulation_part,int *step,double *t,t_state *state);

/* Read everything that can be stored in t_trxframe from a checkpoint file */
extern void read_checkpoint_trxframe(int fp,t_trxframe *fr);

/* Print the complete contents of checkpoint file fn to out */
extern void list_checkpoint(char *fn,FILE *out);

/* Read just the simulation 'generation'. This is necessary already at the beginning of mdrun,
 * to be able to rename the logfile correctly.
 */
int
read_checkpoint_simulation_part(char *filename);


#endif
