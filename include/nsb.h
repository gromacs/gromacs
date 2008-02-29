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

#ifndef _nsb_h
#define _nsb_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern void calc_nsbshift(FILE *fp,t_nsborder *nsb,t_idef *idef);
/* Calculates the shift and bshift variables */

extern void calc_nsb(FILE *fp,t_block *cgs,int nprocs,
		     t_nsborder *nsb,int nstDlb,t_idef *idef);
/* Calculate which blocks of charge groups should be calculated,
 * depending on processor number.
 */

extern void print_nsb(FILE *fp,char *title,t_nsborder *nsb);
/* Print the values to file */

/*extern bool cg_avail(int icg,int jcg,t_nsborder *nsb);*/
/* Determine whether a charge group jcg is in the domain of icg:
 * this is necessary for searching on a grid, to avoid double
 * pairs of interactions.
 */

/*extern int  cg_index(int icg,t_nsborder *nsb);*/
/* Perform a modulo calculation giving the correct cg index */

#endif
