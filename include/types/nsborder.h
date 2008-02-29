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

typedef struct {
  int  nodeid;			/* Node id	         		*/
  int  nnodes;			/* The number of nodes    		*/
  int  cgtotal; 		/* Total number of charge groups	*/
  int  natoms;			/* Total number of atoms		*/
  int  nstDlb;                  /* Every how many steps must we do load */
                                /* balancing                            */
  int  shift,bshift;		/* Coordinates are shifted left for     */
                                /* 'shift' systolic pulses, and right   */
				/* for 'bshift' pulses. Forces are      */
				/* shifted right for 'shift' pulses     */
				/* and left for 'bshift' pulses         */
				/* This way is not necessary to shift   */
				/* the coordinates over the entire ring */
  int  homenr[MAXNODES];       	/* The number of home particles		*/
  int  index[MAXNODES];		/* The starting of the home atoms	*/
  int  cgload[MAXNODES];        /* Division of charge groups over CPUS  */
                                /* This is static, i.e. it does not     */
				/* change during the simulation         */
  int  workload[MAXNODES];      /* This is the load for neighbor-       */
                                /* searching, this is initially the same*/
				/* as cgload, but may change due to     */
				/* dynamic load balancing               */
} t_nsborder;

#define START(nsb)  ((nsb)->index[(nsb)->nodeid])
#define HOMENR(nsb) ((nsb)->homenr[(nsb)->nodeid])
#define CG0(nsb)    (((nsb)->nodeid == 0) ? 0 : (nsb)->cgload[(nsb)->nodeid-1])
#define CG1(nsb)    ((nsb)->cgload[(nsb)->nodeid])
