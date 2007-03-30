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
 * GROwing Monsters And Cloning Shrimps
 */

#ifndef _pull_internal_h
#define _pull_internal_h

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#include "vec.h"
#include "typedefs.h"


/* NOTE: external routines are located in include/pull.h */
  
/* For triclinic boxes only correct for distances
 * shorter than half the smallest diagonal box element.
 */
extern void pull_d_pbc_dx(int npbcdim,
			  matrix box,const dvec x1, const dvec x2, dvec dx);

/* Calculates centers of mass all pull groups */
extern void pull_calc_coms(t_commrec *cr,
			   t_pull *pull,   /* the pull group */
			   t_mdatoms *md,  /* all atoms */
			   rvec x[],       /* local coordinates */
			   rvec *xp,       /* updated x, can be NULL */
			   matrix box);    

/* read the parameter file .ppa and write out what was read in */
extern void read_pullparams(t_pull *pull,
                            char *infile,
                            char *outfile); 

/* Print header at top of pdo file */
extern void print_pull_header(FILE *out,t_pull * pull);

#endif
