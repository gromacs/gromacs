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

#ifndef _pull_h
#define _pull_h

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#include "vec.h"
#include "typedefs.h"

/* This file contains datatypes and function declarations necessary 
   for mdrun to interface with the pull code */

/* Get the distance to the reference and deviation for pull group g */
extern void get_pullgrp_distance(t_pull *pull,int g,matrix box,double t,
				 dvec dr,dvec dev);

/* Determine the COM pull forces and add them to f, return the potential */
extern real pull_potential(int ePull, t_pull *pull,
			   rvec *x, rvec *f, tensor vir, 
			   matrix box, t_topology *top, double t,
			   t_mdatoms *md, t_commrec *cr);

/* Constrain the coordinates xp in the directions in x
 * and also constrain v when v!=NULL.
 */
extern void pull_constraint(t_pull *pull, rvec *x, rvec *xp, rvec *v,
			    tensor vir, matrix box, t_topology *top,
			    real dt, double t,
			    t_mdatoms *md, t_commrec *cr);

/* Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 */
extern void dd_make_local_pull_groups(gmx_domdec_t *dd,
				      t_pull *pull,t_mdatoms *md);

/* get memory and initialize the fields of pull that still need it, and
   do runtype specific initialization */
extern void init_pull(FILE *log,  
                      t_inputrec *ir, /* the inputrec */
                      int nfile,       
                      t_filenm fnm[], /* standard filename struct */
                      rvec *x,        /* all coordinates */
                      t_atoms *atoms, /* masses of all atoms */
                      matrix box,     
                      t_commrec * cr, /* struct for communication info */
		      int start,      /* for particle decomposition */ 
		      int end
                      );

/* Print the pull output (x and/or f) */
extern void pull_print_output(t_pull *pull, int step, real time);

/* In pullutil.c */

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

#endif
