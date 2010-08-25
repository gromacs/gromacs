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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _pull_h
#define _pull_h

#include "vec.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


/* This file contains datatypes and function declarations necessary 
   for mdrun to interface with the pull code */

/* Get the distance to the reference and deviation for pull group g */
void get_pullgrp_distance(t_pull *pull,t_pbc *pbc,int g,double t,
				 dvec dr,dvec dev);

/* Set the all the pull forces to zero */
void clear_pull_forces(t_pull *pull);

/* Determine the COM pull forces and add them to f, return the potential */
real pull_potential(int ePull,t_pull *pull, t_mdatoms *md, t_pbc *pbc,
			   t_commrec *cr, double t, real lambda,
			   rvec *x, rvec *f, tensor vir, real *dvdlambda);

/* Constrain the coordinates xp in the directions in x
 * and also constrain v when v!=NULL.
 */
void pull_constraint(t_pull *pull, t_mdatoms *md, t_pbc *pbc,
			    t_commrec *cr, double dt, double t,
			    rvec *x, rvec *xp, rvec *v, tensor vir);

/* Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 */
void dd_make_local_pull_groups(gmx_domdec_t *dd,
				      t_pull *pull,t_mdatoms *md);

/* get memory and initialize the fields of pull that still need it, and
   do runtype specific initialization */
void init_pull(FILE *fplog,  
                      t_inputrec *ir, /* the inputrec */
                      int nfile,       
                      const t_filenm fnm[], /* standard filename struct */
                      gmx_mtop_t *mtop, /* the topology of the whole system */
                      t_commrec * cr, /* struct for communication info */
                      const output_env_t oenv,  /* output options */
                      gmx_bool bOutFile,   /* open output files */
                      unsigned long Flags);

/* Close the pull output files */
void finish_pull(FILE *fplog,t_pull *pull);

/* Print the pull output (x and/or f) */
void pull_print_output(t_pull *pull, gmx_large_int_t step, double time);

/* In pullutil.c */

/* Calculates centers of mass all pull groups */
void pull_calc_coms(t_commrec *cr,
			   t_pull *pull,   /* the pull group */
			   t_mdatoms *md,  /* all atoms */
			   t_pbc *pbc,
			   double t,       /* only used for cylinder ref. */
			   rvec x[],       /* local coordinates */
			   rvec *xp        /* updated x, can be NULL */
			   );    

#ifdef __cplusplus
}
#endif

#endif
