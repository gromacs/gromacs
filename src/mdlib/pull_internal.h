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


/* NOTE: external routines and datatypes are located in
   include/pull.h
*/
#include "pull.h"

#define NINT(x) (x<0?((int)((x)-0.5)):((int)((x)+0.5)))
#define DEBUG_START fprintf(stderr,"\n\nDEBUG\n");
#define DEBUG_END   fprintf(stderr,"\nEND DEBUG\n\n");
  
/* Not correct for all trilinic boxes !!!
   Should be as pbc_dx and moved to pbc.c
*/
extern void d_pbc_dx(matrix box,const dvec x1, const dvec x2, dvec dx);

extern void put_dvec_in_box(matrix box,dvec v);

/* print to output file for the various types of runs */
extern void print_umbrella(t_pull *pull, int step, real t);
extern void print_afm(t_pull *pull, int step, real t);
extern void print_constraint(t_pull *pull, int step, real t);


/* calculate center of mass of index group, making sure it's inside the box,
   stores it in pg->x_unc. */
extern void calc_com(t_pullgrp *pg,  /* the pull group */
		     rvec x[],       /* coordinates of all atoms in system */ 
                     t_mdatoms *md,  /* all atoms */
                     matrix box);    


/* calculate center of mass of all atoms x[], index needed to get the right
   masses from the atom array,
   stores it in pg->x_unc. */
extern void calc_com2(t_pullgrp *pg,  /* the pull group */
		      rvec x[],       /* coordinates to calc. com from */
                      t_mdatoms *md,  /* all atoms */
                      matrix box); 


/* calculate a running average for center of mass */
extern void calc_running_com(t_pull *pull);


/* calculate the center of mass from the true coordinates, without
   corrections for pbc */
extern void correct_t0_pbc(t_pull *pull, 
                           rvec x[], 
                           t_mdatoms *md,
                           matrix box);


/* read the parameter file .ppa and write out what was read in */
extern void read_pullparams(t_pull *pull,
                            char *infile,
                            char *outfile); 

/* find all atoms in group pull->idx[pull->n] that are inside a cylinder
   with as origin com[i][x],com[i][y] with radius pull->r and possibly
   a switch function pull->rc. Remember their weight. Now each group i
   has its own reference group (HOW?) with com defined as 
   Sum(wi*mi*ri)/(Sum(wi*mi). Basically, build an index structure with
   the reference groups for the groups i, plus an array with the 
   weight factors for each of the atoms in those index groups? 
   */
extern void make_refgrps(t_pull *pull,
                         matrix box,
                         t_mdatoms *md);


/* write a numbered .gro file in procedure to make starting structures */
extern void dump_conf(t_pull *pull,
                      rvec x[],        /* all coordinates */
                      matrix box,      /* box */
                      t_topology *top, /* names and residue info */
                      int nout,        /* sequence number of this file */
                      real time);      /* time in simulation */

/* Print header at top of pdo file */
extern void print_pull_header(t_pull * pull);

#endif
