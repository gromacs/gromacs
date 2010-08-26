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

#ifndef _ns_h
#define _ns_h

#include <stdio.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "pbc.h"
#include "tgroup.h"
#include "network.h"


#ifdef __cplusplus
extern "C" {
#endif

/****************************************************
 *
 *    U T I L I T I E S May be found in ns.c
 *
 ****************************************************/

void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr);
/* 
 * nn is the number of energy terms in the energy matrix
 * (ngener*(ngener-1))/2
 * start is the first atom on this processor
 * homenr is the number of atoms on this processor
 */
 
int calc_naaj(int icg,int cgtot);
/* Calculate the number of charge groups to interact with for icg */

/****************************************************
 *
 *    N E I G H B O R  S E A R C H I N G
 *
 *    Calls either ns5_core (when grid selected in .mdp file)
 *    or ns_simple_core (when simple selected in .mdp file)
 *
 *    Return total number of pairs searched 
 *
 ****************************************************/
void init_ns(FILE *fplog,const t_commrec *cr,
		    gmx_ns_t *ns,t_forcerec *fr,
		    const gmx_mtop_t *mtop,
		    matrix box);

int search_neighbours(FILE *log,t_forcerec *fr,
			     rvec x[],matrix box,
			     gmx_localtop_t *top,
			     gmx_groups_t *groups,
			     t_commrec *cr,
			     t_nrnb *nrnb,t_mdatoms *md,
			     real lambda,real *dvdlambda,
			     gmx_grppairener_t *grppener,
			     gmx_bool bFillGrid,
			     gmx_bool bDoLongRange,
			     gmx_bool bDoForces,rvec *f);
 

/* Debugging routines from wnblist.c */
void dump_nblist(FILE *out,t_commrec *cr,t_forcerec *fr,int nDNL);

int read_nblist(FILE *in,FILE *out,int **mat,int natoms,gmx_bool bSymm);
/* Returns total number of neighbors. If bSymm the matrix is symmetrized. */

int natoms_beyond_ns_buffer(t_inputrec *ir,t_forcerec *fr,t_block *cgs,
				   matrix scale_tot,rvec *x);
/* Returns the number of atoms that moved beyond the ns buffer */

#ifdef __cplusplus
}
#endif


#endif	/* _ns_h */
