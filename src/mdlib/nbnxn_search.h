/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustr
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
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

#ifndef _nbnxn_search_h
#define _nsnxn_search_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Returns the j-cluster size for kernel of type nb_kernel_type */
int nbnxn_kernel_to_cj_size(int nb_kernel_type);

/* Tells if the pair-list corresponding to nb_kernel_type is simple.
 * Returns FALSE for super-sub type pair-list.
 */
gmx_bool nbnxn_kernel_pairlist_simple(int nb_kernel_type);

/* Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 */
real nbnxn_get_rlist_effective_inc(int cluster_size,real atom_density);

/* Allocates and initializes a pair search data structure */
void nbnxn_init_search(nbnxn_search_t * nbs_ptr,
                       ivec *n_dd_cells,
                       gmx_domdec_zones_t *zones,
                       int nthread_max);

/* Put the atoms on the pair search grid.
 * Only atoms a0 to a1 in x are put on the grid.
 * The atom_density is used to determine the grid size.
 * When atom_density=-1, the density is determined from a1-a0 and the corners.
 * With domain decomposition part of the n particles might have migrated,
 * but have not been removed yet. This count is given by nmoved.
 * When move[i] < 0 particle i has migrated and will not be put on the grid.
 * Without domain decomposition move will be NULL.
 */
void nbnxn_put_on_grid(nbnxn_search_t nbs,
                       int ePBC,matrix box,
                       int dd_zone,
                       rvec corner0,rvec corner1,
                       int a0,int a1,
                       real atom_density,
                       const int *atinfo,
                       rvec *x,
                       int nmoved,int *move,
                       int nb_kernel_type,
                       nbnxn_atomdata_t *nbat);

/* As nbnxn_put_on_grid, but for the non-local atoms
 * with domain decomposition. Should be called after calling
 * nbnxn_search_put_on_grid for the local atoms / home zone.
 */
void nbnxn_put_on_grid_nonlocal(nbnxn_search_t nbs,
                                const gmx_domdec_zones_t *zones,
                                const int *atinfo,
                                rvec *x,
				int nb_kernel_type,
                                nbnxn_atomdata_t *nbat);

/* Add simple grid type information to the local super/sub grid */
void nbnxn_grid_add_simple(nbnxn_search_t nbs,
			   nbnxn_atomdata_t *nbat);

/* Return the number of x and y cells in the local grid */
void nbnxn_get_ncells(nbnxn_search_t nbs,int *ncx,int *ncy);

/* Return the order indices *a of the atoms on the ns grid, size n */
void nbnxn_get_atomorder(nbnxn_search_t nbs,int **a,int *n);

/* Renumber the atom indices on the grid to consecutive order */
void nbnxn_set_atomorder(nbnxn_search_t nbs);

/* Initializes a set of pair lists stored in nbnxn_pairlist_set_t */
void nbnxn_init_pairlist_set(nbnxn_pairlist_set_t *nbl_list,
                             gmx_bool simple, gmx_bool combined,
                             nbnxn_alloc_t *alloc,
                             nbnxn_free_t  *free);

/* Make a apir-list with radius rlist, store it in nbl.
 * The parameter min_ci_balanced sets the minimum required
 * number or roughly equally sized ci blocks in nbl.
 * When set >0 ci lists will be chopped up when the estimate
 * for the number of equally sized lists is below min_ci_balanced.
 */
void nbnxn_make_pairlist(const nbnxn_search_t nbs,
			 const nbnxn_atomdata_t *nbat,
			 const t_blocka *excl,
			 real rlist,
			 int min_ci_balanced,
			 nbnxn_pairlist_set_t *nbl_list,
			 int iloc,
			 int nb_kernel_type,
			 t_nrnb *nrnb);

#ifdef __cplusplus
}
#endif

#endif
