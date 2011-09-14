/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
 */

#ifndef _nsbox_h
#define _nsbox_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Allocates and initializes a neighbor searching data structure */
void gmx_nbsearch_init(gmx_nbsearch_t * nbs_ptr,
                       ivec *n_dd_cells,
                       gmx_domdec_zones_t *zones,
                       gmx_bool simple,
                       int natoms_subcell,
                       int nthread_max);

/* Put the atoms on the neighborsearching grid.
 * Only atoms a0 to a1 in x are put on the grid.
 * The atom_density is used to determine the grid size.
 * When atom_density=-1, the density is determined from a1-a0 and the corners.
 * With domain decomposition part of the n particles might have migrated,
 * but have not been removed yet. This count is given by nmoved.
 * When move[i] < 0 particle i has migrated and will not be put on the grid.
 * Without domain decomposition move will be NULL.
 */
void gmx_nbsearch_put_on_grid(gmx_nbsearch_t nbs,
                              int ePBC,matrix box,
                              int dd_zone,
                              rvec corner0,rvec corner1,
                              int a0,int a1,
                              real atom_density,
                              const int *atinfo,
                              rvec *x,
                              int nmoved,int *move,
                              gmx_nb_atomdata_t *nbat);

/* As gmx_nbsearch_put_on_grid, but for the non-local atoms
 * with domain decomposition. Should be called after calling
 * gmx_nbsearch_put_on_grid for the local atoms / home zone.
 */
void gmx_nbsearch_put_on_grid_nonlocal(gmx_nbsearch_t nbs,
                                       const gmx_domdec_zones_t *zones,
                                       const int *atinfo,
                                       rvec *x,
                                       gmx_nb_atomdata_t *nbat);

/* Return the number of x and y cells in the local grid */
void gmx_nbsearch_get_ncells(gmx_nbsearch_t nbs,int *ncx,int *ncy);

/* Return the order indices *a of the atoms on the ns grid.
 * An index >= *moved indicates and atom that moved to another domain.
 */
void gmx_nbsearch_get_atomorder(gmx_nbsearch_t nbs,int **a,int *moved);

/* Renumber the atom indices on the grid to consecutive order */
void gmx_nbsearch_set_atomorder(gmx_nbsearch_t nbs);

/* Initialize a neighbor list data structure */
void gmx_nblist_init(gmx_nblist_t * nbl,
                     gmx_nbat_alloc_t *alloc,
                     gmx_nbat_free_t  *free);

/* Make a neighborlist with radius rlist, store it in nbl.
 * The parameter min_ci_balanced sets the minimum required
 * number or roughly equally sized ci blocks in nbl.
 * When set >0 ci lists will be chopped up when the estimate
 * for the number of equally sized lists is below min_ci_balanced.
 */
void gmx_nbsearch_make_nblist(const gmx_nbsearch_t nbs,
                              const gmx_nb_atomdata_t *nbat,
                              const t_blocka *excl,
                              real rlist,
                              int min_ci_balanced,
                              gmx_bool nonLocal,
                              int nnbl,gmx_nblist_t **nbl,
                              gmx_bool CombineNBLists);

/* Initialize the non-bonded atom data structure.
 * The enum for nbatXFormat is in the file defining gmx_nb_atomdata_t.
 * Copy the ntypes*ntypes*2 sized nbfp non-bonded parameter list
 * to the atom data structure.
 */
void gmx_nb_atomdata_init(gmx_nb_atomdata_t *nbat,
                          int ntype,const real *nbfp,
                          int XFormat,
                          int nout,
                          gmx_nbat_alloc_t *alloc,
                          gmx_nbat_free_t  *free);

/* Copy the atom types to the non-bonded atom data structure */
void gmx_nb_atomdata_set_atomtypes(gmx_nb_atomdata_t *nbat,
                                   const gmx_nbsearch_t nbs,
                                   const int *type);

/* Copy the charges to the non-bonded atom data structure */
void gmx_nb_atomdata_set_charges(gmx_nb_atomdata_t *nbat,
                                 const gmx_nbsearch_t nbs,
                                 const real *charge);

/* Copy the shift vectors to nbat */
void gmx_nb_atomdata_copy_shiftvec(gmx_bool dynamic_box,
                                   rvec *shift_vec,
                                   gmx_nb_atomdata_t *nbat);

enum { enbatATOMSall, enbatATOMSlocal, enbatATOMSnonlocal };

/* Copy x to nbat->x */
void gmx_nb_atomdata_copy_x_to_nbat_x(const gmx_nbsearch_t nbs,
                                      int enbatATOMS,
                                      rvec *x,
                                      gmx_nb_atomdata_t *nbat);

/* Add the forces stored in nbat to f, zeros the forces in nbat */
void gmx_nb_atomdata_add_nbat_f_to_f(const gmx_nbsearch_t nbs,
                                     int enbatATOMS,
                                     const gmx_nb_atomdata_t *nbat,
                                     rvec *f);

/* Add the fshift force stored in nbat to fshift */
void gmx_nb_atomdata_add_nbat_fshift_to_fshift(const gmx_nb_atomdata_t *nbat,
                                               rvec *fshift);

#ifdef __cplusplus
}
#endif

#endif
