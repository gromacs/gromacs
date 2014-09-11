/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _nbnxn_atomdata_h
#define _nbnxn_atomdata_h

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Default nbnxn allocation routine, allocates 32 byte aligned,
 * which works for plain C and aligned SSE and AVX loads/stores.
 */
void nbnxn_alloc_aligned(void **ptr, size_t nbytes);

/* Free function for memory allocated with nbnxn_alloc_aligned */
void nbnxn_free_aligned(void *ptr);

/* Reallocation wrapper function for nbnxn data structures */
void nbnxn_realloc_void(void **ptr,
                        int nbytes_copy, int nbytes_new,
                        nbnxn_alloc_t *ma,
                        nbnxn_free_t  *mf);

/* Reallocate the nbnxn_atomdata_t for a size of n atoms */
void nbnxn_atomdata_realloc(nbnxn_atomdata_t *nbat, int n);

/* Copy na rvec elements from x to xnb using nbatFormat, start dest a0,
 * and fills up to na_round using cx,cy,cz.
 */
void copy_rvec_to_nbat_real(const int *a, int na, int na_round,
                            rvec *x, int nbatFormat, real *xnb, int a0,
                            int cx, int cy, int cz);

enum {
    enbnxninitcombruleDETECT, enbnxninitcombruleGEOM, enbnxninitcombruleLB, enbnxninitcombruleNONE
};

/* Initialize the non-bonded atom data structure.
 * The enum for nbatXFormat is in the file defining nbnxn_atomdata_t.
 * Copy the ntypes*ntypes*2 sized nbfp non-bonded parameter list
 * to the atom data structure.
 * enbnxninitcombrule sets what combination rule data gets stored in nbat.
 */
void nbnxn_atomdata_init(FILE *fp,
                         nbnxn_atomdata_t *nbat,
                         int nb_kernel_type,
                         int enbnxninitcombrule,
                         int ntype, const real *nbfp,
                         int n_energygroups,
                         int nout,
                         nbnxn_alloc_t *alloc,
                         nbnxn_free_t  *free);

/* Copy the atom data to the non-bonded atom data structure */
void nbnxn_atomdata_set(nbnxn_atomdata_t    *nbat,
                        int                  locality,
                        const nbnxn_search_t nbs,
                        const t_mdatoms     *mdatoms,
                        const int           *atinfo);

/* Copy the shift vectors to nbat */
void nbnxn_atomdata_copy_shiftvec(gmx_bool          dynamic_box,
                                  rvec             *shift_vec,
                                  nbnxn_atomdata_t *nbat);

/* Copy x to nbat->x.
 * FillLocal tells if the local filler particle coordinates should be zeroed.
 */
void nbnxn_atomdata_copy_x_to_nbat_x(const nbnxn_search_t nbs,
                                     int                  locality,
                                     gmx_bool             FillLocal,
                                     rvec                *x,
                                     nbnxn_atomdata_t    *nbat);

/* Add the forces stored in nbat to f, zeros the forces in nbat */
void nbnxn_atomdata_add_nbat_f_to_f(const nbnxn_search_t    nbs,
                                    int                     locality,
                                    const nbnxn_atomdata_t *nbat,
                                    rvec                   *f);

/* Add the fshift force stored in nbat to fshift */
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t *nbat,
                                              rvec                   *fshift);

#ifdef __cplusplus
}
#endif

#endif
