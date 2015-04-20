/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _vsite_h
#define _vsite_h

#include <stdio.h>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/pbcutil/ishift.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    t_ilist ilist[F_NRE];     /* vsite ilists for this thread            */
    rvec    fshift[SHIFTS];   /* fshift accumulation buffer              */
    matrix  dxdf;             /* virial dx*df accumulation buffer        */
} gmx_vsite_thread_t;

typedef struct {
    gmx_bool            bHaveChargeGroups;    /* Do we have charge groups?               */
    int                 n_intercg_vsite;      /* The number of inter charge group vsites */
    int                 nvsite_pbc_molt;      /* The array size of vsite_pbc_molt        */
    int              ***vsite_pbc_molt;       /* The pbc atoms for intercg vsites        */
    int               **vsite_pbc_loc;        /* The local pbc atoms                     */
    int                *vsite_pbc_loc_nalloc; /* Sizes of vsite_pbc_loc                  */
    int                 nthreads;             /* Number of threads used for vsites       */
    gmx_vsite_thread_t *tdata;                /* Thread local vsites and work structs    */
    int                *th_ind;               /* Work array                              */
    int                 th_ind_nalloc;        /* Size of th_ind                          */
} gmx_vsite_t;

struct t_graph;

void construct_vsites(gmx_vsite_t *vsite,
                      rvec x[],
                      real dt, rvec v[],
                      t_iparams ip[], t_ilist ilist[],
                      int ePBC, gmx_bool bMolPBC,
                      t_commrec *cr, matrix box);
/* Create positions of vsite atoms based on surrounding atoms
 * for the local system.
 * If v is passed, the velocities of the vsites will be calculated
 * as the new positions minus the old positions divided by dt,
 * thus v should only be passed when the coordinates have been
 * updated with a full time step.
 * Note that velocitis of vsites are completely irrelevant
 * for the integration, they are only useful for analysis.
 */

void construct_vsites_mtop(gmx_vsite_t *vsite,
                           gmx_mtop_t *mtop, rvec x[]);
/* Create positions of vsite atoms based on surrounding atoms
 * for the whole system.
 * This function assumes that all molecules are whole.
 */

void spread_vsite_f(gmx_vsite_t *vsite,
                    rvec x[], rvec f[], rvec *fshift,
                    gmx_bool VirCorr, matrix vir,
                    t_nrnb *nrnb, t_idef *idef,
                    int ePBC, gmx_bool bMolPBC, struct t_graph *g, matrix box,
                    t_commrec *cr);
/* Spread the force operating on the vsite atoms on the surrounding atoms.
 * If fshift!=NULL also update the shift forces.
 * If VirCorr=TRUE add the virial correction for non-linear vsite constructs
 * to vir. This correction is required when the virial is not calculated
 * afterwards from the particle position and forces, but in a different way,
 * as for instance for the PME mesh contribution.
 */

gmx_vsite_t *init_vsite(gmx_mtop_t *mtop, t_commrec *cr,
                        gmx_bool bSerial_NoPBC);
/* Initialize the virtual site struct,
 * returns NULL when there are no virtual sites.
 * bSerial_NoPBC is to generate a simple vsite setup to be
 * used only serial (no MPI or thread parallelization) and without pbc;
 * this is useful for correction vsites of the initial configuration.
 */

void split_vsites_over_threads(const t_ilist   *ilist,
                               const t_iparams *ip,
                               const t_mdatoms *mdatoms,
                               gmx_bool         bLimitRange,
                               gmx_vsite_t     *vsite);
/* Divide the vsite work-load over the threads.
 * Should be called at the end of the domain decomposition.
 */

void set_vsite_top(gmx_vsite_t *vsite, gmx_localtop_t *top, t_mdatoms *md,
                   t_commrec *cr);
/* Set some vsite data for runs without domain decomposition.
 * Should be called once after init_vsite, before calling other routines.
 */

#ifdef __cplusplus
}
#endif

#endif
