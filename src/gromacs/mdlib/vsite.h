/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_VSITE_H
#define GMX_MDLIB_VSITE_H

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_graph;
struct t_ilist;
struct t_mdatoms;
struct t_nrnb;

typedef struct gmx_vsite_t {
    gmx_bool             bHaveChargeGroups;    /* Do we have charge groups?               */
    int                  n_intercg_vsite;      /* The number of inter charge group vsites */
    int                  nvsite_pbc_molt;      /* The array size of vsite_pbc_molt        */
    int               ***vsite_pbc_molt;       /* The pbc atoms for intercg vsites        */
    int                **vsite_pbc_loc;        /* The local pbc atoms                     */
    int                 *vsite_pbc_loc_nalloc; /* Sizes of vsite_pbc_loc                  */
    int                  nthreads;             /* Number of threads used for vsites       */
    struct VsiteThread **tData;                /* Thread local vsites and work structs    */
    int                 *taskIndex;            /* Work array                              */
    int                  taskIndexNalloc;      /* Size of taskIndex                       */
    bool                 useDomdec;            /* Tells whether we use domain decomposition with more than 1 DD rank */

} gmx_vsite_t;

/*! \brief Create positions of vsite atoms based for the local system
 *
 * \param[in]     vsite    The vsite struct, when nullptr is passed, no MPI and no multi-threading is used
 * \param[in,out] x        The coordinates
 * \param[in]     dt       The time step
 * \param[in,out] v        When != nullptr, velocities for vsites are set as displacement/dt
 * \param[in]     ip       Interaction parameters
 * \param[in]     ilist    The interaction list
 * \param[in]     ePBC     The type of periodic boundary conditions
 * \param[in]     bMolPBC  When true, molecules are broken over PBC
 * \param[in]     cr       The communication record
 * \param[in]     box      The box
 */
void construct_vsites(const gmx_vsite_t *vsite,
                      rvec x[],
                      real dt, rvec v[],
                      const t_iparams ip[], const t_ilist ilist[],
                      int ePBC, gmx_bool bMolPBC,
                      t_commrec *cr,
                      const matrix box);

/*! \brief Create positions of vsite atoms for the whole system assuming all molecules are wholex
 *
 * \param[in]     mtop  The global topology
 * \param[in,out] x     The global coordinates
 */
void constructVsitesGlobal(const gmx_mtop_t         &mtop,
                           gmx::ArrayRef<gmx::RVec>  x);

void spread_vsite_f(const gmx_vsite_t *vsite,
                    const rvec x[],
                    rvec f[], rvec *fshift,
                    gmx_bool VirCorr, matrix vir,
                    t_nrnb *nrnb, const t_idef *idef,
                    int ePBC, gmx_bool bMolPBC, const t_graph *g, const matrix box,
                    t_commrec *cr);
/* Spread the force operating on the vsite atoms on the surrounding atoms.
 * If fshift!=NULL also update the shift forces.
 * If VirCorr=TRUE add the virial correction for non-linear vsite constructs
 * to vir. This correction is required when the virial is not calculated
 * afterwards from the particle position and forces, but in a different way,
 * as for instance for the PME mesh contribution.
 */

int count_intercg_vsites(const gmx_mtop_t *mtop);
/* Returns the number of virtual sites that cross charge groups */

/* Initialize the virtual site struct,
 *
 * \param[in] mtop  The global topology
 * \param[in] cr    The communication record
 * \returns A valid vsite struct or nullptr when there are no virtual sites
 */
gmx_vsite_t *initVsite(const gmx_mtop_t &mtop,
                       t_commrec        *cr);

void split_vsites_over_threads(const t_ilist   *ilist,
                               const t_iparams *ip,
                               const t_mdatoms *mdatoms,
                               gmx_vsite_t     *vsite);
/* Divide the vsite work-load over the threads.
 * Should be called at the end of the domain decomposition.
 */

void set_vsite_top(gmx_vsite_t          *vsite,
                   const gmx_localtop_t *top,
                   const t_mdatoms      *md);
/* Set some vsite data for runs without domain decomposition.
 * Should be called once after init_vsite, before calling other routines.
 */

#endif
