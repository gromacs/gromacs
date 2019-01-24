/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief
 * Implements the Nbnxm class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "nbnxm.h"

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/timing/wallcycle.h"

#include "pairlistsets.h"
#include "pairsearch.h"

/*! \cond INTERNAL */

void nbnxn_put_on_grid(nonbonded_verlet_t             *nb_verlet,
                       const matrix                    box,
                       int                             ddZone,
                       const rvec                      lowerCorner,
                       const rvec                      upperCorner,
                       const gmx::UpdateGroupsCog     *updateGroupsCog,
                       int                             atomStart,
                       int                             atomEnd,
                       real                            atomDensity,
                       gmx::ArrayRef<const int>        atomInfo,
                       gmx::ArrayRef<const gmx::RVec>  x,
                       int                             numAtomsMoved,
                       const int                      *move)
{
    nb_verlet->pairSearch_->putOnGrid(box, ddZone, lowerCorner, upperCorner,
                                      updateGroupsCog, atomStart, atomEnd, atomDensity,
                                      atomInfo, x, numAtomsMoved, move,
                                      nb_verlet->nbat.get());
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nonbonded_verlet_t              *nbv,
                                const struct gmx_domdec_zones_t *zones,
                                gmx::ArrayRef<const int>         atomInfo,
                                gmx::ArrayRef<const gmx::RVec>   x)
{
    for (int zone = 1; zone < zones->n; zone++)
    {
        rvec c0, c1;
        for (int d = 0; d < DIM; d++)
        {
            c0[d] = zones->size[zone].bb_x0[d];
            c1[d] = zones->size[zone].bb_x1[d];
        }

        nbnxn_put_on_grid(nbv, nullptr,
                          zone, c0, c1,
                          nullptr,
                          zones->cg_range[zone],
                          zones->cg_range[zone+1],
                          -1,
                          atomInfo,
                          x,
                          0, nullptr);
    }
}

bool nonbonded_verlet_t::isDynamicPruningStepCpu(int64_t step) const
{
    return pairlistSets_->isDynamicPruningStepCpu(step);
}

bool nonbonded_verlet_t::isDynamicPruningStepGpu(int64_t step) const
{
    return pairlistSets_->isDynamicPruningStepGpu(step);
}

gmx::ArrayRef<const int> nonbonded_verlet_t::getLocalAtomOrder() const
{
    /* Return the atom order for the home cell (index 0) */
    const Nbnxm::Grid &grid       = pairSearch_->gridSet().grids()[0];

    const int          numIndices = grid.atomIndexEnd() - grid.firstAtomInColumn(0);

    return gmx::constArrayRefFromArray(pairSearch_->gridSet().atomIndices().data(), numIndices);
}

void nonbonded_verlet_t::setLocalAtomOrder()
{
    pairSearch_->setLocalAtomOrder();
}

void nonbonded_verlet_t::setAtomProperties(const t_mdatoms          &mdatoms,
                                           gmx::ArrayRef<const int>  atomInfo)
{
    nbnxn_atomdata_set(nbat.get(), pairSearch_->gridSet(), &mdatoms, atomInfo.data());
}

void nonbonded_verlet_t::setCoordinates(const Nbnxm::AtomLocality       locality,
                                        const bool                      fillLocal,
                                        gmx::ArrayRef<const gmx::RVec>  x,
                                        bool                            useGpu,
                                        void                           *xPmeDevicePtr,
                                        gmx_wallcycle                  *wcycle)
{
    wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
    wallcycle_sub_start(wcycle, ewcsNB_X_BUF_OPS);
    nbnxn_atomdata_copy_x_to_nbat_x(pairSearch_->gridSet(), locality, fillLocal,
                                    as_rvec_array(x.data()),
                                    nbat.get(), useGpu, gpu_nbv, xPmeDevicePtr);
    wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
    wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
}

gmx::ArrayRef<const int> nonbonded_verlet_t::getGridIndices() const
{
    return pairSearch_->gridSet().cells();
}

void
nonbonded_verlet_t::atomdata_add_nbat_f_to_f(const Nbnxm::AtomLocality  locality,
                                             rvec                      *f,
                                             gmx_wallcycle             *wcycle)
{
    wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
    wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);

    reduceForces(nbat.get(), locality, pairSearch_->gridSet(), f);

    wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
    wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
}

real nonbonded_verlet_t::pairlistInnerRadius() const
{
    return pairlistSets_->params().rlistInner;
}

real nonbonded_verlet_t::pairlistOuterRadius() const
{
    return pairlistSets_->params().rlistOuter;
}

void nonbonded_verlet_t::changePairlistRadii(real rlistOuter,
                                             real rlistInner)
{
    pairlistSets_->changePairlistRadii(rlistOuter, rlistInner);
}

void
nonbonded_verlet_t::atomdata_init_copy_x_to_nbat_x_gpu(const Nbnxm::AtomLocality        locality)
{

    nbnxn_gpu_init_x_to_nbat_x(pairSearch_->gridSet(),
                               gpu_nbv,
                               locality);


}
/*! \endcond */
