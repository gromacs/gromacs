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

#include "internal.h"

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
                       const int                      *atinfo,
                       gmx::ArrayRef<const gmx::RVec>  x,
                       int                             numAtomsMoved,
                       const int                      *move)
{
    nb_verlet->pairSearch_->putOnGrid(box, ddZone, lowerCorner, upperCorner,
                                      updateGroupsCog, atomStart, atomEnd, atomDensity,
                                      atinfo, x, numAtomsMoved, move,
                                      nb_verlet->nbat.get());
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nonbonded_verlet_t              *nbv,
                                const struct gmx_domdec_zones_t *zones,
                                const int                       *atinfo,
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
                          atinfo,
                          x,
                          0, nullptr);
    }
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

void nonbonded_verlet_t::setAtomProperties(const t_mdatoms &mdatoms,
                                           const int       &atinfo)
{
    nbnxn_atomdata_set(nbat.get(), *pairSearch_, &mdatoms, &atinfo);
}

void nonbonded_verlet_t::setCoordinates(const Nbnxm::AtomLocality       locality,
                                        const bool                      fillLocal,
                                        gmx::ArrayRef<const gmx::RVec>  x,
                                        gmx_wallcycle                  *wcycle)
{
    nbnxn_atomdata_copy_x_to_nbat_x(*pairSearch_, locality, fillLocal,
                                    as_rvec_array(x.data()),
                                    nbat.get(), wcycle);
}

void nonbonded_verlet_t::getLocalNumCells(int *numCellsX,
                                          int *numCellsY) const
{
    pairSearch_->gridSet().getLocalNumCells(numCellsX, numCellsY);
}

gmx::ArrayRef<const int> nonbonded_verlet_t::getGridIndices() const
{
    return pairSearch_->gridSet().cells();
}

/*! \endcond */
