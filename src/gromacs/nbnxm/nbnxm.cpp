/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include "gromacs/domdec/domdec_zones.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/message_string_collector.h"

#include "nbnxm_gpu.h"
#include "pairlistsets.h"
#include "pairsearch.h"

/*! \cond INTERNAL */

namespace gmx
{

bool nonbonded_verlet_t::localAtomOrderMatchesNbnxmOrder() const
{
    return pairSearch_->gridSet().localAtomOrderMatchesNbnxmOrder();
}

void nonbonded_verlet_t::putAtomsOnGrid(const matrix            box,
                                        int                     gridIndex,
                                        const RVec&             lowerCorner,
                                        const RVec&             upperCorner,
                                        const UpdateGroupsCog*  updateGroupsCog,
                                        Range<int>              atomRange,
                                        int                     numAtomsWithoutFillers,
                                        real                    atomDensity,
                                        ArrayRef<const int32_t> atomInfo,
                                        ArrayRef<const RVec>    x,
                                        const int*              move)
{
    pairSearch_->putOnGrid(box,
                           gridIndex,
                           lowerCorner,
                           upperCorner,
                           updateGroupsCog,
                           atomRange,
                           numAtomsWithoutFillers,
                           atomDensity,
                           atomInfo,
                           x,
                           move,
                           nbat_.get());
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nonbonded_verlet_t*     nbv,
                                const DomdecZones&      zones,
                                ArrayRef<const int32_t> atomInfo,
                                ArrayRef<const RVec>    x)
{
    for (int zone = 1; zone < zones.numZones(); zone++)
    {
        nbv->putAtomsOnGrid(nullptr,
                            zone,
                            zones.sizes(zone).bb_x0,
                            zones.sizes(zone).bb_x1,
                            nullptr,
                            zones.atomRange(zone),
                            *zones.atomRange(zone).end(),
                            -1,
                            atomInfo,
                            x,
                            nullptr);
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

ArrayRef<const int> nonbonded_verlet_t::getLocalAtomOrder() const
{
    /* Return the atom order for the home cell (index 0) */
    const Grid& grid = pairSearch_->gridSet().grid(0);

    const int numIndices = grid.atomIndexEnd() - grid.firstAtomInColumn(0);

    return constArrayRefFromArray(pairSearch_->gridSet().atomIndices().data(), numIndices);
}

void nonbonded_verlet_t::setLocalAtomOrder() const
{
    pairSearch_->setLocalAtomOrder();
}

void nonbonded_verlet_t::setAtomProperties(ArrayRef<const int>     atomTypes,
                                           ArrayRef<const real>    atomCharges,
                                           ArrayRef<const int32_t> atomInfo) const
{
    nbnxn_atomdata_set(nbat_.get(), pairSearch_->gridSet(), atomTypes, atomCharges, atomInfo);
}

void nonbonded_verlet_t::convertCoordinates(const AtomLocality locality, ArrayRef<const RVec> coordinates)
{
    wallcycle_start(wcycle_, WallCycleCounter::NbXFBufOps);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::NBXBufOps);

    nbnxn_atomdata_copy_x_to_nbat_x(
            pairSearch_->gridSet(), locality, as_rvec_array(coordinates.data()), nbat_.get());

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::NBXBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::NbXFBufOps);
}

void nonbonded_verlet_t::convertCoordinatesGpu(const AtomLocality    locality,
                                               DeviceBuffer<RVec>    d_x,
                                               GpuEventSynchronizer* xReadyOnDevice)
{
    wallcycle_start(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuNBXBufOps);

    nbnxn_atomdata_x_to_nbat_x_gpu(pairSearch_->gridSet(), locality, gpuNbv_, d_x, xReadyOnDevice);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuNBXBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

ArrayRef<const int> nonbonded_verlet_t::getGridIndices() const
{
    return pairSearch_->gridSet().cells();
}

ArrayRef<const int> nonbonded_verlet_t::getLocalGridNumAtomsPerColumn() const
{
    return pairSearch_->gridSet().getLocalGridNumAtomsPerColumn();
}

void nonbonded_verlet_t::atomdata_add_nbat_f_to_f(const AtomLocality locality, ArrayRef<RVec> force)
{

    /* Skip the reduction if there was no short-range GPU work to do
     * (either NB or both NB and bonded work). */
    if (!pairlistIsSimple() && !haveGpuShortRangeWork(gpuNbv_, atomToInteractionLocality(locality)))
    {
        return;
    }

    wallcycle_start(wcycle_, WallCycleCounter::NbXFBufOps);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::NBFBufOps);

    nbat_->reduceForces(locality, pairSearch_->gridSet(), force);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::NBFBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::NbXFBufOps);
}

int nonbonded_verlet_t::getNumAtoms(const AtomLocality locality) const
{
    int numAtoms = 0;

    const GridSet& gridSet = pairSearch_->gridSet();

    if (gridSet.localAtomOrderMatchesNbnxmOrder())
    {
        switch (locality)
        {
            case AtomLocality::All: numAtoms = gridSet.numGridAtomsTotal(); break;
            case AtomLocality::Local: numAtoms = gridSet.numGridAtomsLocal(); break;
            case AtomLocality::NonLocal:
                numAtoms = gridSet.numGridAtomsTotal() - gridSet.numGridAtomsLocal();
                break;
            case AtomLocality::Count:
                GMX_ASSERT(false, "Count is invalid locality specifier");
                break;
        }
    }
    else
    {
        switch (locality)
        {
            case AtomLocality::All: numAtoms = gridSet.numRealAtomsTotal(); break;
            case AtomLocality::Local: numAtoms = gridSet.numRealAtomsLocal(); break;
            case AtomLocality::NonLocal:
                numAtoms = gridSet.numRealAtomsTotal() - gridSet.numRealAtomsLocal();
                break;
            case AtomLocality::Count:
                GMX_ASSERT(false, "Count is invalid locality specifier");
                break;
        }
    }

    return numAtoms;
}

real nonbonded_verlet_t::pairlistInnerRadius() const
{
    return pairlistSets_->params().rlistInner;
}

real nonbonded_verlet_t::pairlistOuterRadius() const
{
    return pairlistSets_->params().rlistOuter;
}

void nonbonded_verlet_t::changePairlistRadii(real rlistOuter, real rlistInner) const
{
    pairlistSets_->changePairlistRadii(rlistOuter, rlistInner);
}

void nonbonded_verlet_t::setupGpuShortRangeWork(const ListedForcesGpu*    listedForcesGpu,
                                                const InteractionLocality iLocality) const
{
    if (useGpu() && !emulateGpu())
    {
        setupGpuShortRangeWorkLow(gpuNbv_, listedForcesGpu, iLocality);
    }
}

void nonbonded_verlet_t::atomdata_init_copy_x_to_nbat_x_gpu() const
{
    nbnxn_gpu_init_x_to_nbat_x(pairSearch_->gridSet(), gpuNbv_);
}

const Grid& nonbonded_verlet_t::localGrid() const
{
    return pairSearch_->gridSet().grid(0);
}

void nonbonded_verlet_t::setNonLocalGrid(const int                           gridIndex,
                                         const int                           ddZone,
                                         const GridDimensions&               gridDimensions,
                                         ArrayRef<const std::pair<int, int>> columns,
                                         ArrayRef<const int32_t>             atomInfo,
                                         ArrayRef<const RVec>                x)
{
    pairSearch_->setNonLocalGrid(gridIndex, ddZone, gridDimensions, columns, atomInfo, x, nbat_.get());
}

std::optional<std::string> nbnxmGpuClusteringDescription()
{
#if GMX_GPU
    return formatString("super-cluster %dx%dx%d / cluster %d (cluster-pair splitting %s)",
                        GMX_GPU_NB_NUM_CLUSTER_PER_CELL_X,
                        GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Y,
                        GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Z,
                        GMX_GPU_NB_CLUSTER_SIZE,
                        GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT ? "off" : "on");
#else
    return std::nullopt;
#endif
}

} // namespace gmx

/*! \endcond */
