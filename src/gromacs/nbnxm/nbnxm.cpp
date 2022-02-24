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

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/message_string_collector.h"

#include "nbnxm_gpu.h"
#include "pairlistsets.h"
#include "pairsearch.h"

/*! \cond INTERNAL */

void nbnxn_put_on_grid(nonbonded_verlet_t*            nb_verlet,
                       const matrix                   box,
                       int                            gridIndex,
                       const rvec                     lowerCorner,
                       const rvec                     upperCorner,
                       const gmx::UpdateGroupsCog*    updateGroupsCog,
                       gmx::Range<int>                atomRange,
                       real                           atomDensity,
                       gmx::ArrayRef<const int64_t>   atomInfo,
                       gmx::ArrayRef<const gmx::RVec> x,
                       int                            numAtomsMoved,
                       const int*                     move)
{
    nb_verlet->pairSearch_->putOnGrid(box,
                                      gridIndex,
                                      lowerCorner,
                                      upperCorner,
                                      updateGroupsCog,
                                      atomRange,
                                      atomDensity,
                                      atomInfo,
                                      x,
                                      numAtomsMoved,
                                      move,
                                      nb_verlet->nbat.get());
}

/* Calls nbnxn_put_on_grid for all non-local domains */
void nbnxn_put_on_grid_nonlocal(nonbonded_verlet_t*              nbv,
                                const struct gmx_domdec_zones_t* zones,
                                gmx::ArrayRef<const int64_t>     atomInfo,
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

        nbnxn_put_on_grid(nbv,
                          nullptr,
                          zone,
                          c0,
                          c1,
                          nullptr,
                          { zones->cg_range[zone], zones->cg_range[zone + 1] },
                          -1,
                          atomInfo,
                          x,
                          0,
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

gmx::ArrayRef<const int> nonbonded_verlet_t::getLocalAtomOrder() const
{
    /* Return the atom order for the home cell (index 0) */
    const Nbnxm::Grid& grid = pairSearch_->gridSet().grids()[0];

    const int numIndices = grid.atomIndexEnd() - grid.firstAtomInColumn(0);

    return gmx::constArrayRefFromArray(pairSearch_->gridSet().atomIndices().data(), numIndices);
}

void nonbonded_verlet_t::setLocalAtomOrder() const
{
    pairSearch_->setLocalAtomOrder();
}

void nonbonded_verlet_t::setAtomProperties(gmx::ArrayRef<const int>     atomTypes,
                                           gmx::ArrayRef<const real>    atomCharges,
                                           gmx::ArrayRef<const int64_t> atomInfo) const
{
    nbnxn_atomdata_set(nbat.get(), pairSearch_->gridSet(), atomTypes, atomCharges, atomInfo);
}

void nonbonded_verlet_t::convertCoordinates(const gmx::AtomLocality        locality,
                                            gmx::ArrayRef<const gmx::RVec> coordinates)
{
    wallcycle_start(wcycle_, WallCycleCounter::NbXFBufOps);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::NBXBufOps);

    nbnxn_atomdata_copy_x_to_nbat_x(
            pairSearch_->gridSet(), locality, as_rvec_array(coordinates.data()), nbat.get());

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::NBXBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::NbXFBufOps);
}

void nonbonded_verlet_t::convertCoordinatesGpu(const gmx::AtomLocality locality,
                                               DeviceBuffer<gmx::RVec> d_x,
                                               GpuEventSynchronizer*   xReadyOnDevice)
{
    wallcycle_start(wcycle_, WallCycleCounter::LaunchGpu);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuNBXBufOps);

    nbnxn_atomdata_x_to_nbat_x_gpu(pairSearch_->gridSet(), locality, gpu_nbv, d_x, xReadyOnDevice);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuNBXBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);
}

gmx::ArrayRef<const int> nonbonded_verlet_t::getGridIndices() const
{
    return pairSearch_->gridSet().cells();
}

void nonbonded_verlet_t::atomdata_add_nbat_f_to_f(const gmx::AtomLocality  locality,
                                                  gmx::ArrayRef<gmx::RVec> force)
{

    /* Skip the reduction if there was no short-range GPU work to do
     * (either NB or both NB and bonded work). */
    if (!pairlistIsSimple() && !Nbnxm::haveGpuShortRangeWork(gpu_nbv, atomToInteractionLocality(locality)))
    {
        return;
    }

    wallcycle_start(wcycle_, WallCycleCounter::NbXFBufOps);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::NBFBufOps);

    reduceForces(nbat.get(), locality, pairSearch_->gridSet(), as_rvec_array(force.data()));

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::NBFBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::NbXFBufOps);
}

int nonbonded_verlet_t::getNumAtoms(const gmx::AtomLocality locality) const
{
    int numAtoms = 0;
    switch (locality)
    {
        case gmx::AtomLocality::All: numAtoms = pairSearch_->gridSet().numRealAtomsTotal(); break;
        case gmx::AtomLocality::Local: numAtoms = pairSearch_->gridSet().numRealAtomsLocal(); break;
        case gmx::AtomLocality::NonLocal:
            numAtoms = pairSearch_->gridSet().numRealAtomsTotal()
                       - pairSearch_->gridSet().numRealAtomsLocal();
            break;
        case gmx::AtomLocality::Count:
            GMX_ASSERT(false, "Count is invalid locality specifier");
            break;
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

void nonbonded_verlet_t::setupGpuShortRangeWork(const gmx::ListedForcesGpu*    listedForcesGpu,
                                                const gmx::InteractionLocality iLocality) const
{
    if (useGpu() && !emulateGpu())
    {
        Nbnxm::setupGpuShortRangeWork(gpu_nbv, listedForcesGpu, iLocality);
    }
}

void nonbonded_verlet_t::atomdata_init_copy_x_to_nbat_x_gpu() const
{
    Nbnxm::nbnxn_gpu_init_x_to_nbat_x(pairSearch_->gridSet(), gpu_nbv);
}

bool buildSupportsNonbondedOnGpu(std::string* error)
{
    gmx::MessageStringCollector errorReasons;
    // Before changing the prefix string, make sure that it is not searched for in regression tests.
    errorReasons.startContext("Nonbonded interactions on GPUs are not supported in:");
    errorReasons.appendIf(GMX_DOUBLE, "Double precision build of GROMACS");
    errorReasons.appendIf(!GMX_GPU, "Non-GPU build of GROMACS.");
    errorReasons.finishContext();
    if (error != nullptr)
    {
        *error = errorReasons.toString();
    }
    return errorReasons.isEmpty();
}

/*! \endcond */
