/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * \brief Translation layer to GROMACS data structures for force calculations.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/gmxsetup.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/smalloc.h"
#include "nblib/exception.h"
#include "nblib/kerneloptions.h"
#include "nblib/nbnxmsetuphelpers.h"
#include "nblib/particletype.h"
#include "nblib/simulationstate.h"

namespace nblib
{

NbvSetupUtil::NbvSetupUtil() : gmxForceCalculator_(std::make_unique<GmxForceCalculator>()) {}

void NbvSetupUtil::setExecutionContext(const NBKernelOptions& options)
{
    setGmxNonBondedNThreads(options.numOpenMPThreads);
}

void NbvSetupUtil::setParticleInfoAllVdv(const size_t numParticles)

{
    particleInfoAllVdw_ = createParticleInfoAllVdw(numParticles);
}

void NbvSetupUtil::setNonBondedParameters(const std::vector<ParticleType>& particleTypes,
                                          const NonBondedInteractionMap&   nonBondedInteractionMap)
{
    nonbondedParameters_ = createNonBondedParameters(particleTypes, nonBondedInteractionMap);
}

void NbvSetupUtil::setAtomProperties(const std::vector<int>&  particleTypeIdOfAllParticles,
                                     const std::vector<real>& charges)
{
    gmxForceCalculator_->nbv_->setAtomProperties(particleTypeIdOfAllParticles, charges, particleInfoAllVdw_);
}

//! Sets up and returns a Nbnxm object for the given options and system
void NbvSetupUtil::setupNbnxmInstance(const size_t numParticleTypes, const NBKernelOptions& options)
{

    gmxForceCalculator_->nbv_ = createNbnxmCPU(numParticleTypes, options, 1, nonbondedParameters_);
}

void NbvSetupUtil::setupStepWorkload(const NBKernelOptions& options)
{
    gmxForceCalculator_->stepWork_ = std::make_unique<gmx::StepWorkload>(createStepWorkload(options));
}

void NbvSetupUtil::setupInteractionConst(const NBKernelOptions& options)
{
    gmxForceCalculator_->interactionConst_ =
            std::make_unique<interaction_const_t>(createInteractionConst(options));
}

void NbvSetupUtil::setupForceRec(const matrix& box)
{
    updateForcerec(gmxForceCalculator_->forcerec_.get(), box);
}

void NbvSetupUtil::setParticlesOnGrid(const std::vector<Vec3>& coordinates, const Box& box)
{
    gmxForceCalculator_->setParticlesOnGrid(particleInfoAllVdw_, coordinates, box);
}

void NbvSetupUtil::constructPairList(ExclusionLists<int> exclusionLists)
{
    gmx::ListOfLists<int> exclusions(std::move(exclusionLists.ListRanges),
                                     std::move(exclusionLists.ListElements));
    gmxForceCalculator_->nbv_->constructPairlist(
            gmx::InteractionLocality::Local, exclusions, 0, gmxForceCalculator_->nrnb_.get());
}


std::unique_ptr<GmxForceCalculator> setupGmxForceCalculator(const Topology&          topology,
                                                            const std::vector<Vec3>& coordinates,
                                                            const Box&               box,
                                                            const NBKernelOptions&   options)
{
    NbvSetupUtil nbvSetupUtil;
    nbvSetupUtil.setExecutionContext(options);
    nbvSetupUtil.setNonBondedParameters(topology.getParticleTypes(),
                                        topology.getNonBondedInteractionMap());
    nbvSetupUtil.setParticleInfoAllVdv(topology.numParticles());

    nbvSetupUtil.setupInteractionConst(options);
    nbvSetupUtil.setupStepWorkload(options);
    nbvSetupUtil.setupNbnxmInstance(topology.getParticleTypes().size(), options);
    nbvSetupUtil.setParticlesOnGrid(coordinates, box);
    nbvSetupUtil.constructPairList(topology.exclusionLists());
    nbvSetupUtil.setAtomProperties(topology.getParticleTypeIdOfAllParticles(), topology.getCharges());
    nbvSetupUtil.setupForceRec(box.legacyMatrix());

    return nbvSetupUtil.getGmxForceCalculator();
}

} // namespace nblib
