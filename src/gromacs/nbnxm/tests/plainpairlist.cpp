/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Tests for the plain pairlist extraction functionality
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/kernel_common.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/vec.h"

#include "testutils/testasserts.h"
#include "testutils/testinit.h"
#include "testutils/testoptions.h"

#include "testsystem.h"

namespace gmx
{

namespace test
{

namespace
{

int g_numOpenMPThreads = 2;

//! \cond
GMX_TEST_OPTIONS(KernelBenchmarkOptions, options)
{
    options->addOption(
            IntegerOption("ntomp").store(&g_numOpenMPThreads).description("Number of OpenMP threads to use"));
}
//! \endcond

//! The options for the kernel
struct KernelOptions
{
    //! Constructor
    KernelOptions(NbnxmKernelType kernelType) :
        useGpu(kernelType == NbnxmKernelType::Gpu8x8x8),
        numThreads(GMX_OPENMP && !useGpu ? g_numOpenMPThreads : 1),
        kernelSetup({ kernelType, EwaldExclusionType::NotSet })
    {
    }

    //! Whether to use a GPU
    bool useGpu;
    //! The number of OpenMP threads to use
    int numThreads;
    //! The kernel setup
    NbnxmKernelSetup kernelSetup;
    //! The normal pairlist range
    real rlist = 0.8;
    //! The plain pairlist range, to avoid rounding issue the cut-off is not exactly 0.7
    real plainPairlistRange = 0.7001;
};

//! Sets up and returns a Nbnxm object for the given benchmark options and system
std::unique_ptr<nonbonded_verlet_t> setupNbnxmForBenchInstance(const KernelOptions& options,
                                                               const TestSystem&    system)
{
    real minBoxSize = norm(system.box[XX]);
    for (int dim = YY; dim < DIM; dim++)
    {
        minBoxSize = std::min(minBoxSize, norm(system.box[dim]));
    }
    GMX_RELEASE_ASSERT(options.rlist <= 0.5 * minBoxSize,
                       "The cut-off should be shorter than half the box size");

    // We don't want to call gmx_omp_nthreads_init(), so we init what we need
    gmx_omp_nthreads_set(ModuleMultiThread::Pairsearch, options.numThreads);
    gmx_omp_nthreads_set(ModuleMultiThread::Nonbonded, options.numThreads);

    // Don't try pinning buffers; we don't need that in this test
    const auto pinPolicy  = PinningPolicy::CannotBePinned;
    const int  numThreads = options.numThreads;

    PairlistParams pairlistParams(
            options.kernelSetup.kernelType, sc_layoutType, false, options.rlist, false);

    GridSet gridSet(
            PbcType::Xyz, false, nullptr, nullptr, pairlistParams.pairlistType, false, false, numThreads, pinPolicy);

    auto pairlistSets = std::make_unique<PairlistSets>(pairlistParams, false, 0, pinPolicy);

    auto pairSearch = std::make_unique<PairSearch>(
            PbcType::Xyz, false, nullptr, nullptr, pairlistParams.pairlistType, false, false, numThreads, pinPolicy);

    auto atomData = std::make_unique<nbnxn_atomdata_t>(pinPolicy,
                                                       MDLogger(),
                                                       options.kernelSetup.kernelType,
                                                       LJCombinationRule::None,
                                                       LJCombinationRule::Geometric,
                                                       system.nonbondedParameters,
                                                       true,
                                                       1,
                                                       numThreads);

    // Put everything together
    auto nbv = std::make_unique<nonbonded_verlet_t>(
            std::move(pairlistSets), std::move(pairSearch), std::move(atomData), options.kernelSetup, nullptr);

    GMX_RELEASE_ASSERT(!TRICLINIC(system.box), "Only rectangular unit-cells are supported here");
    const rvec lowerCorner = { 0, 0, 0 };
    const rvec upperCorner = { system.box[XX][XX], system.box[YY][YY], system.box[ZZ][ZZ] };

    const real atomDensity = system.coordinates.size() / det(system.box);

    nbv->putAtomsOnGrid(system.box,
                        0,
                        lowerCorner,
                        upperCorner,
                        nullptr,
                        { 0, int(system.coordinates.size()) },
                        system.coordinates.size(),
                        atomDensity,
                        system.atomInfo,
                        system.coordinates,
                        nullptr);

    nbv->constructPairlist(gmx::InteractionLocality::Local, system.excls, true, 0, nullptr);

    nbv->setAtomProperties(system.atomTypes, system.charges, system.atomInfo);

    return nbv;
}

int countPairsWithinCutoff(ArrayRef<const PlainPairlist::PairlistEntry> plainPairlist,
                           ArrayRef<const RVec>                         coordinates,
                           ArrayRef<const RVec>                         shiftVecs,
                           const real                                   cutoff)
{
    const real cutoffSquared = cutoff * cutoff;

    int numPairsWithinCutoff  = 0;
    int numPairsOutsideCutoff = 0;

    for (const auto& [atomPair, shiftIndex] : plainPairlist)
    {
        const RVec d = coordinates[atomPair.first] - coordinates[atomPair.second] + shiftVecs[shiftIndex];

        if (norm2(d) < cutoffSquared)
        {
            numPairsWithinCutoff++;
        }
        else
        {
            numPairsOutsideCutoff++;
        }
    }

    EXPECT_EQ(numPairsOutsideCutoff, 0);

    return numPairsWithinCutoff;
}

//! Class that sets up and holds a set of N atoms and a full NxM pairlist
class PlainPairlistTest : public ::testing::TestWithParam<NbnxmKernelType>
{
};

} // namespace

//! Test case that checks that the NBNxM kernel produces correct output.
TEST_P(PlainPairlistTest, ContainsAllPairs)
{
    // The test parameters with which the test case was instantiated
    // TODO rename these in a follow-up change to conform to style
    const NbnxmKernelType kernelType = GetParam();
    const KernelOptions   options(kernelType);

    if (!sc_haveNbnxmSimd4xmKernels && kernelType == NbnxmKernelType::Cpu4xN_Simd_4xN)
    {
        GTEST_SKIP()
                << "Cannot test or generate data for 4xN kernels without suitable SIMD support";
    }

    if (!sc_haveNbnxmSimd2xmmKernels && kernelType == NbnxmKernelType::Cpu4xN_Simd_2xNN)
    {
        GTEST_SKIP()
                << "Cannot test or generate data for 2xNN kernels without suitable SIMD support";
    }

    // TODO rename this in a follow-up change to conform to style
    const TestSystem system_(LJCombinationRule::Geometric, true);

    // Finish setting up data structures
    std::unique_ptr<nonbonded_verlet_t> nbv = setupNbnxmForBenchInstance(options, system_);

    nbv->constructPairlist(InteractionLocality::Local, system_.excls, true, 0, nullptr);

    std::vector<RVec> shiftVecs(c_numShiftVectors);
    calc_shifts(system_.box, shiftVecs);

    const auto& plainPairlist = nbv->plainPairlist(options.plainPairlistRange, shiftVecs);

    const int c_numPairsRef         = 6331;
    const int c_numExcludedPairsRef = 243;

    const int numPairs = countPairsWithinCutoff(
            plainPairlist.pairs, system_.coordinates, shiftVecs, options.plainPairlistRange);

    EXPECT_EQ(numPairs, c_numPairsRef);

    const int numExcludedPairs = countPairsWithinCutoff(
            plainPairlist.excludedPairs, system_.coordinates, shiftVecs, options.plainPairlistRange);

    EXPECT_EQ(numExcludedPairs, c_numExcludedPairsRef);
};

INSTANTIATE_TEST_SUITE_P(WithParameters,
                         PlainPairlistTest,
                         ::testing::Values(NbnxmKernelType::Cpu4x4_PlainC,
                                           NbnxmKernelType::Cpu4xN_Simd_4xN,
                                           NbnxmKernelType::Cpu4xN_Simd_2xNN,
                                           NbnxmKernelType::Gpu8x8x8));

} // namespace test

} // namespace gmx
