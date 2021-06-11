/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * Tests for nbnxm setup utilities
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <cmath>

#include "gromacs/utility/arrayref.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "nblib/box.h"
#include "nblib/nbnxmsetuphelpers.h"

#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NbnxmSetupTest, findNumEnergyGroups)
{
    std::vector<int64_t> v(10);
    int                  arbitraryGid = 7;

    // this sets some bit outside the range of bits used for the group ID
    // having all bits zero except those used for the group ID can otherwise hide bugs
    v[5] |= gmx::sc_atomInfo_HasCharge;
    v[5] = (v[5] & ~gmx::sc_atomInfo_EnergyGroupIdMask) | arbitraryGid;

    int nEnergyGroups = arbitraryGid + 1;
    EXPECT_EQ(nEnergyGroups, findNumEnergyGroups(v));
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnumAuto)
{
    auto kernel = SimdKernels::SimdAuto;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::NotSet);
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnumNo)
{
    auto kernel = SimdKernels::SimdNo;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::Cpu4x4_PlainC);
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnum2XM)
{
    auto kernel = SimdKernels::Simd2XMM;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::Cpu4xN_Simd_2xNN);
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnum4XM)
{
    auto kernel = SimdKernels::Simd4XM;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::Cpu4xN_Simd_4xN);
}

TEST(NbnxmSetupTest, CheckKernelSetupThrowsAuto)
{
    EXPECT_ANY_THROW(checkKernelSetup(SimdKernels::SimdAuto));
}

TEST(NbnxmSetupTest, CheckKernelSetupThrowsCount)
{
    EXPECT_ANY_THROW(checkKernelSetup(SimdKernels::Count));
}

TEST(NbnxmSetupTest, canCreateKernelSetupPlain)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd      = SimdKernels::SimdNo;
    Nbnxm::KernelSetup kernelSetup = createKernelSetupCPU(nbKernelOptions);
    EXPECT_EQ(kernelSetup.kernelType, Nbnxm::KernelType::Cpu4x4_PlainC);
    EXPECT_EQ(kernelSetup.ewaldExclusionType, Nbnxm::EwaldExclusionType::Table);
}

TEST(NbnxmSetupTest, canCreateParticleInfoAllVdv)
{
    size_t  numParticles = 2;
    int64_t mask         = 0;
    mask |= gmx::sc_atomInfo_HasVdw;
    mask |= gmx::sc_atomInfo_HasCharge;
    std::vector<int64_t> refParticles  = { mask, mask };
    std::vector<int64_t> testParticles = createParticleInfoAllVdv(numParticles);
    EXPECT_EQ(refParticles, testParticles);
}

TEST(NbnxmSetupTest, ewaldCoeffWorks)
{
    real                              ewald     = ewaldCoeff(1e-5, 1.0);
    gmx::test::FloatingPointTolerance tolerance = gmx::test::absoluteTolerance(1e-5);
    EXPECT_REAL_EQ_TOL(ewald, 3.12341, tolerance);
}

TEST(NbnxmSetupTest, updateForcerecWorks)
{
    t_forcerec forcerec;
    Box        box(3);
    EXPECT_NO_THROW(updateForcerec(&forcerec, box.legacyMatrix()));
}

// The following tests check if the user is allowed to specify configurations not permitted due
// to conflicting compile time setup flags

TEST(NbnxmSetupTest, canCheckKernelSetup)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd = SimdKernels::SimdNo;
#ifdef GMX_NBNXN_SIMD_4XN
    nbKernelOptions.nbnxmSimd = SimdKernels::Simd4XM;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
    nbKernelOptions.nbnxmSimd = SimdKernels::Simd2XMM;
#endif
    EXPECT_NO_THROW(checkKernelSetup(nbKernelOptions.nbnxmSimd));
}

// check if the user is allowed to ask for SimdKernels::Simd2XMM when NBLIB is not compiled with it
#ifndef GMX_NBNXN_SIMD_2XNN
TEST(NbnxmSetupTest, cannotCreateKernelSetupCPU2XM)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd             = SimdKernels::Simd2XMM;
    nbKernelOptions.useTabulatedEwaldCorr = true;
    EXPECT_ANY_THROW(createKernelSetupCPU(nbKernelOptions));
}
#endif

// check if the user is allowed to ask for SimdKernels::Simd4XM when NBLIB is not compiled with it
#ifndef GMX_NBNXN_SIMD_4XN
TEST(NbnxmSetupTest, cannotCreateKernelSetupCPU4XM)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd             = SimdKernels::Simd4XM;
    nbKernelOptions.useTabulatedEwaldCorr = false;
    EXPECT_ANY_THROW(createKernelSetupCPU(nbKernelOptions));
}
#endif

TEST(NbnxmSetupTest, CanCreateNbnxmCPU)
{
    size_t          numParticles = 1;
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd             = SimdKernels::SimdNo;
    int               numEnergyGroups     = 1;
    std::vector<real> nonbondedParameters = { 1, 1 };
    EXPECT_NO_THROW(createNbnxmCPU(numParticles, nbKernelOptions, numEnergyGroups, nonbondedParameters));
}

} // namespace
} // namespace test
} // namespace nblib
