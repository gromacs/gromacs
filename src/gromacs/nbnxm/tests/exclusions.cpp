/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for exclusions in the Nbnxm CPU pairlists
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "config.h"

#include <utility>

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistwork.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

/*! \brief Sets up and return atom data for \p numAtoms atoms and a \p numAtoms^2 pair list
 *
 * All atom coordinates are zero. All atom pairs are in the list.
 *
 * Note that this function returns a unique pointer, because parts of
 * PairlistSet currently can not be copied.
 */
std::pair<std::unique_ptr<nbnxn_atomdata_t>, std::unique_ptr<PairlistSet>>
diagonalPairlist(const Nbnxm::KernelType kernelType, const int numAtoms)
{
    const gmx::MDLogger emptyLogger;

    t_commrec commRec;
    commRec.duty = (DUTY_PP | DUTY_PME);

    gmx_omp_nthreads_init(emptyLogger, &commRec, 1, 1, 1, 1, false);

    const PairlistParams pairlistParams(kernelType, false, 1, false);

    Nbnxm::GridSet gridSet(
            PbcType::Xyz, false, nullptr, nullptr, pairlistParams.pairlistType, false, 1, gmx::PinningPolicy::CannotBePinned);

    std::vector<real> nbfp({ 0.0_real, 0.0_real });

    std::unique_ptr<nbnxn_atomdata_t> nbat = std::make_unique<nbnxn_atomdata_t>(
            gmx::PinningPolicy::CannotBePinned, emptyLogger, kernelType, 0, 1, nbfp, 1, 1);

    std::vector<gmx::RVec> coords(numAtoms, { 1.0_real, 1.0_real, 1.0_real });

    matrix box         = { { 3.0_real, 0.0_real, 0.0_real },
                   { 0.0_real, 3.0_real, 0.0_real },
                   { 0.0_real, 0.0_real, 3.0_real } };
    rvec   lowerCorner = { 0.0_real, 0.0_real, 0.0_real };
    rvec   upperCorner = { 3.0_real, 3.0_real, 3.0_real };

    std::vector<int64_t> atomInfo(numAtoms, sc_atomInfo_HasVdw);

    gridSet.putOnGrid(box,
                      0,
                      lowerCorner,
                      upperCorner,
                      nullptr,
                      { 0, numAtoms },
                      numAtoms / det(box),
                      atomInfo,
                      coords,
                      0,
                      nullptr,
                      nbat.get());

    std::unique_ptr<PairlistSet> pairlistSet = std::make_unique<PairlistSet>(pairlistParams);

    std::vector<PairsearchWork> searchWork(1);

    gmx::ListOfLists<int> exclusions;

    for (int i = 0; i < numAtoms; i++)
    {
        exclusions.pushBack({});
    }

    pairlistSet->constructPairlists(
            gmx::InteractionLocality::Local, gridSet, searchWork, nbat.get(), exclusions, 0, nullptr, nullptr);

    return std::make_pair(std::move(nbat), std::move(pairlistSet));
}

// Class that sets up and holds a set of N atoms and a full NxM pairlist
class CpuListDiagonalExclusionsTest : public ::testing::TestWithParam<Nbnxm::KernelType>
{
public:
    CpuListDiagonalExclusionsTest()
    {
        const Nbnxm::KernelType kernelType = GetParam();

        const PairlistParams pairlistParams(kernelType, false, 1, false);

        iClusterSize_ = IClusterSizePerListType[pairlistParams.pairlistType];
        jClusterSize_ = JClusterSizePerListType[pairlistParams.pairlistType];

        numAtoms_ = std::max(iClusterSize_, jClusterSize_);

        std::tie(nbat_, pairlistSet_) = diagonalPairlist(kernelType, numAtoms_);
    }

    const NbnxnPairlistCpu& pairlist() const { return pairlistSet_->cpuLists()[0]; }

    int iClusterSize_;
    int jClusterSize_;
    int numAtoms_;

private:
    std::unique_ptr<nbnxn_atomdata_t> nbat_;
    std::unique_ptr<PairlistSet>      pairlistSet_;
};

// Checks that the correct bits are set for avoiding double counting interactions around the diagonal
TEST_P(CpuListDiagonalExclusionsTest, CheckMask)
{
    ASSERT_EQ(numAtoms_ / iClusterSize_, pairlist().ci.size());

    ASSERT_EQ(numAtoms_ * numAtoms_ / (iClusterSize_ * jClusterSize_), pairlist().cj.size());

    for (const auto& iEntry : pairlist().ci)
    {
        const int iCluster = iEntry.ci;

        for (int cjIndex = iEntry.cj_ind_start; cjIndex < iEntry.cj_ind_end; cjIndex++)
        {
            const int          jCluster = pairlist().cj.list_[cjIndex].cj;
            const unsigned int excl     = pairlist().cj.list_[cjIndex].excl;

            for (int iIndex = 0; iIndex < iClusterSize_; iIndex++)
            {
                const int iAtom = iCluster * iClusterSize_ + iIndex;

                for (int jIndex = 0; jIndex < jClusterSize_; jIndex++)
                {
                    const int jAtom = jCluster * jClusterSize_ + jIndex;

                    EXPECT_EQ((excl >> (iIndex * jClusterSize_ + jIndex)) & 1, (jAtom > iAtom ? 1 : 0));
                }
            }
        }
    }
}

const auto testKernelTypes = ::testing::Values(Nbnxm::KernelType::Cpu4x4_PlainC
#if GMX_HAVE_NBNXM_SIMD_4XM
                                               ,
                                               Nbnxm::KernelType::Cpu4xN_Simd_4xN
#endif
#if GMX_HAVE_NBNXM_SIMD_2XMM
                                               ,
                                               Nbnxm::KernelType::Cpu4xN_Simd_2xNN
#endif
);

INSTANTIATE_TEST_SUITE_P(WithParameters, CpuListDiagonalExclusionsTest, testKernelTypes);

} // namespace
} // namespace test
} // namespace gmx
