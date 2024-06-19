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

/*! \internal
 * \ingroup __module_nbnxm
 *
 * \brief Implements classes for accumulating Coulomb and VdW energy contributions
 *
 * \author Berk Hess <hess@kth.se>
 */

#include "simd_energy_accumulator.h"

#include "config.h"

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

//! Returns the smallest p such with 2^p >= n
int log2LargerEq(const int n)
{
    int log2 = 0;
    while (n > (1 << log2))
    {
        log2++;
    }

    return log2;
}

} // namespace

EnergyGroupsPerCluster::EnergyGroupsPerCluster(const int numEnergyGroups, const int iClusterSize) :
    numEnergyGroups2Log_(log2LargerEq(numEnergyGroups)),
    mask_((1 << numEnergyGroups2Log_) - 1),
    iClusterSize_(iClusterSize)
{
}

void EnergyGroupsPerCluster::setEnergyGroups(ArrayRef<const int> energyGroups)
{
    GMX_ASSERT(energyGroups.ssize() % iClusterSize_ == 0,
               "The number of indices should be a multiple of the i-cluster size");
    const int numIClusters = energyGroups.ssize() / iClusterSize_;

    energyGroups_.resize(numIClusters);

    for (int iCluster = 0; iCluster < numIClusters; iCluster++)
    {
        // Store all energy group indices for an i-cluster in one int
        int comb = 0;
        for (int i = 0; i < iClusterSize_; i++)
        {
            comb |= energyGroups[iCluster * iClusterSize_ + i] << (i * numEnergyGroups2Log_);
        }
        energyGroups_[iCluster] = comb;
    }
}

void EnergyAccumulator<false, true>::getEnergies(ArrayRef<real> coulombEnergies,
                                                 ArrayRef<real> vdwEnergies) const
{
    GMX_ASSERT(gmx::ssize(coulombEnergies) == 1, "Buffer should have size 1");
    GMX_ASSERT(gmx::ssize(vdwEnergies) == 1, "Buffer should have size 1");

    coulombEnergies[0] = coulombEnergyReal_;
    vdwEnergies[0]     = vdwEnergyReal_;
}

EnergyAccumulator<true, true>::EnergyAccumulator(const int numEnergyGroups,
                                                 const int iClusterSize,
                                                 const int jClusterSize) :
    numGroups_(numEnergyGroups),
    numGroups2Log_(log2LargerEq(numEnergyGroups)),
    iShift_(log2LargerEq(numEnergyGroups)),
    iMask_((1 << iShift_) - 1),
    jShift_(log2LargerEq(numEnergyGroups) * 2),
#if GMX_SIMD
    jMask_((1 << jShift_) - 1),
#endif
    jStride_((jClusterSize >> 1) * jClusterSize),
    iStride_(numEnergyGroups * (1 << numGroups2Log_) * jStride_),
    energyGroups_(nullptr),
    jClusterSize_(jClusterSize)
{
    const int numBinsUsed =
            square(numEnergyGroups) * (1 << numGroups2Log_) * (jClusterSize / 2) * jClusterSize;

    coulombEnergyGroupPairBins_.resize(numBinsUsed);
    vdwEnergyGroupPairBins_.resize(numBinsUsed);

    coulombBinIAtomPtrs_.resize(iClusterSize);
    vdwBinIAtomPtrs_.resize(iClusterSize);
}

void EnergyAccumulator<true, true>::clearEnergiesAndSetEnergyGroupsForJClusters(
        const EnergyGroupsPerCluster& energyGroupsPerCluster)
{
    std::fill(vdwEnergyGroupPairBins_.begin(), vdwEnergyGroupPairBins_.end(), 0.0_real);
    std::fill(coulombEnergyGroupPairBins_.begin(), coulombEnergyGroupPairBins_.end(), 0.0_real);

    energyGroups_ = energyGroupsPerCluster.energyGroups_.data();
}

template<int jClusterSize>
void EnergyAccumulator<true, true>::getEnergies(ArrayRef<real> coulombEnergies, ArrayRef<real> vdwEnergies) const
{
    // Clear the output buffers
    std::fill(vdwEnergies.begin(), vdwEnergies.end(), 0.0_real);
    std::fill(coulombEnergies.begin(), coulombEnergies.end(), 0.0_real);

    constexpr int c_halfJClusterSize = jClusterSize / 2;
    /* Energies are stored in SIMD registers with size 2^numGroups_2log */
    const int numGroupsStorage = (1 << numGroups2Log_);

    const real* gmx_restrict vVdwSimd     = vdwEnergyGroupPairBins_.data();
    const real* gmx_restrict vCoulombSimd = coulombEnergyGroupPairBins_.data();
    real* gmx_restrict       vVdw         = vdwEnergies.data();
    real* gmx_restrict       vCoulomb     = coulombEnergies.data();

    /* The size of the SIMD energy group buffer array is:
     * numGroups_ * numGroups_ * numGroupsStorage * halfJClusterSize * simdWidth
     */
    for (int i = 0; i < numGroups_; i++)
    {
        for (int j1 = 0; j1 < numGroups_; j1++)
        {
            for (int j0 = 0; j0 < numGroups_; j0++)
            {
                int c = ((i * numGroups_ + j1) * numGroupsStorage + j0) * c_halfJClusterSize * jClusterSize;
                for (int s = 0; s < c_halfJClusterSize; s++)
                {
                    vVdw[i * numGroups_ + j0] += vVdwSimd[c + 0];
                    vVdw[i * numGroups_ + j1] += vVdwSimd[c + 1];
                    vCoulomb[i * numGroups_ + j0] += vCoulombSimd[c + 0];
                    vCoulomb[i * numGroups_ + j1] += vCoulombSimd[c + 1];
                    c += jClusterSize + 2;
                }
            }
        }
    }
}

void EnergyAccumulator<true, true>::getEnergies(ArrayRef<real> coulombEnergies, ArrayRef<real> vdwEnergies) const
{
    GMX_ASSERT(gmx::ssize(coulombEnergies) == numGroups_ * numGroups_,
               "Buffer size should be #energy-groups^2");
    GMX_ASSERT(gmx::ssize(vdwEnergies) == numGroups_ * numGroups_,
               "Buffer size should be #energy-groups^2");

    switch (jClusterSize_)
    {
        case 2: getEnergies<2>(coulombEnergies, vdwEnergies); break;
        case 4: getEnergies<4>(coulombEnergies, vdwEnergies); break;
        case 8: getEnergies<8>(coulombEnergies, vdwEnergies); break;
        default: GMX_ASSERT(false, "j-cluster size not implemented");
    }
}

} // namespace gmx
