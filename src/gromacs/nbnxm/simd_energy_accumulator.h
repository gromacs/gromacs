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

/*! \internal
 * \ingroup __module_nbnxm
 *
 * \brief Defines classes for accumulating Coulomb and VdW energy contributions
 *
 * There is class for no energy output which does nothing, but is useful for
 * avoiding lots of conditionals in the kernel. The classes for a single and
 * multiple energy groups should be used with the following call sequence:
 *
 * - clearEnergies...()
 * - for (i-cluster ...) {
 * -   initICluster()
 * -   for (j-cluster ...) {
 * -     addEnergy/ies()
 * -   }
 * -   reduceIEnergies()
 * - }
 * - getEnergies()
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef GMX_NBNXM_ENERGY_ACCUMULATOR_H
#define GMX_NBNXM_ENERGY_ACCUMULATOR_H

#include <cstddef>
#include <cstdint>

#include <array>
#include <vector>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "nbnxm_simd.h"

namespace gmx
{

class EnergyGroupsPerCluster;

#if GMX_SIMD

//! Adds energies to temporary energy group pair buffers for the 4xM kernel layout
template<std::size_t offsetJJSize>
inline void accumulateGroupPairEnergies4xM(SimdReal energies,
                                           real*    groupPairEnergyBuffersPtr,
                                           const std::array<int, offsetJJSize>& offsetJJ)
{
    static_assert(offsetJJSize == GMX_SIMD_REAL_WIDTH / 2);

    /* We need to balance the number of store operations with
     * the rapidly increasing number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for (int jj = 0; jj < GMX_SIMD_REAL_WIDTH / 2; jj++)
    {
        SimdReal groupPairEnergyBuffers =
                load<SimdReal>(groupPairEnergyBuffersPtr + offsetJJ[jj] + jj * GMX_SIMD_REAL_WIDTH);

        store(groupPairEnergyBuffersPtr + offsetJJ[jj] + jj * GMX_SIMD_REAL_WIDTH,
              groupPairEnergyBuffers + energies);
    }
}

//! Adds energies to temporary energy group pair buffers for the 2xMM kernel layout
template<std::size_t offsetJJSize>
inline void accumulateGroupPairEnergies2xMM(SimdReal energies,
                                            real*    groupPairEnergyBuffersPtr0,
                                            real*    groupPairEnergyBuffersPtr1,
                                            const std::array<int, offsetJJSize>& offsetJJ)
{
    static_assert(offsetJJSize == GMX_SIMD_REAL_WIDTH / 4);

    for (int jj = 0; jj < GMX_SIMD_REAL_WIDTH / 4; jj++)
    {
        incrDualHsimd(groupPairEnergyBuffersPtr0 + offsetJJ[jj] + jj * GMX_SIMD_REAL_WIDTH / 2,
                      groupPairEnergyBuffersPtr1 + offsetJJ[jj] + jj * GMX_SIMD_REAL_WIDTH / 2,
                      energies);
    }
}

#endif // GMX_SIMD

//! Base energy accumulator class, only specializations are used
template<bool, bool>
class EnergyAccumulator;

//! Specialized energy accumulator class for no energy calculation
template<bool useEnergyGroups>
class EnergyAccumulator<useEnergyGroups, false>
{
public:
    //! Does nothing
    template<int iClusterSize>
    inline void initICluster(const int gmx_unused iCluster)
    {
    }

#if GMX_SIMD
    //! Does nothing
    template<int nRCoulomb, int nRVdw, KernelLayout kernelLayout, int iClusterSize, std::size_t cSize, std::size_t vdwSize>
    inline void addEnergies(const int gmx_unused              jCluster,
                            const std::array<SimdReal, cSize> gmx_unused& coulombEnergy,
                            const std::array<SimdReal, vdwSize> gmx_unused& vdwEnergy)
    {
    }
#endif // GMX_SIMD

    //! Does nothing
    inline void addCoulombEnergy(const int gmx_unused iAtomInCluster, const real gmx_unused energy)
    {
    }

    //! Does nothing
    inline void addVdwEnergy(const int gmx_unused iAtomInCluster, const real gmx_unused energy) {}

    //! Does nothing
    inline void reduceIEnergies(const bool gmx_unused calculateCoulomb) {}
};

/*! \brief Specialized energy accumulator class for energy accumulation without energy groups
 *
 * Note that this specialization accumulates over each j-list to internal buffers with an entry
 * per i-particle and then reduces to the final buffers. This is done as to mimimize the rounding
 * errors in the reductions.
 */
template<>
class EnergyAccumulator<false, true>
{
public:
    //! Clears all energy buffers and sets the energy group indices for the j-clusters
    inline void clearEnergies()
    {
        coulombEnergyReal_ = 0;
        vdwEnergyReal_     = 0;
    }

#if GMX_SIMD
    //! Clears the accumulation buffer which is used per i-cluster
    template<int iClusterSize>
    inline void initICluster(const int gmx_unused iCluster)
    {
        coulombEnergySum_ = setZero();
        vdwEnergySum_     = setZero();
    }

    //! Adds a single Coulomb energy contribution
    inline void addCoulombEnergy(const int gmx_unused iAtomInCluster, const real energy)
    {
        coulombEnergyReal_ += energy;
    }

    //! Adds a single VdW energy contribution
    inline void addVdwEnergy(const int gmx_unused iAtomInCluster, const real energy)
    {
        vdwEnergyReal_ += energy;
    }

    /*! \brief Adds Coulomb and/or VdW contributions for interactions of a j-cluster with an i-cluster
     *
     * \tparam nRCoulomb     The number of elements to reduce in \p coulombEnergy
     * \tparam nRVdw         The number of elements to reduce in \p vdwEnergy
     * \tparam kernelLayout  The kernel layout
     * \tparam cSize         The size of \p coulombEnergy should be >= \p nRCoulomb
     * \tparam vdwSize       The size of \p vdwEnergy should be >= \p nRVdw
     *
     * \param jCluster       The index of the j-cluster to add energies for
     * \param coulombEnergy  List of Coulomb energies per i-register
     * \param vdwEnergy      List of Van der Waals energies per i-register
     */
    template<int nRCoulomb, int nRVdw, KernelLayout kernelLayout, int iClusterSize, std::size_t cSize, std::size_t vdwSize>
    inline void addEnergies(const int gmx_unused                 jCluster,
                            const std::array<SimdReal, cSize>&   coulombEnergy,
                            const std::array<SimdReal, vdwSize>& vdwEnergy)
    {
        static_assert(cSize >= nRCoulomb);
        static_assert(vdwSize >= nRVdw);

        for (int i = 0; i < nRCoulomb; i++)
        {
            coulombEnergySum_ = coulombEnergySum_ + coulombEnergy[i];
        }
        for (int i = 0; i < nRVdw; i++)
        {
            vdwEnergySum_ = vdwEnergySum_ + vdwEnergy[i];
        }
    }

    //! Reduces the accumulated energies to the final output buffer
    inline void reduceIEnergies(const bool calculateCoulomb)
    {
        if (calculateCoulomb)
        {
            coulombEnergyReal_ += reduce(coulombEnergySum_);
        }

        vdwEnergyReal_ += reduce(vdwEnergySum_);
    }
#endif // GMX_SIMD

    /*! \brief Returns the internally stored energies in the output buffers.
     *
     * \param  coulombEnergies  Buffer of Coulomb energies to accumulate to, should have size 1
     * \param  vdwEnergies      Buffer of VdW energies to accumulate to, should have size 1
     */
    void getEnergies(ArrayRef<real> coulombEnergies, ArrayRef<real> vdwEnergies) const;

private:
#if GMX_SIMD
    //! Coulomb energy accumulation buffers for a j-list for one i-cluster
    SimdReal coulombEnergySum_;
    //! VdW energy accumulation buffers for a j-list for one i-cluster
    SimdReal vdwEnergySum_;
#endif // GMX_SIMD
    //! Single Coulomb energy accumulation buffer
    real coulombEnergyReal_;
    //! Single VdW energy accumulation buffer
    real vdwEnergyReal_;
};

/*! \brief Specialized energy accumulator class for energy accumulation with energy groups
 *
 * Sums energies into a temporary buffer with bins for each combination of an i-atom energy group
 * with a pair of energy groups for two j-atoms. Reduction of this list of bins into the final
 * energy group pair matrix is done outside the non-bonded kernel.
 */
template<>
class EnergyAccumulator<true, true>
{
public:
    /*! \brief Constructor
     *
     * \param numEnergyGroups  The number of energy groups
     * \param iClusterSize     The i-cluster size
     * \param jClusterSize     The j-cluster size
     */
    EnergyAccumulator(int numEnergyGroups, int iClusterSize, int jClusterSize);

    //! Clears all energy buffers and sets the energy group indices for the j-clusters
    void clearEnergiesAndSetEnergyGroupsForJClusters(const EnergyGroupsPerCluster& energyGroupsPerCluster);

    //! Sets (internal) parameters for the atom in i-cluster \p iCluster
    template<int iClusterSize>
    inline void initICluster(const int iCluster)
    {
        energyGroupsICluster_ = energyGroups_[iCluster];
        for (int iAtom = 0; iAtom < iClusterSize; iAtom++)
        {
            const int iAtomIndex        = (energyGroupsICluster_ >> (iAtom * iShift_)) & iMask_;
            coulombBinIAtomPtrs_[iAtom] = coulombEnergyGroupPairBins_.data() + iAtomIndex * iStride_;
            vdwBinIAtomPtrs_[iAtom]     = vdwEnergyGroupPairBins_.data() + iAtomIndex * iStride_;
        }
    }

    //! Adds a single Coulomb energy contribution for atom with index in cluster: \p iAtomInCluster
    inline void addCoulombEnergy(const int iAtomInCluster, const real energy)
    {
        const int pairIndex = ((energyGroupsICluster_ >> (iAtomInCluster * iShift_)) & iMask_) * jStride_;

        coulombBinIAtomPtrs_[iAtomInCluster][pairIndex] += energy;
    }

    //! Adds a single VdW energy contribution for atom with index in cluster: \p iAtomInCluster
    inline void addVdwEnergy(const int iAtomInCluster, const real energy)
    {
        const int pairIndex = ((energyGroupsICluster_ >> (iAtomInCluster * iShift_)) & iMask_) * jStride_;

        vdwBinIAtomPtrs_[iAtomInCluster][pairIndex] += energy;
    }

#if GMX_SIMD
    /*! \brief Adds Coulomb and/or VdW contributions for interactions of a j-cluster with an i-cluster
     *
     * \tparam nRCoulomb     The number of elements to reduce in \p coulombEnergy
     * \tparam nRVdw         The number of elements to reduce in \p vdwEnergy
     * \tparam kernelLayout  The kernel layout
     * \tparam cSize         The size of \p coulombEnergy should be >= \p nRCoulomb
     * \tparam vdwSize       The size of \p vdwEnergy should be >= \p nRVdw
     *
     * \param jCluster       The index of the j-cluster to add energies for
     * \param coulombEnergy  List of Coulomb energies per i-register
     * \param vdwEnergy      List of Van der Waals energies per i-register
     */
    template<int nRCoulomb, int nRVdw, KernelLayout kernelLayout, int iClusterSize, std::size_t cSize, std::size_t vdwSize>
    inline void addEnergies(const int                            jCluster,
                            const std::array<SimdReal, cSize>&   coulombEnergy,
                            const std::array<SimdReal, vdwSize>& vdwEnergy)
    {
        static_assert(cSize >= nRCoulomb);
        static_assert(vdwSize >= nRVdw);

        constexpr int jClusterSize =
                (kernelLayout == KernelLayout::r4xM ? GMX_SIMD_REAL_WIDTH : GMX_SIMD_REAL_WIDTH / 2);

        /* Energy group indices for two atom pairs packed into one int, one int for each i-atom */
        std::array<int, jClusterSize / 2> ijGroupPair;

        /* Extract the group index per j pair.
         * Energy groups are stored per i-cluster, so things get
         * complicated when the i- and j-cluster sizes don't match.
         */
        if constexpr (jClusterSize == 2)
        {
            const int jPairGroups = energyGroups_[jCluster >> 1];
            ijGroupPair[0] = ((jPairGroups >> ((jCluster & 1) * jShift_)) & jMask_) * jStride_;
        }
        else
        {
            static_assert(iClusterSize <= jClusterSize);

            for (int jdi = 0; jdi < jClusterSize / iClusterSize; jdi++)
            {
                const int jPairGroups = energyGroups_[jCluster * (jClusterSize / iClusterSize) + jdi];
                for (int jj = 0; jj < (iClusterSize / 2); jj++)
                {
                    ijGroupPair[jdi * (iClusterSize / 2) + jj] =
                            ((jPairGroups >> (jj * jShift_)) & jMask_) * jStride_;
                }
            }
        }

        for (int i = 0; i < nRCoulomb; i++)
        {
            if constexpr (jClusterSize == GMX_SIMD_REAL_WIDTH)
            {
                accumulateGroupPairEnergies4xM(coulombEnergy[i], coulombBinIAtomPtrs_[i], ijGroupPair);
            }
            else
            {
                accumulateGroupPairEnergies2xMM(coulombEnergy[i],
                                                coulombBinIAtomPtrs_[i * 2],
                                                coulombBinIAtomPtrs_[i * 2 + 1],
                                                ijGroupPair);
            }
        }

        for (int i = 0; i < nRVdw; i++)
        {
            if constexpr (jClusterSize == GMX_SIMD_REAL_WIDTH)
            {
                accumulateGroupPairEnergies4xM(vdwEnergy[i], vdwBinIAtomPtrs_[i], ijGroupPair);
            }
            else
            {
                accumulateGroupPairEnergies2xMM(
                        vdwEnergy[i], vdwBinIAtomPtrs_[i * 2], vdwBinIAtomPtrs_[i * 2 + 1], ijGroupPair);
            }
        }
    }
#endif

    //! Nothing do to here, reduction happens after the kernel call
    inline void reduceIEnergies(const bool gmx_unused calculateCoulomb) {}

    /*! \brief Reduce the group-pair energy buffers produced by a SIMD kernels
     * and return the results in the output buffers.
     *
     * The SIMD kernels produce a large number of energy buffer in SIMD registers
     * to avoid scattered reads and writes.
     *
     * \param coulombEnergies  Buffer of Coulomb energies to accumulate to
     * \param vdwEnergies      Buffer of VdW energies to accumulate to
     */
    void getEnergies(ArrayRef<real> coulombEnergies, ArrayRef<real> vdwEnergies) const;

private:
    //! As getEnergies above, but templated on jClusterSize
    template<int jClusterSize>
    void getEnergies(ArrayRef<real> coulombEnergies, ArrayRef<real> vdwEnergies) const;

    //! The number of energy groups
    const int numGroups_;
    //! The base 2 log of number of energy groups, rounded up
    const int numGroups2Log_;
    //! i-cluster shift
    const int iShift_;
    //! i-cluster mask
    const int iMask_;
    //! j-cluster shift
    const int jShift_;
#if GMX_SIMD
    //! j-cluster mask
    const int jMask_;
#endif
    //! j-cluster stride
    const int jStride_;
    //! Major division is over i-particle energy groups, this gives the bin stride for an i-atom
    const int iStride_;

    //! Pointer to a list of energy groups for j-clusters, packed into an int
    const int* energyGroups_;

    //! The complete list of Coulomb energy bins for all energy group pair combinations
    AlignedVector<real> coulombEnergyGroupPairBins_;
    //! The complete list of VdW energy bins for all energy group pair combinations
    AlignedVector<real> vdwEnergyGroupPairBins_;

    //! Energy groups for the i-cluster, packed into an int
    int energyGroupsICluster_;
    //! Pointers to the Coulomb energy bins for the atoms in the current i-cluster
    std::vector<real*> coulombBinIAtomPtrs_;
    //! Pointers to the VdW energy bins for the atoms in the current i-cluster
    std::vector<real*> vdwBinIAtomPtrs_;

    //! The size of the j-clusters, only needed for the reduction at the end
    const int jClusterSize_;
};

//! Holds energy group indices for use in EnergyAccumulator<true, true>
class EnergyGroupsPerCluster
{
public:
    /*! \brief Constructor
     *
     * \param numEnergyGroups  The number of non-bonded energy groups
     * \param iClusterSize     The size in atoms of an i-cluster
     */
    EnergyGroupsPerCluster(int numEnergyGroups, int iClusterSize);

    /*! Sets all energy groups for a direct list of energy groups
     *
     * \param[in] energyGroups  List of energy groups, size should a multiple of the i-cluster size
     */
    void setEnergyGroups(ArrayRef<const int> energyGroups);

    /*! \brief Sets (part of) the energy groups starting at \p offset using indices
     *
     * Note that entries of indices are allowed to be -1, those get assigned energy group 0.
     *
     * \param[in] indices   Atom indices to set the energy groups for, size should be a multiple of the i-cluster size
     * \param[in] atomInfo  List of atom info which should cover at least all indexed atoms
     * \param[in] mask      Mask to extract the energy group from an entry in \p atomInfo
     * \param[in] clusterOffset  The cluster offset from where to start storing the energy groups
     */
    void setEnergyGroups(ArrayRef<const int>     indices,
                         ArrayRef<const int32_t> atomInfo,
                         const int               mask,
                         const int               clusterOffset)
    {
        GMX_ASSERT(indices.ssize() % iClusterSize_ == 0,
                   "The number of indices should be a multiple of the i-cluster size");
        const int numIClusters = indices.ssize() / iClusterSize_;

        if (clusterOffset + numIClusters > gmx::ssize(energyGroups_))
        {
            energyGroups_.resize(clusterOffset + numIClusters);
        }

        for (int iCluster = 0; iCluster < numIClusters; iCluster++)
        {
            // Store all energy group indices for an i-cluster in one int
            int comb = 0;
            for (int i = 0; i < iClusterSize_; i++)
            {
                const int index = indices[iCluster * iClusterSize_ + i];
                if (index >= 0)
                {
                    comb |= (atomInfo[index] & mask) << (i * numEnergyGroups2Log_);
                }
            }
            energyGroups_[clusterOffset + iCluster] = comb;
        }
    }

    //! Returns the energy group for atom \p atomIndexInCluster in i-cluster \p iCLuster
    int getEnergyGroup(int iCluster, int atomIndexInCluster)
    {
        return (energyGroups_[iCluster] >> (atomIndexInCluster * numEnergyGroups2Log_)) & mask_;
    }

private:
    //! 2log(nenergrp)
    int numEnergyGroups2Log_;
    //! The energy groups, one int entry per cluster, only set when needed
    AlignedVector<int> energyGroups_;
    //! The mask to obtain an energy group for one atom from the combined groups
    int mask_;
    //! The size in atoms of an i-cluster
    int iClusterSize_;

    friend void EnergyAccumulator<true, true>::clearEnergiesAndSetEnergyGroupsForJClusters(
            const EnergyGroupsPerCluster& energyGroupsPerCluster);
};

} // namespace gmx

#endif // GMX_NBNXM_ENERGY_ACCUMULATOR_H
