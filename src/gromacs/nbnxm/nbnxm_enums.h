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

/*! \libinternal \file
 *
 * \brief Nbnxm internal enumerations.
 *
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_nbnxm
 */


#ifndef GMX_NBNXM_NBNXM_ENUMS_H
#define GMX_NBNXM_NBNXM_ENUMS_H

#include "config.h"

namespace gmx
{

/*! \brief Nbnxm electrostatic GPU kernel flavors.
 *
 *  Types of electrostatics implementations available in the GPU non-bonded
 *  force kernels. These represent both the electrostatics types implemented
 *  by the kernels (cut-off, RF, and Ewald - a subset of what's defined in
 *  enums.h) as well as encode implementation details analytical/tabulated
 *  and single or twin cut-off (for Ewald kernels).
 *  Note that the cut-off and RF kernels have only analytical flavor and unlike
 *  in the CPU kernels, the tabulated kernels are ATM Ewald-only.
 *
 *  The row-order of pointers to different electrostatic kernels defined in
 *  nbnxn_cuda.cu by the nb_*_kfunc_ptr function pointer table
 *  should match the order of enumerated types below.
 */
enum class ElecType : int
{
    Cut,          //!< Plain cut-off
    RF,           //!< Reaction field
    EwaldTab,     //!< Tabulated Ewald with single cut-off
    EwaldTabTwin, //!< Tabulated Ewald with twin cut-off
    EwaldAna,     //!< Analytical Ewald with single cut-off
    EwaldAnaTwin, //!< Analytical Ewald with twin cut-off
    Count         //!< Number of valid values
};

//! Number of possible \ref ElecType values.
constexpr int c_numElecTypes = static_cast<int>(ElecType::Count);

/*! \brief Nbnxm VdW GPU kernel flavors.
 *
 * The enumerates values correspond to the LJ implementations in the GPU non-bonded
 * kernels.
 *
 * The column-order of pointers to different electrostatic kernels defined in
 * nbnxn_cuda_ocl.cpp/.cu by the nb_*_kfunc_ptr function pointer table
 * should match the order of enumerated types below.
 */
enum class VdwType : int
{
    Cut,         //!< Plain cut-off
    CutCombGeom, //!< Cut-off with geometric combination rules
    CutCombLB,   //!< Cut-off with Lorentz-Berthelot combination rules
    FSwitch,     //!< Smooth force switch
    PSwitch,     //!< Smooth potential switch
    EwaldGeom,   //!< Ewald with geometric combination rules
    EwaldLB,     //!< Ewald with Lorentz-Berthelot combination rules
    Count        //!< Number of valid values
};

//! Number of possible \ref VdwType values.
constexpr int c_numVdwTypes = static_cast<int>(VdwType::Count);

/*! \brief Nonbonded NxN kernel types: plain C, CPU SIMD, GPU, GPU emulation */
enum class NbnxmKernelType : int
{
    NotSet = 0,       //<! Legacy leftover
    Cpu4x4_PlainC,    //<! Plain C CPU kernels, only for comparison
    Cpu4xN_Simd_4xN,  //<! SIMD 4N CPU kernels
    Cpu4xN_Simd_2xNN, //<! SIMD 2NN CPU kernels
    Gpu8x8x8,         //<! GPU kernels, will be specialized further later
    Cpu8x8x8_PlainC,
    Count
};

static constexpr bool isGpuKernelType(const NbnxmKernelType kernelType)
{
    return kernelType == NbnxmKernelType::Gpu8x8x8;
}

/*! \brief Ewald exclusion types */
enum class EwaldExclusionType : int
{
    NotSet = 0,
    Table,
    Analytical,
    DecidedByGpuModule
};

const char* nbnxmKernelTypeToName(NbnxmKernelType kernelType);

//! The available pair list types
enum class PairlistType : int
{
    Simple4x2,
    Simple4x4,
    Simple4x8,
    HierarchicalNxN,
    Count
};

//! \brief Kinds of electrostatic treatments in SIMD Verlet kernels
enum class CoulombKernelType : int
{
    ReactionField,
    Table,
    TableTwin,
    Ewald,
    EwaldTwin,
    Count
};

//! The i- and j-cluster size for GPU lists, 8 atoms for CUDA, set at configure time for OpenCL and SYCL
#if GMX_GPU_OPENCL || GMX_GPU_SYCL
constexpr int c_nbnxnGpuClusterSize = GMX_GPU_NB_CLUSTER_SIZE;
#else
constexpr int        c_nbnxnGpuClusterSize      = 8;
#endif

/*! \brief The number of clusters along a direction in a pair-search grid cell for GPU lists
 *
 * Typically all 2, but X can be 1 when targeting Intel Ponte Vecchio */
//! \{
constexpr int c_gpuNumClusterPerCellZ = GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Z;
constexpr int c_gpuNumClusterPerCellY = GMX_GPU_NB_NUM_CLUSTER_PER_CELL_Y;
constexpr int c_gpuNumClusterPerCellX = GMX_GPU_NB_NUM_CLUSTER_PER_CELL_X;
//! \}

/*! \brief The number of sub-parts used for data storage for a GPU cluster pair
 *
 * In CUDA the number of threads in a warp is 32 and we have cluster pairs
 * of 8*8=64 atoms, so it's convenient to store data for cluster pair halves,
 * i.e. split in 2.
 *
 * On architectures with 64-wide execution however it is better to avoid splitting
 * (e.g. AMD GCN, CDNA and later).
 */
#if GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
static constexpr int c_nbnxnGpuClusterpairSplit = 1;
#else
static constexpr int c_nbnxnGpuClusterpairSplit = 2;
#endif

//! The fixed size of the exclusion mask array for a half GPU cluster pair
static constexpr int c_nbnxnGpuExclSize =
        c_nbnxnGpuClusterSize * c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit;

//! The number of clusters in a pair-search grid cell for GPU lists
static constexpr int c_gpuNumClusterPerCell =
        c_gpuNumClusterPerCellZ * c_gpuNumClusterPerCellY * c_gpuNumClusterPerCellX;

/*! \brief The number of clusters in a super-cluster, used for GPU
 *
 * Configured via GMX_GPU_NB_NUM_CLUSTER_PER_CELL_[XYZ] CMake options.
 * Typically 8 (2*2*2), but can be 4 (1*2*2) when targeting Intel Ponte Vecchio. */
constexpr int c_nbnxnGpuNumClusterPerSupercluster =
        c_gpuNumClusterPerCellX * c_gpuNumClusterPerCellY * c_gpuNumClusterPerCellZ;

/*! \brief With GPU kernels we group cluster pairs in 4 to optimize memory usage
 * of integers containing 32 bits.
 */
constexpr int c_nbnxnGpuJgroupSize = (32 / c_nbnxnGpuNumClusterPerSupercluster);

} // namespace gmx

//! \brief Whether have a separate cut-off check for VDW interactions
enum class VdwCutoffCheck : int
{
    No,
    Yes
};

//! \brief Kind of Lennard-Jones Ewald treatments in NBNxM SIMD kernels
enum class LJEwald : int
{
    None,
    CombGeometric
};

//! \brief What kind of energies are output, whole system or for energy groups
enum class EnergyOutput : int
{
    None,
    System,
    GroupPairs
};

#endif // GMX_NBNXM_NBNXM_H
