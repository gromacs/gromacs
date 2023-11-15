/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * \brief Implements common internal types for different NBNXN GPU implementations
 *
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \ingroup module_nbnxm
 */

#ifndef GMX_MDLIB_NBNXN_GPU_COMMON_TYPES_H
#define GMX_MDLIB_NBNXN_GPU_COMMON_TYPES_H

#include "config.h"

#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/enumerationhelpers.h"

#include "nbnxm.h"
#include "pairlist.h"

#if GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gpuregiontimer_ocl.h"
#endif

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gpuregiontimer.cuh"
#endif

#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gpuregiontimer_sycl.h"
#endif

/*! \brief Number of separate bins used during sorting of plist on gpu
 *
 * Ideally this number would be increased for very large system sizes (the cpu version of sorting
 * uses 2 x avg(num cjPacked) but as sorting has negligible impact for very large system sizes we
 * use a constant here for simplicity. On H100 sorting begins to have negligible effect for
 * system sizes greater than ~400k atoms.
 */
static constexpr int c_sciHistogramSize = 8192;

/*! \brief Number of threads per block used by the gpu sorting kernel
 *
 * TODO this is a reasonable default but the number has not been tuned
 */
static constexpr int c_sciSortingThreadsPerBlock = 256;

/*! \brief Macro definining default for the prune kernel's jPacked processing concurrency.
 *
 *  The GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY macro allows compile-time override with the default value of 4.
 */
#ifndef GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY
#    define GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY 4
#endif
//! Default for the prune kernel's jPacked processing concurrency.
static constexpr int c_pruneKernelJPackedConcurrency = GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY;

/*! \internal
 * \brief Staging area for temporary data downloaded from the GPU.
 *
 * Since SYCL buffers already have host-side storage, this is a bit redundant.
 * But it allows prefetching of the data from GPU, and brings GPU backends closer together.
 */
struct NBStagingData
{
    //! LJ energy
    float* eLJ = nullptr;
    //! electrostatic energy
    float* eElec = nullptr;
    //! shift forces
    Float3* fShift = nullptr;
};

/** \internal
 * \brief Nonbonded atom data - both inputs and outputs.
 */
struct NBAtomDataGpu
{
    //! number of atoms
    int numAtoms;
    //! number of local atoms
    int numAtomsLocal;
    //! allocation size for the atom data (xq, f)
    int numAtomsAlloc;

    //! atom coordinates + charges, size \ref numAtoms
    DeviceBuffer<Float4> xq;
    //! force output array, size \ref numAtoms
    DeviceBuffer<Float3> f;

    //! LJ energy output, size 1
    DeviceBuffer<float> eLJ;
    //! Electrostatics energy input, size 1
    DeviceBuffer<float> eElec;

    //! shift forces
    DeviceBuffer<Float3> fShift;

    //! number of atom types
    int numTypes;
    //! atom type indices, size \ref numAtoms
    DeviceBuffer<int> atomTypes;
    //! sqrt(c6),sqrt(c12) size \ref numAtoms
    DeviceBuffer<Float2> ljComb;

    //! shifts
    DeviceBuffer<Float3> shiftVec;
    //! true if the shift vector has been uploaded
    bool shiftVecUploaded;
};

/** \internal
 * \brief Parameters required for the GPU nonbonded calculations.
 */
struct NBParamGpu
{

    //! type of electrostatics
    enum Nbnxm::ElecType elecType;
    //! type of VdW impl.
    enum Nbnxm::VdwType vdwType;

    //! charge multiplication factor
    float epsfac;
    //! Reaction-field/plain cutoff electrostatics const.
    float c_rf;
    //! Reaction-field electrostatics constant
    float two_k_rf;
    //! Ewald/PME parameter
    float ewald_beta;
    //! Ewald/PME correction term subtracted from the direct-space potential
    float sh_ewald;
    //! LJ-Ewald/PME correction term added to the correction potential
    float sh_lj_ewald;
    //! LJ-Ewald/PME coefficient
    float ewaldcoeff_lj;

    //! Coulomb cut-off squared
    float rcoulomb_sq;

    //! VdW cut-off squared
    float rvdw_sq;
    //! VdW switched cut-off
    float rvdw_switch;
    //! Full, outer pair-list cut-off squared
    float rlistOuter_sq;
    //! Inner, dynamic pruned pair-list cut-off squared
    float rlistInner_sq;
    //! True if we use dynamic pair-list pruning
    bool useDynamicPruning;

    //! VdW shift dispersion constants
    shift_consts_t dispersion_shift;
    //! VdW shift repulsion constants
    shift_consts_t repulsion_shift;
    //! VdW switch constants
    switch_consts_t vdw_switch;

    /* LJ non-bonded parameters - accessed through texture memory */
    //! nonbonded parameter table with 6*C6/12*C12 pairs per atom type-pair, ntype^2 elements
    DeviceBuffer<Float2> nbfp{};
    //! texture object bound to nbfp
    DeviceTexture nbfp_texobj;
    //! nonbonded parameter table per atom type, ntype elements
    DeviceBuffer<Float2> nbfp_comb{};
    //! texture object bound to nbfp_comb
    DeviceTexture nbfp_comb_texobj;

    /* Ewald Coulomb force table data - accessed through texture memory */
    //! table scale/spacing
    float coulomb_tab_scale;
    //! pointer to the table in the device memory
    DeviceBuffer<float> coulomb_tab{};
    //! texture object bound to coulomb_tab
    DeviceTexture coulomb_tab_texobj;
};

namespace Nbnxm
{

using gmx::AtomLocality;
using gmx::InteractionLocality;

/*! \internal
 * \brief GPU region timers used for timing GPU kernels and H2D/D2H transfers.
 *
 * The two-sized arrays hold the local and non-local values and should always
 * be indexed with eintLocal/eintNonlocal.
 */
struct GpuTimers
{
    /*! \internal
     * \brief Timers for local or non-local coordinate/force transfers
     */
    struct XFTransfers
    {
        //! timer for x/q H2D transfers (l/nl, every step)
        GpuRegionTimer nb_h2d;
        //! timer for f D2H transfer (l/nl, every step)
        GpuRegionTimer nb_d2h;
    };

    /*! \internal
     * \brief Timers for local or non-local interaction related operations
     */
    struct Interaction
    {
        //! timer for pair-list H2D transfers (l/nl, every PS step)
        GpuRegionTimer pl_h2d;
        //! true when a pair-list transfer has been done at this step
        bool didPairlistH2D = false;
        //! timer for non-bonded kernels (l/nl, every step)
        GpuRegionTimer nb_k;
        //! timer for the 1st pass list pruning kernel (l/nl, every PS step)
        GpuRegionTimer prune_k;
        //! true when we timed pruning and the timings need to be accounted for
        bool didPrune = false;
        //! timer for rolling pruning kernels (l/nl, frequency depends on chunk size)
        GpuRegionTimer rollingPrune_k;
        //! true when we timed rolling pruning (at the previous step) and the timings need to be accounted for
        bool didRollingPrune = false;
    };

    //! timer for atom data transfer (every PS step)
    GpuRegionTimer atdat;
    //! timers for coordinate/force transfers (every step)
    gmx::EnumerationArray<AtomLocality, XFTransfers> xf;
    //! timers for interaction related transfers
    gmx::EnumerationArray<InteractionLocality, Nbnxm::GpuTimers::Interaction> interaction;
};


/*! \internal
 * \brief Sorted pair list on GPU and data required for performing the sorting */
struct gpuPlistSorting
{
    //! size of scanTemporary, working array used for exclusive prefix sum calculation
    int nscanTemporary;

    //! allocation size of scanTemporary
    int scanTemporaryNalloc;

    // //! Temporary data of scan algorithm
    DeviceBuffer<char> scanTemporary;

    //! number of buckets in histogram
    int nsciHistogram;

    //! allocation size of sciHistogram
    int sciHistogramNalloc;

    //! Histogram of sci nsp
    DeviceBuffer<int> sciHistogram;

    //! size of sciOffset, number of histogram buckets
    int nsciOffset;

    //! allocation size of sciOffset
    int sciOffsetNalloc;

    //! Sci offset, the exclusive prefix sum of sciHistogram
    DeviceBuffer<int> sciOffset;

    //! size of sci, # of i clusters in the list
    int nsciCounted;

    //! allocation size of sci
    int sciCountedNalloc;

    //! list of imask counts of sorted i-cluster ("super-clusters")
    DeviceBuffer<int> sciCount;

    //! size of sci, # of i clusters in the list
    int nsciSorted;
    //! allocation size of sci
    int sciSortedNalloc;

    //! list of sorted i-cluster ("super-clusters")
    DeviceBuffer<nbnxn_sci_t> sciSorted;
};

/*! \internal
 * \brief GPU pair list structure */
struct gpu_plist
{
    //! number of atoms per cluster
    int na_c;

    //! size of sci, # of i clusters in the list
    int nsci;
    //! allocation size of sci
    int sci_nalloc;
    //! list of i-cluster ("super-clusters")
    DeviceBuffer<nbnxn_sci_t> sci;

    //! sorted pair list and data used for sorting
    gpuPlistSorting sorting;

    //! total # of packed j clusters
    int ncjPacked;
    //! allocation size of cjPacked
    int cjPacked_nalloc;
    //! Packed j cluster list, contains j cluster number and index into the i cluster list
    DeviceBuffer<nbnxn_cj_packed_t> cjPacked;
    //! # of packed j clusters * # of warps
    int nimask;
    //! allocation size of imask
    int imask_nalloc;
    //! imask for 2 warps for each 4*j cluster group
    DeviceBuffer<unsigned int> imask;
    //! atom interaction bits
    DeviceBuffer<nbnxn_excl_t> excl;
    //! count for excl
    int nexcl;
    //! allocation size of excl
    int excl_nalloc;

    /* parameter+variables for normal and rolling pruning */
    //! true after search, indicates that initial pruning with outer pruning is needed
    bool haveFreshList;
    //! the number of parts/steps over which one cycle of rolling pruning takes places
    int rollingPruningNumParts;
    //! the next part to which the rolling pruning needs to be applied
    int rollingPruningPart;
    //! device memory buffer (1 value per thread block) for next part to which the rolling pruning needs to be applied
    DeviceBuffer<int> d_rollingPruningPart;
    //! size of rolling pruning part buffer on device
    int d_rollingPruningPart_size = -1;
    //! allocated size of rolling pruning part buffer on device
    int d_rollingPruningPart_size_alloc = -1;
};

} // namespace Nbnxm

#endif
