/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team.
 * Copyright (c) 2013-2019,2020, by the GROMACS development team, led by
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
 *  \brief
 *  Data types used internally in the nbnxn_cuda module.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_nbnxm
 */

#ifndef NBNXM_CUDA_TYPES_H
#define NBNXM_CUDA_TYPES_H

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/enumerationhelpers.h"

/*! \brief Macro definining default for the prune kernel's j4 processing concurrency.
 *
 *  The GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY macro allows compile-time override.
 */
#ifndef GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY
#    define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY 4
#endif
/*! \brief Default for the prune kernel's j4 processing concurrency.
 *
 *  Initialized using the #GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY macro which allows compile-time override.
 */
const int c_cudaPruneKernelJ4Concurrency = GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY;

/* TODO: consider moving this to kernel_utils */
/* Convenience defines */
/*! \brief cluster size = number of atoms per cluster. */
static constexpr int c_clSize = c_nbnxnGpuClusterSize;

/*! \brief Electrostatic CUDA kernel flavors.
 *
 *  Types of electrostatics implementations available in the CUDA non-bonded
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
enum eelCu
{
    eelCuCUT,
    eelCuRF,
    eelCuEWALD_TAB,
    eelCuEWALD_TAB_TWIN,
    eelCuEWALD_ANA,
    eelCuEWALD_ANA_TWIN,
    eelCuNR
};

/*! \brief VdW CUDA kernel flavors.
 *
 * The enumerates values correspond to the LJ implementations in the CUDA non-bonded
 * kernels.
 *
 * The column-order of pointers to different electrostatic kernels defined in
 * nbnxn_cuda.cu by the nb_*_kfunc_ptr function pointer table
 * should match the order of enumerated types below.
 */
enum evdwCu
{
    evdwCuCUT,
    evdwCuCUTCOMBGEOM,
    evdwCuCUTCOMBLB,
    evdwCuFSWITCH,
    evdwCuPSWITCH,
    evdwCuEWALDGEOM,
    evdwCuEWALDLB,
    evdwCuNR
};

/* All structs prefixed with "cu_" hold data used in GPU calculations and
 * are passed to the kernels, except cu_timers_t. */
/*! \cond */
typedef struct cu_atomdata cu_atomdata_t;
typedef struct cu_nbparam  cu_nbparam_t;
/*! \endcond */


/** \internal
 * \brief Staging area for temporary data downloaded from the GPU.
 *
 *  The energies/shift forces get downloaded here first, before getting added
 *  to the CPU-side aggregate values.
 */
struct nb_staging_t
{
    //! LJ energy
    float* e_lj = nullptr;
    //! electrostatic energy
    float* e_el = nullptr;
    //! shift forces
    float3* fshift = nullptr;
};

/** \internal
 * \brief Nonbonded atom data - both inputs and outputs.
 */
struct cu_atomdata
{
    //! number of atoms
    int natoms;
    //! number of local atoms
    int natoms_local;
    //! allocation size for the atom data (xq, f)
    int nalloc;

    //! atom coordinates + charges, size natoms
    float4* xq;
    //! force output array, size natoms
    float3* f;

    //! LJ energy output, size 1
    float* e_lj;
    //! Electrostatics energy input, size 1
    float* e_el;

    //! shift forces
    float3* fshift;

    //! number of atom types
    int ntypes;
    //! atom type indices, size natoms
    int* atom_types;
    //! sqrt(c6),sqrt(c12) size natoms
    float2* lj_comb;

    //! shifts
    float3* shift_vec;
    //! true if the shift vector has been uploaded
    bool bShiftVecUploaded;
};

/** \internal
 * \brief Parameters required for the CUDA nonbonded calculations.
 */
struct cu_nbparam
{

    //! type of electrostatics, takes values from #eelCu
    int eeltype;
    //! type of VdW impl., takes values from #evdwCu
    int vdwtype;

    //! charge multiplication factor
    float epsfac;
    //! Reaction-field/plain cutoff electrostatics const.
    float c_rf;
    //! Reaction-field electrostatics constant
    float two_k_rf;
    //! Ewald/PME parameter
    float ewald_beta;
    //! Ewald/PME correction term substracted from the direct-space potential
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
    //! nonbonded parameter table with C6/C12 pairs per atom type-pair, 2*ntype^2 elements
    float* nbfp;
    //! texture object bound to nbfp
    cudaTextureObject_t nbfp_texobj;
    //! nonbonded parameter table per atom type, 2*ntype elements
    float* nbfp_comb;
    //! texture object bound to nbfp_texobj
    cudaTextureObject_t nbfp_comb_texobj;

    /* Ewald Coulomb force table data - accessed through texture memory */
    //! table scale/spacing
    float coulomb_tab_scale;
    //! pointer to the table in the device memory
    float* coulomb_tab;
    //! texture object bound to coulomb_tab
    cudaTextureObject_t coulomb_tab_texobj;
};

/** \internal
 * \brief Pair list data.
 */
using cu_plist_t = Nbnxm::gpu_plist;

/** \internal
 * \brief Typedef of actual timer type.
 */
typedef struct Nbnxm::gpu_timers_t cu_timers_t;

class GpuEventSynchronizer;

/*! \internal
 * \brief Main data structure for CUDA nonbonded force calculations.
 */
struct NbnxmGpu
{
    /*! \brief GPU device context.
     *
     * \todo Make it constant reference, once NbnxmGpu is a proper class.
     */
    const DeviceContext* deviceContext_;
    /*! \brief true if doing both local/non-local NB work on GPU */
    bool bUseTwoStreams = false;
    /*! \brief atom data */
    cu_atomdata_t* atdat = nullptr;
    /*! \brief f buf ops cell index mapping */
    int* cell = nullptr;
    /*! \brief number of indices in cell buffer */
    int ncell = 0;
    /*! \brief number of indices allocated in cell buffer */
    int ncell_alloc = 0;
    /*! \brief array of atom indices */
    int* atomIndices = nullptr;
    /*! \brief size of atom indices */
    int atomIndicesSize = 0;
    /*! \brief size of atom indices allocated in device buffer */
    int atomIndicesSize_alloc = 0;
    /*! \brief x buf ops num of atoms */
    int* cxy_na = nullptr;
    /*! \brief number of elements in cxy_na */
    int ncxy_na = 0;
    /*! \brief number of elements allocated allocated in device buffer */
    int ncxy_na_alloc = 0;
    /*! \brief x buf ops cell index mapping */
    int* cxy_ind = nullptr;
    /*! \brief number of elements in cxy_ind */
    int ncxy_ind = 0;
    /*! \brief number of elements allocated allocated in device buffer */
    int ncxy_ind_alloc = 0;
    /*! \brief parameters required for the non-bonded calc. */
    cu_nbparam_t* nbparam = nullptr;
    /*! \brief pair-list data structures (local and non-local) */
    gmx::EnumerationArray<Nbnxm::InteractionLocality, cu_plist_t*> plist = { { nullptr } };
    /*! \brief staging area where fshift/energies get downloaded */
    nb_staging_t nbst;
    /*! \brief local and non-local GPU streams */
    gmx::EnumerationArray<Nbnxm::InteractionLocality, const DeviceStream*> deviceStreams;

    /*! \brief Events used for synchronization */
    /*! \{ */
    /*! \brief Event triggered when the non-local non-bonded
     * kernel is done (and the local transfer can proceed) */
    cudaEvent_t nonlocal_done = nullptr;
    /*! \brief Event triggered when the tasks issued in the local
     * stream that need to precede the non-local force or buffer
     * operation calculations are done (e.g. f buffer 0-ing, local
     * x/q H2D, buffer op initialization in local stream that is
     * required also by nonlocal stream ) */
    cudaEvent_t misc_ops_and_local_H2D_done = nullptr;
    /*! \} */

    /*! \brief True if there is work for the current domain in the
     * respective locality.
     *
     * This includes local/nonlocal GPU work, either bonded or
     * nonbonded, scheduled to be executed in the current
     * domain. As long as bonded work is not split up into
     * local/nonlocal, if there is bonded GPU work, both flags
     * will be true. */
    gmx::EnumerationArray<Nbnxm::InteractionLocality, bool> haveWork = { { false } };

    /*! \brief Pointer to event synchronizer triggered when the local
     * GPU buffer ops / reduction is complete
     *
     * \note That the synchronizer is managed outside of this module
     * in StatePropagatorDataGpu.
     */
    GpuEventSynchronizer* localFReductionDone = nullptr;

    /*! \brief Event triggered when non-local coordinate buffer
     * has been copied from device to host. */
    GpuEventSynchronizer* xNonLocalCopyD2HDone = nullptr;

    /* NOTE: With current CUDA versions (<=5.0) timing doesn't work with multiple
     * concurrent streams, so we won't time if both l/nl work is done on GPUs.
     * Timer init/uninit is still done even with timing off so only the condition
     * setting bDoTime needs to be change if this CUDA "feature" gets fixed. */
    /*! \brief True if event-based timing is enabled. */
    bool bDoTime = false;
    /*! \brief CUDA event-based timers. */
    cu_timers_t* timers = nullptr;
    /*! \brief Timing data. TODO: deprecate this and query timers for accumulated data instead */
    gmx_wallclock_gpu_nbnxn_t* timings = nullptr;
};

#endif /* NBNXN_CUDA_TYPES_H */
