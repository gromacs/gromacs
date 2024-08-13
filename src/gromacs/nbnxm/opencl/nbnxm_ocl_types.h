/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief
 *  Data types used internally in the nbnxm_ocl module.
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Szilárd Páll <pszilard@kth.se>
 *  \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_NBNXM_OPENCL_TYPES_H
#define GMX_NBNXM_NBNXM_OPENCL_TYPES_H

#include <memory>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

struct gmx_wallclock_gpu_nbnxn_t;

/*! \brief Pruning kernel flavors.
 *
 * The values correspond to the first call of the pruning post-list generation
 * and the rolling pruning, respectively.
 */
enum ePruneKind
{
    epruneFirst,
    epruneRolling,
    ePruneNR
};

namespace gmx
{

/*! \internal
 * \brief Data structure shared between the OpenCL device code and OpenCL host code
 *
 * Must not contain OpenCL objects (buffers)
 * TODO: review, improve */
typedef struct cl_nbparam_params
{

    //! type of electrostatics
    enum ElecType elecType;
    //! type of VdW impl.
    enum VdwType vdwType;

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

    //! VdW shift dispersion constants
    shift_consts_t dispersion_shift;
    //! VdW shift repulsion constants
    shift_consts_t repulsion_shift;
    //! VdW switch constants
    switch_consts_t vdw_switch;

    /* Ewald Coulomb force table data - accessed through texture memory */
    //! table scale/spacing
    float coulomb_tab_scale;
} cl_nbparam_params_t;


/*! \internal
 * \brief Main data structure for OpenCL nonbonded force calculations.
 */
struct NbnxmGpu
{
    /* \brief OpenCL device context
     *
     * \todo Make it constant reference, once NbnxmGpu is a proper class.
     */
    const DeviceContext* deviceContext_;
    //! OpenCL runtime data (context, kernels)
    struct gmx_device_runtime_data_t* dev_rundata = nullptr;

    /**< Pointers to non-bonded kernel functions
     * organized similar with nb_kfunc_xxx arrays in nbnxn_ocl.cpp */
    ///@{
    cl_kernel kernel_noener_noprune_ptr[c_numElecTypes][c_numVdwTypes] = { { nullptr } };
    cl_kernel kernel_ener_noprune_ptr[c_numElecTypes][c_numVdwTypes]   = { { nullptr } };
    cl_kernel kernel_noener_prune_ptr[c_numElecTypes][c_numVdwTypes]   = { { nullptr } };
    cl_kernel kernel_ener_prune_ptr[c_numElecTypes][c_numVdwTypes]     = { { nullptr } };
    ///@}
    //! prune kernels, ePruneKind defined the kernel kinds
    cl_kernel kernel_pruneonly[ePruneNR] = { nullptr };

    //! true if prefetching fg i-atom LJ parameters should be used in the kernels
    bool bPrefetchLjParam = false;

    /**< auxiliary kernels implementing memset-like functions */
    ///@{
    cl_kernel kernel_memset_f  = nullptr;
    cl_kernel kernel_memset_f2 = nullptr;
    cl_kernel kernel_memset_f3 = nullptr;
    ///@}

    //! true if doing both local/non-local NB work on GPU
    bool bUseTwoStreams = false;
    //! true indicates that the nonlocal_done event was marked
    bool bNonLocalStreamDoneMarked = false;

    //! atom data
    NBAtomDataGpu* atdat = nullptr;
    //! parameters required for the non-bonded calc.
    NBParamGpu* nbparam = nullptr;
    //! pair-list data structures (local and non-local)
    gmx::EnumerationArray<InteractionLocality, std::unique_ptr<GpuPairlist>> plist = { nullptr };
    //! staging area where fshift/energies get downloaded
    NBStagingData nbst;

    // Data for GPU-side coordinate conversion between integrator and NBNXM
    /*! \brief array of atom indices */
    DeviceBuffer<int> atomIndices;
    /*! \brief size of atom indices */
    int atomIndicesSize = 0;
    /*! \brief size of atom indices allocated in device buffer */
    int atomIndicesSize_alloc = 0;
    /*! \brief x buf ops num of atoms */
    DeviceBuffer<int> cxy_na;
    /*! \brief number of elements in cxy_na */
    int ncxy_na = 0;
    /*! \brief number of elements allocated allocated in device buffer */
    int ncxy_na_alloc = 0;
    /*! \brief x buf ops cell index mapping */
    DeviceBuffer<int> cxy_ind;
    /*! \brief number of elements in cxy_ind */
    int ncxy_ind = 0;
    /*! \brief number of elements allocated allocated in device buffer */
    int ncxy_ind_alloc = 0;

    //! local and non-local GPU queues
    gmx::EnumerationArray<InteractionLocality, const DeviceStream*> deviceStreams;

    /*! \brief Events used for synchronization */
    /*! \{ */
    /*! \brief Event triggered when the non-local non-bonded
     * kernel is done (and the local transfer can proceed) */
    GpuEventSynchronizer nonlocal_done;
    /*! \brief Event triggered when the tasks issued in the local
     * stream that need to precede the non-local force or buffer
     * operation calculations are done (e.g. f buffer 0-ing, local
     * x/q H2D, buffer op initialization in local stream that is
     * required also by nonlocal stream ) */
    GpuEventSynchronizer misc_ops_and_local_H2D_done;
    /*! \} */

    //! True if there has been local/nonlocal GPU work, either bonded or nonbonded, scheduled
    //  to be executed in the current domain. As long as bonded work is not split up into
    //  local/nonlocal, if there is bonded GPU work, both flags will be true.
    gmx::EnumerationArray<InteractionLocality, bool> haveWork;


    //! True if event-based timing is enabled.
    bool bDoTime = false;
    //! OpenCL event-based timers.
    GpuTimers* timers = nullptr;
    //! Timing data. TODO: deprecate this and query timers for accumulated data instead
    gmx_wallclock_gpu_nbnxn_t* timings = nullptr;
};

} // namespace gmx

#endif /* NBNXN_OPENCL_TYPES_H */
