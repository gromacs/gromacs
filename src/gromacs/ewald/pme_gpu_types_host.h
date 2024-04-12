/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Defines the host-side PME GPU data structures.
 * \todo Some renaming/refactoring, which does not impair the performance:
 * -- bringing the function names up to guidelines
 * -- PmeGpuSettings -> PmeGpuTasks
 * -- refining GPU notation application (#2053)
 * -- renaming coefficients to charges (?)
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_TYPES_HOST_H
#define GMX_EWALD_PME_GPU_TYPES_HOST_H

#include "config.h"

#include <memory>
#include <vector>

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_gpu_program.h"
#include "gromacs/gpu_utils/clfftinitializer.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"

#include "pme_force_sender_gpu.h"
#include "pme_gpu_settings.h"
#include "pme_gpu_staging.h"

#if GMX_GPU && !GMX_GPU_HIP
struct PmeGpuSpecific;
struct PmeGpuHaloExchange;
#else
/*! \brief A dummy typedef for the GPU host data placeholder on non-GPU builds */
typedef int PmeGpuSpecific;
typedef int PmeGpuHaloExchange;
#endif

#if GMX_GPU_CUDA
struct PmeGpuCudaKernelParams;
/*! \brief A typedef for including the GPU kernel arguments data by pointer */
typedef PmeGpuCudaKernelParams PmeGpuKernelParams;
#elif GMX_GPU_OPENCL || GMX_GPU_SYCL
struct PmeGpuKernelParamsBase;
/*! \brief A typedef for including the GPU kernel arguments data by pointer */
typedef PmeGpuKernelParamsBase PmeGpuKernelParams;
#else
/*! \brief A dummy typedef for the GPU kernel arguments data placeholder on non-GPU builds */
typedef int PmeGpuKernelParams;
#endif

struct DeviceInformation;

/*! \internal \brief
 * Contains information about the PP ranks that partners this PME rank.
 * used in the pme gather kernel for nvshmem purpose*/
struct PpRanksSendFInfo
{
    //! The MPI rank ID of this partner PP rank.
    int rankId = -1;
    //! The number of atoms to communicate with this partner PP rank.
    int numAtoms = -1;
    //! starting offset of this PP rank in the larger force buffer
    int startAtomOffset = -1;
};

/*! \internal \brief
 * The PME GPU structure for all the data copied directly from the CPU PME structure.
 * The copying is done when the CPU PME structure is already (re-)initialized
 * (pme_gpu_reinit is called at the end of gmx_pme_init).
 * All the variables here are named almost the same way as in gmx_pme_t.
 * The types are different: pointers are replaced by vectors.
 * TODO: use the shared data with the PME CPU.
 * Included in the main PME GPU structure by value.
 */
struct PmeShared
{
    /*! \brief Grid count */
    int ngrids;
    /*! \brief Grid dimensions - nkx, nky, nkz */
    int nk[DIM];
    /*! \brief PME interpolation order */
    int pme_order;
    /*! \brief Ewald splitting coefficient for Coulomb */
    real ewaldcoeff_q;
    /*! \brief Electrostatics parameter */
    real epsilon_r;
    /*! \brief Gridline indices - nnx, nny, nnz */
    std::vector<int> nn;
    /*! \brief Fractional shifts - fshx, fshy, fshz */
    std::vector<real> fsh;
    /*! \brief Precomputed B-spline values */
    std::vector<real> bsp_mod[DIM];
    /*! \brief The PME codepath being taken */
    PmeRunMode runMode;
    /*! \brief  Whether PME execution is happening on a PME-only rank (from gmx_pme_t.bPPnode). */
    bool isRankPmeOnly;
    /*! \brief The box scaler based on inputrec - created in pme_init and managed by CPU structure */
    class EwaldBoxZScaler* boxScaler;
    /*! \brief The previous computation box to know if we even need to update the current box params.
     * \todo Manage this on higher level.
     * \todo Alternatively, when this structure is used by CPU PME code, make use of this field there as well.
     */
    matrix previousBox;

    /*! \brief The The number of decomposition dimensions */
    int ndecompdim;
    /*! \brief MPI rank within communicator */
    int nodeid;
    /*! \brief Number of MPI ranks doing PME */
    int nnodes;
    /*! \brief MPI rank within communicator for PME X-decomposition */
    int nodeidX;
    /*! \brief MPI rank within communicator for PME Y-decomposition */
    int nodeidY;
    /*! \brief Number of MPI ranks in X-decomposition */
    int nnodesX;
    /*! \brief Number of MPI ranks in Y-decomposition */
    int nnodesY;
    /*! \brief MPI communicator for PME ranks */
    MPI_Comm mpiComm;
    /*! \brief MPI communicator for ranks in X-decomposition */
    MPI_Comm mpiCommX;
    /*! \brief MPI communicator for ranks in Y-decomposition */
    MPI_Comm mpiCommY;
    /*! \brief local interpolation grid start values in x-dimension*/
    std::vector<int> s2g0X;
    /*! \brief local interpolation grid end values in x-dimension*/
    std::vector<int> s2g1X;
    /*! \brief local interpolation grid start values in y-dimension*/
    std::vector<int> s2g0Y;
    /*! \brief local interpolation grid end values in y-dimension*/
    std::vector<int> s2g1Y;
    /*! \brief local grid size*/
    std::array<int, DIM> pmegridNk;
    /*! \brief Size of the grid halo region */
    int gridHalo;
};

struct PmeNvshmemHost
{
    /*! \brief symmetric value across PME+PP ranks of max num atoms to allocate */
    int nAtomsAlloc_symmetric;
    /*! \brief PpRanks struct reference */
    gmx::ArrayRef<PpRanks> ppRanksRef;
    /*! \brief PpRanksSendFInfo struct containing info about each PP rank offsets */
    gmx::HostVector<PpRanksSendFInfo> ppRanksFInfo;
    /*! \brief PpRanksSendFInfo struct allocation size tracker */
    int ppRanksFInfoSize = 0;
    /*! \brief PpRanksSendFInfo struct allocation size tracker */
    int ppRanksFInfoSizeAlloc = 0;
    /*! \brief per pp rank atomic counter allocation size tracker */
    int lastProcessedBlockPerPpRankSize = 0;
    /*! \brief per pp rank atomic counter allocation size tracker */
    int lastProcessedBlockPerPpRankSizeAlloc = 0;
    /*! \brief sync object for nvshmem based pme-pp force comm allocation size tracker */
    int forcesReadyNvshmemFlagsSize = 0;
    /*! \brief sync object for nvshmem based pme-pp force comm allocation size tracker */
    int forcesReadyNvshmemFlagsSizeAlloc = 0;
};

/*! \internal \brief
 * The main PME GPU host structure, included in the PME CPU structure by pointer.
 */
struct PmeGpu
{
    /*! \brief The information copied once per reinit from the CPU structure. */
    std::shared_ptr<PmeShared> common; // TODO: make the CPU structure use the same type

    //! A handle to the program created by buildPmeGpuProgram()
    const PmeGpuProgram* programHandle_;

    //! Handle that ensures the clFFT library has been initialized once per process.
    std::unique_ptr<gmx::ClfftInitializer> initializedClfftLibrary_;

    /*! \brief The settings. */
    PmeGpuSettings settings;

    /*! \brief The host-side buffers.
     * The device-side buffers are buried in kernelParams, but that will have to change.
     */
    PmeGpuStaging staging;

    /*! \brief Number of local atoms, padded to be divisible by c_pmeAtomDataAlignment.
     *
     * Used only as a basic size for almost all the atom data allocations
     * (spline parameter data is also aligned by PME_SPREADGATHER_PARTICLES_PER_WARP).
     * kernelParams.atoms.nAtoms is the actual atom count to be used for most data copying.
     *
     * TODO: memory allocation/padding properties should be handled by
     * something like a container
     */
    int nAtomsAlloc;

    /*! \brief Kernel scheduling grid width limit in X - derived from deviceinfo compute capability in CUDA.
     * Declared as very large int to make it useful in computations with type promotion, to avoid overflows.
     * OpenCL seems to not have readily available global work size limit, so we just assign a large arbitrary constant to this instead.
     * TODO: this should be in PmeGpuProgram(Impl)
     */
    std::intmax_t maxGridWidthX;

    /*! \brief Minimum particle count to prefer recalculating splines.
     *
     * The gather kernel can either recalculate the splines or load
     * those saved during the spline (and spread)
     * kernel. Recalculating is advantageous when there are enough
     * particles. When so doing, it is best to use fewer threads per
     * atom in the spline and spread.
     *
     * This feature is supported by CUDA and SYCL. Spread pipelining
     * requires spline recalculation. */
    int minParticleCountToRecalculateSplines = 23000;

    /*! \brief A single structure encompassing all the PME data used on GPU.
     * Its value is the only argument to all the PME GPU kernels.
     * \todo Test whether this should be copied to the constant GPU memory once for each computation
     * (or even less often with no box updates) instead of being an argument.
     */
    std::shared_ptr<PmeGpuKernelParams> kernelParams;

    /*! \brief The pointer to GPU-framework specific host-side data, such as CUDA streams and events. */
    std::shared_ptr<PmeGpuSpecific> archSpecific; /* FIXME: make it an unique_ptr */

    /*! \brief The pointer to PME halo-exchange specific host-side data */
    std::unique_ptr<PmeGpuHaloExchange> haloExchange;

    bool useNvshmem = false;

    std::unique_ptr<PmeNvshmemHost> nvshmemParams;
};

#endif
