/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

/*! \libinternal \file
 * \brief Defines the GPU-agnostic PME GPU data structures
 * (the host-side PME GPU data, and the GPU function parameters).
 * \todo Due to Gerrit workflow and time constraints, some renaming/refactoring
 * which does not impair the performance will be performed once
 * most of the initial PME CUDA implementation is merged
 * into the master branch (likely, after release 2017).
 * This should include:
 * -- bringing the structure/function names up to guidelines
 * ---- pme_gpu_settings_t -> PmeGpuTasks
 * -- refining GPU notation application (#2053)
 * -- renaming coefficients to charges (?)
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_TYPES_H
#define GMX_EWALD_PME_GPU_TYPES_H

#include "config.h"

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_hw_info;
struct gmx_device_info_t;

/*! \brief Possible PME codepaths
 * \todo: make this enum class with gmx_pme_t C++ refactoring
 */
enum PmeRunMode
{
    CPU,     //!< Whole PME step is done on CPU
    GPU,     //!< Whole PME step is done on GPU
    Hybrid,  //!< Mixed mode: only spread and gather run on GPU; FFT and solving are done on CPU.
};

//! PME gathering output forces treatment
enum class PmeForceOutputHandling
{
    Set,             /**< Gather simply writes into provided force buffer */
    ReduceWithInput, /**< Gather adds its output to the buffer.
                        On GPU, that means additional H2D copy before the kernel launch. */
};

#if GMX_GPU == GMX_GPU_CUDA

struct pme_gpu_cuda_t;
/*! \brief A typedef for including the GPU host data by pointer */
typedef pme_gpu_cuda_t pme_gpu_specific_t;

struct pme_gpu_cuda_kernel_params_t;
/*! \brief A typedef for including the GPU kernel arguments data by pointer */
typedef pme_gpu_cuda_kernel_params_t pme_gpu_kernel_params_t;

#else

/*! \brief A dummy typedef for the GPU host data placeholder on non-GPU builds */
typedef int pme_gpu_specific_t;
/*! \brief A dummy typedef for the GPU kernel arguments data placeholder on non-GPU builds */
typedef int pme_gpu_kernel_params_t;

#endif

/* What follows is all the PME GPU function arguments,
 * sorted into several device-side structures depending on the update rate.
 * This is GPU agnostic (float3 replaced by float[3], etc.).
 * The GPU-framework specifics (e.g. cudaTextureObject_t handles) are described
 * in the larger structure pme_gpu_cuda_kernel_params_t in the pme.cuh.
 */

/*! \internal \brief
 * A GPU data structure for storing the constant PME data.
 * This only has to be initialized once.
 */
struct pme_gpu_const_params_t
{
    /*! \brief Electrostatics coefficient = ONE_4PI_EPS0 / pme->epsilon_r */
    float elFactor;
    /*! \brief Virial and energy GPU array. Size is PME_GPU_ENERGY_AND_VIRIAL_COUNT (7) floats.
     * The element order is virxx, viryy, virzz, virxy, virxz, viryz, energy. */
    float *d_virialAndEnergy;
};

/*! \internal \brief
 * A GPU data structure for storing the PME data related to the grid sizes and cut-off.
 * This only has to be updated at every DD step.
 */
struct pme_gpu_grid_params_t
{
    /* Grid sizes */
    /*! \brief Real-space grid data dimensions. */
    int   realGridSize[DIM];
    /*! \brief Real-space grid dimensions, only converted to floating point. */
    float realGridSizeFP[DIM];
    /*! \brief Real-space grid dimensions (padded). The padding as compared to realGridSize includes the (order - 1) overlap. */
    int   realGridSizePadded[DIM]; /* Is major dimension of this ever used in kernels? */
    /*! \brief Fourier grid dimensions. This counts the complex numbers! */
    int   complexGridSize[DIM];
    /*! \brief Fourier grid dimensions (padded). This counts the complex numbers! */
    int   complexGridSizePadded[DIM];

    /* Grid pointers */
    /*! \brief Real space grid. */
    float *d_realGrid;
    /*! \brief Complex grid - used in FFT/solve. If inplace cuFFT is used, then it is the same pointer as realGrid. */
    float *d_fourierGrid;

    /*! \brief Ewald solving factor = (M_PI / pme->ewaldcoeff_q)^2 */
    float ewaldFactor;

    /*! \brief Grid spline values as in pme->bsp_mod
     * (laid out sequentially (XXX....XYYY......YZZZ.....Z))
     */
    float              *d_splineModuli;
    /*! \brief Offsets for X/Y/Z components of d_splineModuli */
    int                 splineValuesOffset[DIM];

    /*! \brief Fractional shifts lookup table as in pme->fshx/fshy/fshz, laid out sequentially (XXX....XYYY......YZZZ.....Z) */
    float               *d_fractShiftsTable;
    /*! \brief Gridline indices lookup table
     * (modulo lookup table as in pme->nnx/nny/nnz, laid out sequentially (XXX....XYYY......YZZZ.....Z)) */
    int                *d_gridlineIndicesTable;
    /*! \brief Offsets for X/Y/Z components of d_fractShiftsTable and d_gridlineIndicesTable */
    int                 tablesOffsets[DIM];
};

/*! \internal \brief
 * A GPU data structure for storing the PME data of the atoms, local to this process' domain partition.
 * This only has to be updated every DD step.
 */
struct pme_gpu_atom_params_t
{
    /*! \brief Number of local atoms */
    int    nAtoms;
    /*! \brief Pointer to the global GPU memory with input rvec atom coordinates.
     * The coordinates themselves change and need to be copied to the GPU every MD step,
     * but reallocation happens only at DD.
     */
    float *d_coordinates;
    /*! \brief Pointer to the global GPU memory with input atom charges.
     * The charges only need to be reallocated and copied to the GPU at DD step.
     */
    float  *d_coefficients;
    /*! \brief Pointer to the global GPU memory with input/output rvec atom forces.
     * The forces change and need to be copied from (and possibly to) the GPU every MD step,
     * but reallocation happens only at DD.
     */
    float  *d_forces;
    /*! \brief Pointer to the global GPU memory with ivec atom gridline indices.
     * Computed on GPU in the spline calculation part.
     */
    int *d_gridlineIndices;

    /* B-spline parameters are computed entirely on GPU every MD step, not copied.
     * Unless we want to try something like GPU spread + CPU gather?
     */
    /*! \brief Pointer to the global GPU memory with B-spline values */
    float  *d_theta;
    /*! \brief Pointer to the global GPU memory with B-spline derivative values */
    float  *d_dtheta;
};

/*! \internal \brief
 * A GPU data structure for storing the PME data which might change every MD step.
 */
struct pme_gpu_step_params_t
{
    /* The box parameters. The box only changes size each step with pressure coupling enabled. */
    /*! \brief
     * Reciprocal (inverted unit cell) box.
     *
     * The box is transposed as compared to the CPU pme->recipbox.
     * Basically, spread uses matrix columns (while solve and gather use rows).
     * This storage format might be not the most optimal since the box is always triangular so there are zeroes.
     */
    float  recipBox[DIM][DIM];
    /*! \brief The unit cell volume for solving. */
    float  boxVolume;
};

/*! \internal \brief
 * A single structure encompassing almost all the PME data used in GPU kernels on device.
 * This is inherited by the GPU framework-specific structure
 * (pme_gpu_cuda_kernel_params_t in pme.cuh).
 * This way, most code preparing the kernel parameters can be GPU-agnostic by casting
 * the kernel parameter data pointer to pme_gpu_kernel_params_base_t.
 */
struct pme_gpu_kernel_params_base_t
{
    /*! \brief Constant data that is set once. */
    pme_gpu_const_params_t constants;
    /*! \brief Data dependent on the grid size/cutoff. */
    pme_gpu_grid_params_t  grid;
    /*! \brief Data dependent on the DD and local atoms. */
    pme_gpu_atom_params_t  atoms;
    /*! \brief Data that possibly changes on every MD step. */
    pme_gpu_step_params_t  step;
};

/* Here are the host-side structures */

/*! \internal \brief
 * The PME GPU settings structure, included in the main PME GPU structure by value.
 */
struct pme_gpu_settings_t
{
    /* Permanent settings set on initialization */
    /*! \brief A boolean which tells if the solving is performed on GPU. Currently always true */
    bool performGPUSolve;
    /*! \brief A boolean which tells if the gathering is performed on GPU. Currently always true */
    bool performGPUGather;
    /*! \brief A boolean which tells if the FFT is performed on GPU. Currently true for a single MPI rank. */
    bool performGPUFFT;
    /*! \brief A convenience boolean which tells if PME decomposition is used. */
    bool useDecomposition;
    /*! \brief A boolean which tells if any PME GPU stage should copy all of its outputs to the host.
     * Only intended to be used by the test framework.
     */
    bool copyAllOutputs;
    /*! \brief Various computation flags for the current step, corresponding to the GMX_PME_ flags in pme.h. */
    int  stepFlags;
};

/*! \internal \brief
 * The PME GPU intermediate buffers structure, included in the main PME GPU structure by value.
 * Buffers are managed by the PME GPU module.
 */
struct pme_gpu_staging_t
{
    /*! \brief Virial and energy intermediate host-side buffer. Size is PME_GPU_VIRIAL_AND_ENERGY_COUNT. */
    float  *h_virialAndEnergy;
    /*! \brief B-spline values intermediate host-side buffer. */
    float  *h_splineModuli;

    /*! \brief Pointer to the host memory with B-spline values. Only used for host-side gather, or unit tests */
    float  *h_theta;
    /*! \brief Pointer to the host memory with B-spline derivative values. Only used for host-side gather, or unit tests */
    float  *h_dtheta;
    /*! \brief Pointer to the host memory with ivec atom gridline indices. Only used for host-side gather, or unit tests */
    int    *h_gridlineIndices;
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
struct pme_shared_t
{
    /*! \brief Grid count - currently always 1 on GPU */
    int ngrids;
    /*! \brief Grid dimensions - nkx, nky, nkz */
    int nk[DIM];
    /*! \brief Padded grid dimensions - pmegrid_nx, pmegrid_ny, pmegrid_nz
     * TODO: find out if these are really needed for the CPU FFT compatibility.
     */
    int                    pmegrid_n[DIM];
    /*! \brief PME interpolation order */
    int                    pme_order;
    /*! \brief Ewald splitting coefficient for Coulomb */
    real                   ewaldcoeff_q;
    /*! \brief Electrostatics parameter */
    real                   epsilon_r;
    /*! \brief Gridline indices - nnx, nny, nnz */
    std::vector<int>       nn;
    /*! \brief Fractional shifts - fshx, fshy, fshz */
    std::vector<real>      fsh;
    /*! \brief Precomputed B-spline values */
    std::vector<real>      bsp_mod[DIM];
    /*! \brief The PME codepath being taken */
    PmeRunMode             runMode;
    /*! \brief The box scaler based on inputrec - created in pme_init and managed by CPU structure */
    class EwaldBoxZScaler *boxScaler;
    /*! \brief The previous step box to know if we even need to update the current step box params.
     * \todo Manage this on higher level.
     * \todo Alternatively, when this structure is used by CPU PME code, make use of this field there as well.
     */
    matrix previousBox;
};

/*! \internal \brief
 * The main PME GPU host structure, included in the PME CPU structure by pointer.
 */
struct pme_gpu_t
{
    /*! \brief The information copied once per reinit from the CPU structure. */
    std::shared_ptr<pme_shared_t> common; // TODO: make the CPU structure use the same type

    /*! \brief The settings. */
    pme_gpu_settings_t settings;

    /*! \brief The host-side buffers.
     * The device-side buffers are buried in kernelParams, but that will have to change.
     */
    pme_gpu_staging_t staging;

    /*! \brief Number of local atoms, padded to be divisible by PME_ATOM_DATA_ALIGNMENT.
     * Used for kernel scheduling.
     * kernelParams.atoms.nAtoms is the actual atom count to be used for data copying.
     * TODO: this and the next member represent a memory allocation/padding properties -
     * what a container type should do ideally.
     */
    int nAtomsPadded;
    /*! \brief Number of local atoms, padded to be divisible by PME_ATOM_DATA_ALIGNMENT
     * if c_usePadding is true.
     * Used only as a basic size for almost all the atom data allocations
     * (spline parameter data is also aligned by PME_SPREADGATHER_PARTICLES_PER_WARP).
     * This should be the same as (c_usePadding ? nAtomsPadded : kernelParams.atoms.nAtoms).
     * kernelParams.atoms.nAtoms is the actual atom count to be used for most data copying.
     */
    int nAtomsAlloc;

    /*! \brief A pointer to the device used during the execution. */
    gmx_device_info_t *deviceInfo;

    /*! \brief A single structure encompassing all the PME data used on GPU.
     * Its value is the only argument to all the PME GPU kernels.
     * \todo Test whether this should be copied to the constant GPU memory once per MD step
     * (or even less often with no box updates) instead of being an argument.
     */
    std::shared_ptr<pme_gpu_kernel_params_t> kernelParams;

    /*! \brief The pointer to GPU-framework specific host-side data, such as CUDA streams and events. */
    std::shared_ptr<pme_gpu_specific_t> archSpecific; /* FIXME: make it an unique_ptr */
};

#endif
