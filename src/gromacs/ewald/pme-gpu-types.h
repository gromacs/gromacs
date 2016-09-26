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
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef PMEGPUTYPES_H
#define PMEGPUTYPES_H

#include "config.h"

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_hw_info;
struct gmx_device_info_t;

#if GMX_GPU == GMX_GPU_CUDA

struct pme_gpu_cuda_t;
/*! \brief A typedef for including the GPU framework-specific host data by pointer */
typedef pme_gpu_cuda_t pme_gpu_specific_t;

struct pme_gpu_cuda_kernel_params_t;
/*! \brief A typedef for including the GPU framework-specific kernel arguments data by pointer */
typedef pme_gpu_cuda_kernel_params_t pme_gpu_kernel_params_t;

#else

/*! \brief A dummy typedef */
typedef int pme_gpu_specific_t;
/*! \brief A dummy typedef */
typedef int pme_gpu_kernel_params_t;

#endif

/* What follows is all the PME GPU function arguments,
 * sorted into several device-side structures depending on the update rate.
 * This is almost entirely GPU agnostic (float3 replaced by float[3], etc.).
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
 * This only has to be updated at every DLB step.
 */
struct pme_gpu_grid_params_t
{
    /* Grid sizes */
    /*! \brief Real-space grid data dimensions. */
    int   realGridSize[DIM];
    /*! \brief Real-space grid dimensions - floating point. */
    float realGridSizeFP[DIM];
    /*! \brief Real-space grid dimensions (padded). The padding as compared to realGridSize includes the (order - 1) overlap. */
    int   realGridSizePadded[DIM]; /* Is major dimension of this ever used in kernels? */
    /*! \brief Fourier grid dimensions.*/
    int   complexGridSize[DIM];
    /*! \brief Fourier grid dimensions (padded). */
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
     * but reallocation happens only on DD.
     */
    float *d_coordinates;
    /*! \brief Pointer to the global GPU memory with input atom charges.
     * The charges only need to be reallocated and copied to the GPU on DD step.
     */
    float  *d_coefficients;
    /*! \brief Pointer to the global GPU memory with input/output rvec atom forces.
     * The forces change and need to be copied from (and possibly to) the GPU every MD step,
     * but reallocation happens only on DD.
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
    /*! \brief A convenience boolean which tells if there is only one PME GPU process. */
    bool useDecomposition;
    /*! \brief A boolean which tells the PME GPU to call pme_gpu_reinit_atoms() at the beginning of the run.
     * Set to true initially, then to false after the first MD step.
     * The pme_reinit_atoms() after the DD gets called directly in gmx_pmeonly.
     */
    bool needToUpdateAtoms;
};

/*! \internal \brief
 * The PME GPU intermediate buffers structure, included in the main PME GPU structure by value.
 * Buffers are managed by the PME GPU.
 */
struct pme_gpu_staging_t
{
    /*! \brief Virial and energy intermediate host-side buffer. Size is PME_GPU_VIRIAL_AND_ENERGY_COUNT. */
    float  *h_virialAndEnergy;
    /*! \brief B-spline values intermediate host-side buffers. Sizes are the grid sizes. */
    float  *h_splineModuli[DIM];
    /*! \brief Sizes of the corresponding h_splineValues arrays in bytes */
    size_t  splineModuliSizes[DIM];

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
    int               pmegrid_n[DIM];
    /*! \brief PME interpolation order */
    int               pme_order;
    /*! \brief Ewald splitting coefficient for Coulomb */
    real              ewaldcoeff_q;
    /*! \brief Electrostatics parameter */
    real              epsilon_r;
    /*! \brief Gridline indices - nnx, nny, nnz */
    std::vector<int>  nn[DIM];
    /*! \brief Fractional shifts - fshx, fshy, fshz */
    std::vector<real> fsh[DIM];
    /*! \brief Precomputed B-spline values */
    std::vector<real> bsp_mod[DIM];
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

    /*! \brief The unit cell box from the previous step.
     * Only used to know if the kernelParams.step needs to be updated.
     * \todo Figure out if this can go away with the once-per-step constant memory update.
     */
    matrix previousBox;
    /*! \brief Number of local atoms, padded to be divisible by PME_SPREADGATHER_ATOMS_PER_BLOCK.
     * Used for kernel scheduling.
     * kernelParams.atoms.nAtoms is the actual atom count to be used for data copying.
     */
    int nAtomsPadded;
    /*! \brief Number of local atoms, padded to be divisible by PME_SPREADGATHER_ATOMS_PER_BLOCK
     * if (PME_GPU_USE_PADDING == 1).
     * Used only as a basic size for almost all the atom data allocations
     * (spline parameter data is also aligned by PME_SPREADGATHER_PARTICLES_PER_WARP).
     * This should be the same as (PME_GPU_USE_PADDING ? nAtomsPadded : kernelParams.atoms.nAtoms).
     * kernelParams.atoms.nAtoms is the actual atom count to be used for data copying.
     */
    int nAtomsAlloc;

    /*! \brief A pointer to the device used during the execution. */
    gmx_device_info_t *deviceInfo;

    /*! \brief A single structure encompassing all the PME data used on GPU.
     * This should be the only parameter to all the PME GPU kernels (FIXME: pme_solve_kernel).
     * Can probably be copied to the constant GPU memory once per MD step (or even less often) instead of being a parameter.
     */
    std::shared_ptr<pme_gpu_kernel_params_t> kernelParams;

    /*! \brief The pointer to GPU-framework specific host-side data, such as CUDA streams and events. */
    std::shared_ptr<pme_gpu_specific_t> archSpecific; /* FIXME: make it an unique_ptr */
};

#endif // PMEGPUTYPES_H
