/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 *  \brief Defines the GPU-agnostic PME GPU data structures:
 *  (the host-side PME GPU data, and the GPU function parameters).
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef PMEGPUTYPES_H
#define PMEGPUTYPES_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif


#if GMX_GPU == GMX_GPU_CUDA
struct pme_gpu_cuda_t;
/*! \brief  A typedef for including the GPU framework-specific data by pointer */
typedef pme_gpu_cuda_t pme_gpu_specific_t;
#else
/*! \brief A dummy typedef */
typedef int pme_gpu_specific_t;
#endif

/* What follows is all the PME GPU function arguments,
 * sorted into several device-side structures depending on the update rate.
 * This is almost entirely GPU agnostic (float3 replaced by float[3], etc.).
 * The only exception are 2 cudaTextureObject_t disguised as unsigned long long.
 */

/*! \brief \internal
 * A GPU data structure for storing the constant PME data.
 *
 * This only has to be initialized once.
 */
struct pme_gpu_const_params_t
{
    /*! \brief Electrostatics coefficient = ONE_4PI_EPS0 / pme->epsilon_r */
    float elFactor;
    /*! \brief Virial and energy GPU array. Size is PME_GPU_ENERGY_AND_VIRIAL_COUNT (7) floats.
     * The element order is virxx, viryy, virzz, virxy, virxz, viryz, energy. */
    float *virialAndEnergy;
};

/*! \brief \internal
 * A GPU data structure for storing the PME data related to the grid size and cut-off.
 *
 * This only has to be updated every DLB step.
 */
struct pme_gpu_grid_params_t
{
    /* Grid sizes */
    /*! \brief Grid data dimensions - integer. */
    int   localGridSize[DIM];
    /*! \brief Grid data dimensions - floating point. */
    float localGridSizeFP[DIM];
    /*! \brief Grid size dimensions - integer. The padding as compared to localGridSize includes the (order - 1) overlap. */
    int   localGridSizePadded[DIM]; /* Is major dimension of this ever used in kernels? */

    /* Grid pointers */
    /*! \brief Real space grid. */
    float  *realGrid;
    /*! \brief Complex grid - used in FFT/solve. If inplace cuFFT is used, then it is the same pointer as realGrid. */
    float  *fourierGrid;

    /*! \brief Count of the overlap zones */
#define PME_GPU_OVERLAP_ZONES_COUNT 7  /* can go away with a better rewrite of wrap/unwrap */
    /*! \brief The Y and Z sizes of the overlap zones */
    int  overlapSizes[2 * PME_GPU_OVERLAP_ZONES_COUNT];
    /*! \brief The total cell counts of the overlap zones */
    int  overlapCellCounts[PME_GPU_OVERLAP_ZONES_COUNT];

    /*! \brief Ewald solving factor = (M_PI / pme->ewaldcoeff_q)^2 */
    float ewaldFactor;

    /*! \brief Grid spline values as in pme->bsp_mod
     * (laid out sequentially (XXX....XYYY......YZZZ.....Z))
     */
    float              *splineValuesArray;
    /*! \brief Offsets for X/Y/Z components of splineValuesArray */
    int                 splineValuesOffset[DIM];

    /*! \brief Fractional shifts as in pme->fshx/fshy/fshz, laid out sequentially (XXX....XYYY......YZZZ.....Z) */
    float               *fshArray;
    /*! \brief Fractional shifts gridline indices
     * (modulo lookup table as in pme->nnx/nny/nnz, laid out sequentially (XXX....XYYY......YZZZ.....Z)) */
    int                *nnArray;
    /*! \brief Offsets for X/Y/Z components of fshArray and nnArray */
    int                 fshOffset[DIM];

    /* These are CUDA-specific kernel parameters.
     * Their actual type is cudaTextureObject_t (which is typedef'd as unsigned long long by CUDA itself).
     * Please don't use them outside of CUDA code.
     */
    /*! \brief Fractional shifts - a CUDA texture object for accessing fshArray. */
    unsigned long long fshTexture;
    /*! \brief Fractional shifts gridline indices - a CUDA texture object for accessing nnArray */
    unsigned long long nnTexture;
};

/*! \brief \internal
 * A GPU data structure for storing the PME data of the atoms, local to this process' domain partition.
 *
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
    float *coordinates;
    /*! \brief Pointer to the global GPU memory with input atom charges.
     * The charges only need to be reallocated and copied to the GPU on DD step.
     */
    float  *coefficients;
    /*! \brief Pointer to the global GPU memory with input/output rvec atom forces.
     * The forces change and need to be copied from (and possibly to) the GPU every MD step,
     * but reallocation happens only on DD.
     */
    float  *forces;
    /*! \brief Pointer to the global GPU memory with ivec atom gridline indices.
     * Computed on GPU in the spline calculation part.
     */
    int *gridlineIndices;

    /* B-spline parameters are computed entirely on GPU every MD step, not copied.
     * Unless we want to try something like GPU spread + CPU gather?
     */
    /*! \brief Pointer to the global GPU memory with B-spline values */
    float  *theta;
    /*! \brief Pointer to the global GPU memory with B-spline derivative values */
    float  *dtheta;
};

/*! \brief \internal
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

/*! \brief \internal
 * A single structure encompassing all the PME data used in GPU kernels.
 */
struct pme_gpu_kernel_params_t
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

/* Here is the host-side structures;
 */

/*! \brief \internal
 * The PME GPU settings structure, included in the main PME GPU structure by value.
 */
struct pme_gpu_settings_t
{
    /* Permanent settings set on initialization */
    /*! \brief A boolean which tells if the solving is performed on GPU. Currently always TRUE */
    gmx_bool bGPUSolve;
    /*! \brief A boolean which tells if the gathering is performed on GPU. Currently always TRUE */
    gmx_bool bGPUGather;
    /*! \brief A boolean which tells if the FFT is performed on GPU. Currently TRUE for a single MPI rank. */
    gmx_bool bGPUFFT;
    /*! \brief A convenience boolean which tells if there is only one PME GPU process. */
    gmx_bool bGPUSingle;
    /*! \brief A boolean which tells the PME to call pme_gpu_reinit_atoms() at the beginning of the run.
     * Set to TRUE initially, then to FALSE after the first MD step.
     * The pme_gpu_reinit_atoms() after the DD gets called directly in gmx_pmeonly.
     */
    gmx_bool bNeedToUpdateAtoms;
};

/*! \brief \internal
 * The PME GPU host-side I/O buffers structure, included in the main PME GPU structure by value.
 * Intermediate internal host buffers live here as well.
 * And what will happen with the introduction of the external device-side I/O pointers?
 */
struct pme_gpu_io_t
{
    /*! \brief Input coordinates (XYZ rvec) */
    float  *h_coordinates;
    /*! \brief Input charges */
    float  *h_coefficients;
    /*! \brief Output forces (and possibly the input if pme_kernel_gather does the reduction) */
    float  *h_forces;      /* rvec/float3 */
    /*! \brief Virial and energy intermediate host-side buffer, managed and pinned by PME GPU entirely. Size is PME_GPU_VIRIAL_AND_ENERGY_COUNT. */
    float  *h_virialAndEnergy;
    /*! \brief B-spline values intermediate host-side buffers, managed and pinned by PME GPU entirely. Sizes are the grid sizes. */
    float  *h_splineValues[DIM];
    /*! \brief Sizes of the corresponding h_splineValues arrays in bytes */
    size_t  splineValuesSizes[DIM];
};

/*! \brief \internal
 * The main PME GPU host structure, included in the PME CPU structure by pointer.
 */
struct pme_gpu_t
{
    /*! \brief The settings. */
    pme_gpu_settings_t settings;

    /*! \brief The host-side buffers.
     * The device-side buffers are buried in kernelParams, but that will have to change.
     */
    pme_gpu_io_t io;

    /*! \brief The unit cell box from the previous step.
     * Only used to know if the kernelParams.step needs to be updated.
     * Does not really fit anywhere else, does it?
     */
    matrix previousBox;

    /*! \brief A pointer to the device used during the execution. */
    struct gmx_device_info_t *deviceInfo;

    /*! \brief A single structure encompassing all the PME data used on GPU.
     * This should be the only parameter to all the PME GPU kernels (FIXME: pme_solve_kernel).
     * Can probably be copied to the constant GPU memory once per MD step (or even less often) instead of being a parameter.
     */
    pme_gpu_kernel_params_t kernelParams;

    /*! \brief The pointer to the GPU-framework specific host-side data, such as CUDA streams and events.*/
    pme_gpu_specific_t *archSpecific;
};

#ifdef __cplusplus
}
#endif

#endif // PMEGPUTYPES_H
