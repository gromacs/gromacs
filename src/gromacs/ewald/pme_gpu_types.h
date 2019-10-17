/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief Defines the PME GPU data structures
 * (the GPU function parameters used both on host and device sides).
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_TYPES_H
#define GMX_EWALD_PME_GPU_TYPES_H

/*
 * In OpenCL, the structures must be laid out on the host and device exactly the same way.
 * If something is off, one might get an error CL_INVALID_ARG_SIZE if any structure's sizes don't
 * match. What's worse, structures might be of same size but members might be aligned differently,
 * resulting in wrong kernel results. The structures below are aligned manually.
 * The pattern is ordering the members of structs from smallest to largest sizeof
 * (arrays behave the same way as sequences of separate fields),
 * as described in "The Lost Art of C Structure Packing".
 *
 * However, if the need arises at some point, they can all be aligned forcefully:
 *
 * #define GMX_GPU_ALIGNED __attribute__ ((aligned(8)))
 * struct GMX_GPU_ALIGNED PmeGpuConstParams
 * struct GMX_GPU_ALIGNED PmeGpuGridParams
 * etc...
 *
 * One might also try __attribute__ ((packed)), but it doesn't work with DeviceBuffer,
 * as it appears to not be POD.
 */


/*! \brief A workaround to hide DeviceBuffer template from OpenCL kernel compilation
 * - to turn it into a dummy of the same size as host implementation of device buffer.
 * As we only care about 64-bit, 8 bytes is fine.
 * TODO: what we should be doing is providing separate device-side views of the same structures -
 * then there would be no need for macro.
 */
#ifndef __OPENCL_C_VERSION__
#    include "gromacs/gpu_utils/devicebuffer.h"
#    define HIDE_FROM_OPENCL_COMPILER(x) x
static_assert(sizeof(DeviceBuffer<float>) == 8,
              "DeviceBuffer is defined as an 8 byte stub for OpenCL C");
static_assert(sizeof(DeviceBuffer<int>) == 8,
              "DeviceBuffer is defined as an 8 byte stub for OpenCL C");
#else
#    define HIDE_FROM_OPENCL_COMPILER(x) char8
#endif

/* What follows is all the PME GPU function arguments,
 * sorted into several device-side structures depending on the update rate.
 * This is GPU agnostic (float3 replaced by float[3], etc.).
 * The GPU-framework specifics (e.g. cudaTextureObject_t handles) are described
 * in the larger structure PmeGpuCudaKernelParams in the pme.cuh.
 */

/*! \internal \brief
 * A GPU data structure for storing the constant PME data.
 * This only has to be initialized once.
 */
struct PmeGpuConstParams
{
    /*! \brief Electrostatics coefficient = ONE_4PI_EPS0 / pme->epsilon_r */
    float elFactor;
    /*! \brief Virial and energy GPU array. Size is PME_GPU_ENERGY_AND_VIRIAL_COUNT (7) floats.
     * The element order is virxx, viryy, virzz, virxy, virxz, viryz, energy. */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_virialAndEnergy;
};

/*! \internal \brief
 * A GPU data structure for storing the PME data related to the grid sizes and cut-off.
 * This only has to be updated at every DD step.
 */
struct PmeGpuGridParams
{
    /*! \brief Ewald solving factor = (M_PI / pme->ewaldcoeff_q)^2 */
    float ewaldFactor;

    /* Grid sizes */
    /*! \brief Real-space grid data dimensions. */
    int realGridSize[DIM];
    /*! \brief Real-space grid dimensions, only converted to floating point. */
    float realGridSizeFP[DIM];
    /*! \brief Real-space grid dimensions (padded). The padding as compared to realGridSize includes the (order - 1) overlap. */
    int realGridSizePadded[DIM]; /* Is major dimension of this ever used in kernels? */
    /*! \brief Fourier grid dimensions. This counts the complex numbers! */
    int complexGridSize[DIM];
    /*! \brief Fourier grid dimensions (padded). This counts the complex numbers! */
    int complexGridSizePadded[DIM];

    /*! \brief Offsets for X/Y/Z components of d_splineModuli */
    int splineValuesOffset[DIM];
    /*! \brief Offsets for X/Y/Z components of d_fractShiftsTable and d_gridlineIndicesTable */
    int tablesOffsets[DIM];

    /* Grid arrays */
    /*! \brief Real space grid. */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_realGrid;
    /*! \brief Complex grid - used in FFT/solve. If inplace cu/clFFT is used, then it is the same handle as realGrid. */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_fourierGrid;

    /*! \brief Grid spline values as in pme->bsp_mod
     * (laid out sequentially (XXX....XYYY......YZZZ.....Z))
     */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_splineModuli;
    /*! \brief Fractional shifts lookup table as in pme->fshx/fshy/fshz, laid out sequentially (XXX....XYYY......YZZZ.....Z) */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_fractShiftsTable;
    /*! \brief Gridline indices lookup table
     * (modulo lookup table as in pme->nnx/nny/nnz, laid out sequentially (XXX....XYYY......YZZZ.....Z)) */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<int>) d_gridlineIndicesTable;
};

/*! \internal \brief
 * A GPU data structure for storing the PME data of the atoms, local to this process' domain
 * partition. This only has to be updated every DD step.
 */
struct PmeGpuAtomParams
{
    /*! \brief Number of local atoms */
    int nAtoms;
    /*! \brief Global GPU memory array handle with input rvec atom coordinates.
     * The coordinates themselves change and need to be copied to the GPU for every PME computation,
     * but reallocation happens only at DD.
     */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_coordinates;
    /*! \brief Global GPU memory array handle with input atom charges.
     * The charges only need to be reallocated and copied to the GPU at DD step.
     */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_coefficients;
    /*! \brief Global GPU memory array handle with input/output rvec atom forces.
     * The forces change and need to be copied from (and possibly to) the GPU for every PME
     * computation, but reallocation happens only at DD.
     */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_forces;
    /*! \brief Global GPU memory array handle with ivec atom gridline indices.
     * Computed on GPU in the spline calculation part.
     */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<int>) d_gridlineIndices;
    /* B-spline parameters are computed entirely on GPU for every PME computation, not copied.
     * Unless we want to try something like GPU spread + CPU gather?
     */
    /*! \brief Global GPU memory array handle with B-spline values */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_theta;
    /*! \brief Global GPU memory array handle with B-spline derivative values */
    HIDE_FROM_OPENCL_COMPILER(DeviceBuffer<float>) d_dtheta;
};

/*! \internal \brief
 * A GPU data structure for storing the PME data which might change for each new PME computation.
 */
struct PmeGpuDynamicParams
{
    /* The box parameters. The box only changes size with pressure coupling enabled. */
    /*! \brief
     * Reciprocal (inverted unit cell) box.
     *
     * The box is transposed as compared to the CPU pme->recipbox.
     * Basically, spread uses matrix columns (while solve and gather use rows).
     * This storage format might be not the most optimal since the box is always triangular so there are zeroes.
     */
    float recipBox[DIM][DIM];
    /*! \brief The unit cell volume for solving. */
    float boxVolume;
};

/*! \internal \brief
 * A single structure encompassing almost all the PME data used in GPU kernels on device.
 * This is inherited by the GPU framework-specific structure
 * (PmeGpuCudaKernelParams in pme.cuh).
 * This way, most code preparing the kernel parameters can be GPU-agnostic by casting
 * the kernel parameter data pointer to PmeGpuKernelParamsBase.
 */
struct PmeGpuKernelParamsBase
{
    /*! \brief Constant data that is set once. */
    struct PmeGpuConstParams constants;
    /*! \brief Data dependent on the grid size/cutoff. */
    struct PmeGpuGridParams grid;
    /*! \brief Data dependent on the DD and local atoms. */
    struct PmeGpuAtomParams atoms;
    /*! \brief Data that possibly changes for every new PME computation.
     * This should be kept up-to-date by calling pme_gpu_prepare_computation(...)
     * before launching spreading.
     */
    struct PmeGpuDynamicParams current;
};

#endif
