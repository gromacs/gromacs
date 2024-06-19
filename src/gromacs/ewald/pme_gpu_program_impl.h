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
/*! \internal \file
 * \brief
 * Declares PmeGpuProgramImpl, which stores PME GPU (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_PME_PME_GPU_PROGRAM_IMPL_H
#define GMX_EWALD_PME_PME_GPU_PROGRAM_IMPL_H

#include "config.h"

#include <cstddef>

#include <memory>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/utility/classhelpers.h"

class ISyclKernelFunctor;
class DeviceContext;
struct DeviceInformation;

/*! \internal
 * \brief
 * PME GPU persistent host program/kernel data, which should be initialized once for the whole execution.
 *
 * Primary purpose of this is to not recompile GPU kernels for each OpenCL unit test,
 * while the relevant GPU context (e.g. cl_context) instance persists.
 * In CUDA, this just assigns the kernel function pointers.
 * This also implicitly relies on the fact that reasonable share of the kernels are always used.
 * If there were more template parameters, even smaller share of all possible kernels would be used.
 *
 * \todo In future if we would need to react to either user input or
 * auto-tuning to compile different kernels, then we might wish to
 * revisit the number of kernels we pre-compile, and/or the management
 * of their lifetime.
 *
 * This also doesn't manage cuFFT/clFFT kernels, which depend on the PME grid dimensions.
 */
struct PmeGpuProgramImpl
{
    /*! \brief A handle to the device context (managed externally).
     *
     * Only used in OpenCL for JIT compilation.
     */
    const DeviceContext& deviceContext_;

    //! Conveniently all the PME kernels use the same single argument type
#if GMX_GPU_CUDA
    using PmeKernelHandle = void (*)(const struct PmeGpuCudaKernelParams);
#elif GMX_GPU_OPENCL
    using PmeKernelHandle = cl_kernel;
#else
    using PmeKernelHandle = ISyclKernelFunctor*;
#endif

    /*! \brief
     * Maximum synchronous GPU thread group execution width.
     * "Warp" is a CUDA term which we end up reusing in OpenCL kernels as well.
     * For CUDA, this is a static value that comes from gromacs/gpu_utils/cuda_arch_utils.cuh;
     * for OpenCL, we have to query it dynamically.
     */
    size_t warpSize_;

    //@{
    /**
     * Spread/spline kernels are compiled only for order of 4.
     * There are multiple versions of each kernel, paramaretized according to
     *   Number of threads per atom. Using either order(4) or order*order (16) threads per atom is
     * supported If the spline data is written in the spline/spread kernel and loaded in the gather
     *   or recalculated in the gather.
     * Spreading kernels also have hardcoded X/Y indices wrapping parameters,
     * as a placeholder for implementing 1/2D decomposition.
     * The kernels are templated separately for spreading on one grid (one or
     * two sets of coefficients) or on two grids (required for energy and virial
     * calculations).
     */
    size_t spreadWorkGroupSize;

    PmeKernelHandle splineKernelSingle;
    PmeKernelHandle splineKernelThPerAtom4Single;
    PmeKernelHandle spreadKernelSingle;
    PmeKernelHandle spreadKernelThPerAtom4Single;
    PmeKernelHandle splineAndSpreadKernelSingle;
    PmeKernelHandle splineAndSpreadKernelThPerAtom4Single;
    PmeKernelHandle splineAndSpreadKernelWriteSplinesSingle;
    PmeKernelHandle splineAndSpreadKernelWriteSplinesThPerAtom4Single;
    PmeKernelHandle splineKernelDual;
    PmeKernelHandle splineKernelThPerAtom4Dual;
    PmeKernelHandle spreadKernelDual;
    PmeKernelHandle spreadKernelThPerAtom4Dual;
    PmeKernelHandle splineAndSpreadKernelDual;
    PmeKernelHandle splineAndSpreadKernelThPerAtom4Dual;
    PmeKernelHandle splineAndSpreadKernelWriteSplinesDual;
    PmeKernelHandle splineAndSpreadKernelWriteSplinesThPerAtom4Dual;
    //@}

    //@{
    /** Same for gather: hardcoded X/Y unwrap parameters, order of 4, plus
     * it can either reduce with previous forces in the host buffer, or ignore them.
     * Also similarly to the gather we can use either order(4) or order*order (16) threads per atom
     * and either recalculate the splines or read the ones written by the spread
     * The kernels are templated separately for using one or two grids (required for
     * calculating energies and virial).
     */
    size_t gatherWorkGroupSize;

    PmeKernelHandle gatherKernelSingle;
    PmeKernelHandle gatherKernelThPerAtom4Single;
    PmeKernelHandle gatherKernelReadSplinesSingle;
    PmeKernelHandle gatherKernelReadSplinesThPerAtom4Single;
    PmeKernelHandle gatherKernelDual;
    PmeKernelHandle gatherKernelThPerAtom4Dual;
    PmeKernelHandle gatherKernelReadSplinesDual;
    PmeKernelHandle gatherKernelReadSplinesThPerAtom4Dual;
    //@}

    //@{
    /** Solve kernel doesn't care about the interpolation order, but can optionally
     * compute energy and virial, and supports XYZ and YZX grid orderings.
     * The kernels are templated separately for grids in state A and B.
     */
    size_t solveMaxWorkGroupSize;

    PmeKernelHandle solveYZXKernelA;
    PmeKernelHandle solveXYZKernelA;
    PmeKernelHandle solveYZXEnergyKernelA;
    PmeKernelHandle solveXYZEnergyKernelA;
    PmeKernelHandle solveYZXKernelB;
    PmeKernelHandle solveXYZKernelB;
    PmeKernelHandle solveYZXEnergyKernelB;
    PmeKernelHandle solveXYZEnergyKernelB;
    //@}
    PmeKernelHandle nvshmemSignalKern;

    PmeGpuProgramImpl() = delete;
    //! Constructor for the given device
    explicit PmeGpuProgramImpl(const DeviceContext& deviceContext);
    // NOLINTNEXTLINE(performance-trivially-destructible)
    ~PmeGpuProgramImpl();
    GMX_DISALLOW_COPY_AND_ASSIGN(PmeGpuProgramImpl);

    //! Return the warp size for which the kernels were compiled
    int warpSize() const { return warpSize_; }

private:
    // Compiles kernels, if supported. Called by the constructor.
    void compileKernels(const DeviceInformation& deviceInfo);
};

#endif
