/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *  \brief Declare utility routines for SYCL
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_SYCLUTILS_H
#define GMX_GPU_UTILS_SYCLUTILS_H

#include <string>

#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

class DeviceStream;

#ifndef DOXYGEN

//! \brief Interface for SYCL kernel function objects.
class ISyclKernelFunctor
{
public:
    //! \brief Virtual destructor.
    virtual ~ISyclKernelFunctor() = default;
    /*! \brief Set the kernel argument number \p argIndex to \p arg.
     *
     * \param argIndex Index of the argument. Maximum allowed value depends
     *                 on the specific concrete class implementing this interface.
     * \param arg      Pointer to the argument value.
     *
     * \note Valid values of \p argIndex and types of \p arg depend on the
     *       specific concrete class implementing this interface. Passing
     *       illegal values is undefined behavior.
     * \note Similar to \c clSetKernelArg, it is not safe to call this
     *       function on the same kernel object from multiple host threads.
     */
    virtual void setArg(size_t argIndex, void* arg) = 0;
    /*! \brief Launch the kernel.
     *
     * \param config       Work-group configuration.
     * \param deviceStream \c DeviceStream to use.
     */
    virtual void launch(const KernelLaunchConfig& /*config*/, const DeviceStream& /*deviceStream*/) = 0;
};


/*! \brief
 * A function for setting up a single SYCL kernel argument.
 * This is the tail of the compile-time recursive function below.
 * It has to be seen by the compiler first.
 *
 * \param[in]     kernel          Kernel function handle
 * \param[in]     argIndex        Index of the current argument
 */
void inline prepareGpuKernelArgument(ISyclKernelFunctor* /*kernel*/, size_t /*argIndex*/) {}

/*! \brief
 * Compile-time recursive function for setting up a single SYCL kernel argument.
 * This function uses one kernel argument pointer \p argPtr to call
 * \c ISyclKernelFunctor::setArg, and calls itself on the next argument, eventually
 * calling the tail function above.
 *
 * \tparam        CurrentArg      Type of the current argument
 * \tparam        RemainingArgs   Types of remaining arguments after the current one
 * \param[in]     kernel          Kernel function handle
 * \param[in]     argIndex        Index of the current argument
 * \param[in]     argPtr          Pointer to the current argument
 * \param[in]     otherArgsPtrs   Pack of pointers to arguments remaining to process after the current one
 */
template<typename CurrentArg, typename... RemainingArgs>
void prepareGpuKernelArgument(ISyclKernelFunctor* kernel,
                              size_t              argIndex,
                              const CurrentArg*   argPtr,
                              const RemainingArgs*... otherArgsPtrs)
{
    kernel->setArg(argIndex, const_cast<void*>(reinterpret_cast<const void*>(argPtr)));
    prepareGpuKernelArgument(kernel, argIndex + 1, otherArgsPtrs...);
}

/*! \brief
 * A wrapper function for setting up all the SYCL kernel arguments.
 * Calls the recursive functions above.
 *
 * \tparam    Args            Types of all the kernel arguments
 * \param[in] kernel          Kernel function handle
 * \param[in] config          Kernel configuration for launching
 * \param[in] argsPtrs        Pointers to all the kernel arguments
 * \returns A dummy value to be used with launchGpuKernel() as the last argument.
 */
template<typename... Args>
void* prepareGpuKernelArguments(void* kernel, const KernelLaunchConfig& /*config*/, const Args*... argsPtrs)
{
    auto* kernelFunctor = reinterpret_cast<ISyclKernelFunctor*>(kernel);
    prepareGpuKernelArgument(kernelFunctor, 0, argsPtrs...);
    return nullptr;
}

/*! \brief Launches the SYCL kernel and handles the errors.
 *
 * \param[in] kernel          Kernel function handle
 * \param[in] config          Kernel configuration for launching
 * \param[in] deviceStream    GPU stream to launch kernel in
 * \param[in] timingEvent     Timing event, fetched from GpuRegionTimer. Unused.
 * \param[in] kernelName      Human readable kernel description, for error handling only. Unused.
 * \param[in] kernelArgs      Unused.
 * \throws gmx::InternalError on kernel launch failure
 */
inline void launchGpuKernel(void*                     kernel,
                            const KernelLaunchConfig& config,
                            const DeviceStream&       deviceStream,
                            CommandEvent* /*timingEvent*/,
                            const char* /*kernelName*/,
                            const void* /*kernelArgs*/)
{
    auto* kernelFunctor = reinterpret_cast<ISyclKernelFunctor*>(kernel);
    kernelFunctor->launch(config, deviceStream);
}

/* To properly mark function as [[noreturn]], we must do it everywhere it is declared, which
 * will pollute common headers.*/
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")

/*! \brief Pretend to check a SYCL stream for unfinished work (dummy implementation).
 *
 *  \returns  Not implemented in SYCL.
 */
static inline bool haveStreamTasksCompleted(const DeviceStream& /* deviceStream */)
{
    GMX_THROW(gmx::NotImplementedError("Not implemented on SYCL yet"));
}

CLANG_DIAGNOSTIC_RESET

#endif // !DOXYGEN

#endif
