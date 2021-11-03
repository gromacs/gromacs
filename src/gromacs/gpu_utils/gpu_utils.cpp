/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017,2018,2019, by the GROMACS development team.
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 *  \brief Function definitions for non-GPU builds
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gpu_utils.h"

#include "config.h"

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/stringutil.h"

#ifdef _MSC_VER
#    pragma warning(disable : 6237)
#endif

const char* enumValueToString(GpuApiCallBehavior enumValue)
{
    static constexpr gmx::EnumerationArray<GpuApiCallBehavior, const char*> s_gpuApiCallBehaviorNames = {
        "Synchronous", "Asynchronous"
    };
    return s_gpuApiCallBehaviorNames[enumValue];
}

bool decideGpuTimingsUsage()
{
    if (GMX_GPU_CUDA || GMX_GPU_SYCL)
    {
        /* CUDA: timings are incorrect with multiple streams.
         * This is the main reason why they are disabled by default.
         * TODO: Consider turning on by default when we can detect nr of streams.
         *
         * SYCL: compilers and runtimes change rapidly, so we disable timings by default
         * to avoid any possible overhead. */
        return (getenv("GMX_ENABLE_GPU_TIMING") != nullptr);
    }
    else if (GMX_GPU_OPENCL)
    {
        return (getenv("GMX_DISABLE_GPU_TIMING") == nullptr);
    }
    else
    {
        // CPU-only build
        return false;
    }
}
