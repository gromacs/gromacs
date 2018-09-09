/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief
 * Implements PmeGpuProgramImpl, which stores permanent PME GPU context-derived data,
 * such as (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/ocl_compiler.h"
#include "gromacs/utility/stringutil.h"

#include "pme-gpu-internal.h" // for GridOrdering enum
#include "pme-gpu-program-impl.h"
#include "pme-gpu-types-host.h"

namespace gmx
{

PmeGpuProgramImpl::PmeGpuProgramImpl(const gmx_device_info_t *deviceInfo)
{
    // Context creation (which should happen outside of this class: #2522
    cl_platform_id        platformId = deviceInfo->ocl_gpu_id.ocl_platform_id;
    cl_device_id          deviceId   = deviceInfo->ocl_gpu_id.ocl_device_id;
    cl_context_properties contextProperties[3];
    contextProperties[0] = CL_CONTEXT_PLATFORM;
    contextProperties[1] = (cl_context_properties) platformId;
    contextProperties[2] = 0; /* Terminates the list of properties */

    cl_int  clError;
    context = clCreateContext(contextProperties, 1, &deviceId, nullptr, nullptr, &clError);
    if (clError != CL_SUCCESS)
    {
        const std::string errorString = gmx::formatString("Failed to create context for PME on GPU #%s:\n OpenCL error %d: %s",
                                                          deviceInfo->device_name, clError, ocl_get_error_string(clError).c_str());
        GMX_THROW(gmx::InternalError(errorString));
    }

    warpSize = gmx::ocl::getWarpSize(context, deviceId);

    //TODO: OpenCL kernel compilation should be here.
}

PmeGpuProgramImpl::~PmeGpuProgramImpl()
{
    // TODO: log releasing errors
    clReleaseContext(context);
}

} //namespace gmx
