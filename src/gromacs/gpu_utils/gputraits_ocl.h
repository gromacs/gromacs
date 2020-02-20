/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_GPUTRAITS_OCL_H
#define GMX_GPU_UTILS_GPUTRAITS_OCL_H

/*! \libinternal \file
 *  \brief Declares the OpenCL type traits.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_gpu_utils
 */

#include "gromacs/gpu_utils/gmxopencl.h"

//! OpenCL device vendors
enum class DeviceVendor : int
{
    Unknown = 0, //!< No data
    Nvidia  = 1, //!< NVIDIA
    Amd     = 2, //!< Advanced Micro Devices
    Intel   = 3, //!< Intel
    Count   = 4
};

/*! \internal
 * \brief OpenCL device information.
 *
 * The OpenCL device information is queried and set at detection and contains
 * both information about the device/hardware returned by the runtime as well
 * as additional data like support status.
 */
struct DeviceInformation
{
    cl_platform_id oclPlatformId;       //!< OpenCL Platform ID.
    cl_device_id   oclDeviceId;         //!< OpenCL Device ID.
    char           device_name[256];    //!< Device name.
    char           device_version[256]; //!< Device version.
    char           vendorName[256];     //!< Device vendor name.
    int            compute_units;       //!< Number of compute units.
    int            adress_bits;         //!< Number of address bits the device is capable of.
    int            stat;                //!< Device status takes values of e_gpu_detect_res_t.
    DeviceVendor   deviceVendor;        //!< Device vendor.
    size_t         maxWorkItemSizes[3]; //!< Workgroup size limits (CL_DEVICE_MAX_WORK_ITEM_SIZES).
    size_t maxWorkGroupSize; //!< Workgroup total size limit (CL_DEVICE_MAX_WORK_GROUP_SIZE).
};

//! \brief Single GPU call timing event
using CommandEvent = cl_event;

/*! \internal \brief
 * GPU kernels scheduling description. This is same in OpenCL/CUDA.
 * Provides reasonable defaults, one typically only needs to set the GPU stream
 * and non-1 work sizes.
 */
struct KernelLaunchConfig
{
    //! Work groups (CUDA blocks) counts
    size_t gridSize[3] = { 1, 1, 1 };
    //! Per work group (CUDA block) thread counts
    size_t blockSize[3] = { 1, 1, 1 };
    //! Shared memory size in bytes
    size_t sharedMemorySize = 0;
};

/*! \brief Sets whether device code can use arrays that are embedded in structs.
 * Note that OpenCL 2.x might be able to do this, but we use 1.2.
 */
#define c_canEmbedBuffers false

#endif
