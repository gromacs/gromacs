/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_DEVICE_CONTEXT_OCL_H
#define GMX_GPU_UTILS_DEVICE_CONTEXT_OCL_H

/*! \libinternal \file
 *
 * \brief Declarations for DeviceContext class.
 *
 * Only needed for OpenCL builds. Other platforms will be given a stub class.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 * \inlibraryapi
 */

#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/utility/classhelpers.h"

struct DeviceInformation;

// OpenCL device context class
class DeviceContext
{
public:
    /*! \brief Constructor that creates the \c cl_context
     *
     * \param[in] deviceInfo Platform-specific device information.
     *
     * \throws InternalError if context creation failed.
     */
    DeviceContext(const DeviceInformation& deviceInfo);
    //! Destructor
    ~DeviceContext();

    //! Get the associated device information
    const DeviceInformation& deviceInfo() const { return deviceInfo_; }
    //! Getter
    cl_context context() const;

private:
    //! A reference to the device information used upon context creation
    const DeviceInformation& deviceInfo_;
    //! OpenCL context object
    cl_context context_ = nullptr;

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(DeviceContext);
};

#endif // GMX_GPU_UTILS_DEVICE_CONTEXT_OCL_H
