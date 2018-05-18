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
#ifndef GMX_EWALD_PME_PME_GPU_CONTEXT_H
#define GMX_EWALD_PME_PME_GPU_CONTEXT_H

#include <memory>

#include "gromacs/gpu_utils/gpu_macros.h"

struct PmeGpuContextImpl;
struct gmx_device_info_t;

class PmeGpuContext
{
    public:
        PmeGpuContext(const gmx_device_info_t *deviceInfo);
        const gmx_device_info_t *getDeviceInfo() const {return deviceInfo_; }
        PmeGpuContext() = default;

        // TODO: design getters for information inside, if needed for PME, and make this private?
        std::unique_ptr<PmeGpuContextImpl> impl_;
    private:
        const gmx_device_info_t           *deviceInfo_;
};

/*! \brief This is a handle for storing own PME GPU context data.
 * TODO: it should be unique_ptr, but careless copying of TestHardwareContext's in unit tests prevents it.
 */
using PmeGpuContextStorage = std::shared_ptr<PmeGpuContext>;

/*! \brief This is a handle for passing references to PME GPU context data.
 * TODO: it should be a const reference, but for that the PmeGpu types need to be C++
 */
using PmeGpuContextHandle = const PmeGpuContext *;

/*! \brief
 * Factory function used to build persistent data for the device at once.
 * \todo This should shortly become GPU_FUNC to support OpenCL.
 */
GPU_FUNC_QUALIFIER PmeGpuContextStorage buildPmeGpuContext(const gmx_device_info_t *) GPU_FUNC_TERM_WITH_RETURN(nullptr)

#endif
